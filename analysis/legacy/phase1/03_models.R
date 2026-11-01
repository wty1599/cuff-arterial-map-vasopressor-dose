# Private-workspace guard: run only after workspace preparation.
.workspace <- "@CUFF_WORKSPACE@"
if (grepl("@", .workspace, fixed=TRUE) || !grepl("^(/|[A-Za-z]:/)", .workspace)) stop("Prepare an absolute private workspace before running analysis templates.")
.workspace <- normalizePath(.workspace, winslash="/", mustWork=FALSE)
.script_arg <- grep("^--file=", commandArgs(), value=TRUE)
.candidate <- if (length(.script_arg)) dirname(normalizePath(sub("^--file=", "", .script_arg[1]), winslash="/", mustWork=FALSE)) else getwd()
repeat {
 if (file.exists(file.path(.candidate, "analysis", "prepare_workspace.py")) && dir.exists(file.path(.candidate, "outputs", "aggregate"))) {
  if (.workspace == .candidate || startsWith(.workspace, paste0(.candidate, "/"))) stop("Patient-data workspace must remain outside repository checkout.")
 }
 .parent <- dirname(.candidate)
 if (.parent == .candidate) break
 .candidate <- .parent
}
# End private-workspace guard.

options(stringsAsFactors=FALSE,warn=1)
suppressPackageStartupMessages({library(data.table);library(lme4);library(metafor);library(sandwich)})
root <- '@CUFF_WORKSPACE@/phase1'
write_out <- function(x,name) fwrite(x,file.path(root,name),na='NA')
set.seed(20261027)
rates <- fread(file.path(root,'phase1_all_rates.csv'))
levs <- c('known','<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0')
meta <- list()
for(ep in c('missed','false_positive')) for(lab in levs) {
  q <- rates[endpoint==ep & stratum==lab]
  z <- q[denominator>0 & is.finite(logit_estimate) & is.finite(logit_se) & logit_se>0]
  row <- data.table(endpoint=ep,stratum=lab,databases=nrow(z),included=paste(z$database,collapse=';'),excluded=paste(setdiff(q$database,z$database),collapse=';'),numerator=sum(z$numerator),denominator=sum(z$denominator),denominator_patients=sum(z$denominator_patients))
  if(nrow(z)>=2) {
    fit <- rma.uni(yi=z$logit_estimate,sei=z$logit_se,method='REML',test='adhoc')
    row[,`:=`(estimate=plogis(as.numeric(coef(fit))),ci_lower=plogis(fit$ci.lb),ci_upper=plogis(fit$ci.ub),I2=fit$I2,tau2=fit$tau2,Q=fit$QE,Q_p=fit$QEp,status='ok')]
  } else row[,`:=`(estimate=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,I2=NA_real_,tau2=NA_real_,Q=NA_real_,Q_p=NA_real_,status='fewer_than_two_databases')]
  row[,precision:=ifelse(ci_upper-ci_lower>.30,'insufficient',ifelse(is.finite(estimate),'adequate','not_estimable'))]
  meta[[length(meta)+1L]] <- row
}
write_out(rbindlist(meta,fill=TRUE),'phase1_pooled.csv')
cat('META_COMPLETE\n')

d0 <- fread(file.path(root,'_restricted/model_input.csv'),na.strings=c('','NA'))
vars <- c('nee','art_map','age','sex','bmi40','af_model')
fl <- list()
for(db in unique(d0$database)) {
 q <- d0[database==db]
 fl[[length(fl)+1L]] <- data.table(database=db,criterion='low_arterial_input',pairs=nrow(q),patients=uniqueN(q$patient))
 for(v in vars) fl[[length(fl)+1L]] <- data.table(database=db,criterion=paste0('missing_',v),pairs=sum(is.na(q[[v]])),patients=uniqueN(q$patient[is.na(q[[v]])]))
 cc <- q[complete.cases(q[,..vars])]
 fl[[length(fl)+1L]] <- data.table(database=db,criterion='complete_case_model',pairs=nrow(cc),patients=uniqueN(cc$patient))
}
write_out(rbindlist(fl),'phase1_model_flow.csv')
d <- as.data.frame(d0[complete.cases(d0[,..vars])]);d$patient<-factor(d$patient);d$database<-factor(d$database);d$sex<-factor(d$sex)
d$art_scaled <- (d$art_map-60)/10;d$age_scaled<-(d$age-65)/10
knots <- as.numeric(quantile(d$nee,c(.1,.5,.9),names=FALSE,type=7))
rcs <- function(x) ((pmax(x-knots[1],0)^3) -(pmax(x-knots[2],0)^3)*(knots[3]-knots[1])/(knots[3]-knots[2]) +(pmax(x-knots[3],0)^3)*(knots[2]-knots[1])/(knots[3]-knots[2]))/(knots[3]-knots[1])^2
dose_center<-median(d$nee);dose_scale<-sd(d$nee)
curve_center<-0;curve_scale<-1
spline_ok<-length(unique(knots))==3
if(spline_ok) { raw_curve<-rcs(d$nee);curve_center<-mean(raw_curve);curve_scale<-sd(raw_curve);spline_ok<-is.finite(curve_scale)&&curve_scale>0 }
d$dose_linear<-(d$nee-dose_center)/dose_scale
if(spline_ok) d$dose_curve<-(rcs(d$nee)-curve_center)/curve_scale
write_out(data.table(parameter=c('knot_10','knot_50','knot_90','dose_center','dose_scale','curve_center','curve_scale'),value=c(knots,dose_center,dose_scale,curve_center,curve_scale)),'phase1_model_parameters.csv')
versions<-if(spline_ok)c('spline_random_intercept','linear_random_intercept','linear_cluster_robust') else c('linear_random_intercept','linear_cluster_robust')
attempts<-list();fit<-NULL
for(version in versions) {
  term <- if(version=='spline_random_intercept') 'dose_linear + dose_curve' else 'dose_linear'
  ff <- as.formula(paste('missed ~',term,'+ art_scaled + age_scaled + sex + bmi40 + af_model + database'))
  mixed <- version!='linear_cluster_robust'
  warnings<-character()
  out<-tryCatch(withCallingHandlers({
    if(mixed) glmer(update(ff,.~.+(1|patient)),data=d,family=binomial(),nAGQ=1,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000)))
    else glm(ff,data=d,family=binomial())
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(out,'error')) { ok<-FALSE;msg<-conditionMessage(out);sing<-NA }
  else if(mixed) { msg<-paste(unlist(out@optinfo$conv$lme4$messages),collapse='; ');ok<-is.null(out@optinfo$conv$lme4$messages)&&all(out@optinfo$conv$opt==0);sing<-isSingular(out,tol=1e-4) }
  else { ok<-isTRUE(out$converged);msg<-if(ok)'' else 'glm_nonconvergence';sing<-NA }
  attempts[[length(attempts)+1L]]<-data.table(version=version,converged=ok,singular=sing,message=msg,warnings=paste(unique(warnings),collapse='; '))
  write_out(rbindlist(attempts,fill=TRUE),'phase1_model_attempts.csv')
  cat('MODEL_ATTEMPT',version,ok,'\n');flush.console()
  if(ok){fit<-out;break}
}
if(is.null(fit)) {write_out(data.table(status='not_completed',reason='all_requested_model_fallbacks_failed'),'phase1_model_status.csv');quit(status=0)}
saveRDS(fit,file.path(root,'_restricted/dose_association_model.rds'))
if(mixed){b<-fixef(fit);V<-as.matrix(vcov(fit));random_sd<-sqrt(as.numeric(VarCorr(fit)$patient[1]))}else{b<-coef(fit);V<-vcovCL(fit,cluster=d$patient,type='HC1');random_sd<-0}
se<-sqrt(diag(V))
write_out(data.table(term=names(b),estimate=as.numeric(b),standard_error=se,odds_ratio=exp(b),or_ci_lower=exp(b-1.96*se),or_ci_upper=exp(b+1.96*se),p_value=2*pnorm(-abs(b/se)),model=version),'phase1_secondary_model_coefficients.csv')
id<-which(names(b)%in%c('dose_linear','dose_curve'));wald<-as.numeric(t(b[id])%*%solve(V[id,id,drop=FALSE],b[id]))
write_out(data.table(test=c('joint_dose','nonlinearity'),df=c(length(id),if('dose_curve'%in%names(b))1 else NA),statistic=c(wald,if('dose_curve'%in%names(b))(b['dose_curve']/se['dose_curve'])^2 else NA),p_value=c(pchisq(wald,length(id),lower.tail=FALSE),if('dose_curve'%in%names(b))2*pnorm(-abs(b['dose_curve']/se['dose_curve'])) else NA)),'phase1_secondary_dose_tests.csv')
# Gauss-Hermite nodes for a standard normal distribution by a Jacobi eigendecomposition.
nq<-30;J<-matrix(0,nq,nq);J[cbind(1:(nq-1),2:nq)]<-sqrt(1:(nq-1));J<-J+t(J);ee<-eigen(J,symmetric=TRUE);nodes<-ee$values;qw<-ee$vectors[1,]^2
covarff <- delete.response(terms(ff));Xbase<-model.matrix(covarff,d)
estimate_at <- function(x) {
  X<-Xbase;X[,'dose_linear']<-(x-dose_center)/dose_scale
  if('dose_curve'%in%colnames(X)) X[,'dose_curve']<-(rcs(x)-curve_center)/curve_scale
  eta<-as.vector(X%*%b);p<-numeric(nrow(d));deriv<-p
  for(j in seq_len(nq)){pj<-plogis(eta+random_sd*nodes[j]);p<-p+qw[j]*pj;deriv<-deriv+qw[j]*pj*(1-pj)}
  mu<-mean(p);grad<-colMeans(X*deriv);vg<-as.numeric(t(grad)%*%V%*%grad);sel<-sqrt(max(vg,0))/(mu*(1-mu))
  list(row=data.table(nee=x,estimate=mu,ci_lower=plogis(qlogis(mu)-1.96*sel),ci_upper=plogis(qlogis(mu)+1.96*sel),model=version,ci_method='delta_fixed_effect_covariance_random_variance_plugin',random_intercept_sd=random_sd),gradient=grad)
}
points<-c(.1,.25,.5,1,2);pr<-lapply(points,estimate_at)
write_out(rbindlist(lapply(pr,'[[','row')),'phase1_secondary_model_predictions.csv')
grid<-sort(unique(c(seq(max(0,min(d$nee)),min(5,max(d$nee)),length.out=121),points)))
curve<-rbindlist(lapply(grid,function(x)estimate_at(x)$row));curve[,outside_observed_range:=nee<min(d$nee)|nee>max(d$nee)]
write_out(curve,'phase1_secondary_model_curve.csv')
contr<-list()
for(x in c(.25,.5,1,2)) {a<-estimate_at(x);z<-estimate_at(.1);g<-a$gradient-z$gradient;rd<-a$row$estimate-z$row$estimate;s<-sqrt(as.numeric(t(g)%*%V%*%g));contr[[length(contr)+1L]]<-data.table(dose=x,reference=.1,risk_difference=rd,ci_lower=rd-1.96*s,ci_upper=rd+1.96*s,p_value=2*pnorm(-abs(rd/s)))}
write_out(rbindlist(contr),'phase1_secondary_model_contrasts.csv')
bins<-cut(d$nee,breaks=c(-Inf,.1,.25,.5,1,Inf),labels=c('<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0'))
write_out(data.table(database=d$database,stratum=bins,patient=d$patient)[,.(pairs=.N,patients=uniqueN(patient)),by=.(database,stratum)],'phase1_model_dose_support.csv')
write_out(data.table(status='complete',model=version,pairs=nrow(d),patients=nlevels(d$patient),random_intercept_sd=random_sd,dose_min=min(d$nee),dose_max=max(d$nee),warnings=paste(unique(warnings),collapse='; ')),'phase1_model_status.csv')
capture.output(sessionInfo(),file=file.path(root,'phase1_R_session.txt'))
cat('MODELS_COMPLETE\n')
