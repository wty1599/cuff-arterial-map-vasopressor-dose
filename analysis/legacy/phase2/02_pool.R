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

options(warn=1)
suppressPackageStartupMessages({library(data.table);library(metafor)})
root<-'@CUFF_WORKSPACE@/phase2'
d<-fread(file.path(root,'phase2_A_deep_miss.csv'))
out<-list()
for(met in c('deep_miss_art_lt55','deep_miss_art_lt60'))for(s in c('known','<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0')){
  a<-d[d$metric==met & d$stratum==s]
  z<-a[denominator>0 & is.finite(logit_estimate)&is.finite(logit_se)&logit_se>0]
  row<-data.table(metric=met,stratum=s,databases=nrow(z),included=paste(z$database,collapse=';'),excluded=paste(setdiff(a$database,z$database),collapse=';'),numerator=sum(z$numerator),denominator=sum(z$denominator))
  if(nrow(z)>=2){m<-rma.uni(yi=z$logit_estimate,sei=z$logit_se,method='REML',test='adhoc');row[,`:=`(estimate=plogis(as.numeric(coef(m))),ci_lower=plogis(m$ci.lb),ci_upper=plogis(m$ci.ub),I2=m$I2,tau2=m$tau2,Q=m$QE,Q_p=m$QEp,status='ok')]}
  else row[,`:=`(estimate=NA_real_,ci_lower=NA_real_,ci_upper=NA_real_,I2=NA_real_,tau2=NA_real_,Q=NA_real_,Q_p=NA_real_,status='fewer_than_two_usable_databases')]
  row[,precision:=ifelse(ci_upper-ci_lower>.30,'insufficient',ifelse(is.finite(estimate),'adequate','not_estimable'))]
  out[[length(out)+1L]]<-row
}
fwrite(rbindlist(out,fill=TRUE),file.path(root,'phase2_A_deep_pooled.csv'),na='NA')
capture.output(sessionInfo(),file=file.path(root,'phase2_R_session.txt'))
cat('PHASE2_POOL_COMPLETE\n')
