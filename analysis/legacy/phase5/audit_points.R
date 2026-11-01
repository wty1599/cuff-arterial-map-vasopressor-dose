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

root <- '@CUFF_WORKSPACE@'
out <- file.path(root, 'phase5')
res <- file.path(out, '_restricted')
results <- list(); model_checks <- list()
labels <- c('<=0.1','>0.1-0.25','>0.25-0.5','>0.5-1.0','>1.0')
for (db in c('mimic','eicu')) {
 for(replicate in c(0,1)) {
  d <- read.csv(file.path(res, paste0('audit_original_pairs_', db, '.csv')))
  if(replicate==1) {
    copies <- read.csv(file.path(res,paste0('audit_first_cluster_weights_',db,'.csv')))$copies
    stopifnot(length(copies)==nrow(d), all(copies>=0))
    d <- d[rep(seq_len(nrow(d)),times=copies),]
  }
  d$outcome <- as.integer(d$cuff_map >= 65)
  d$low <- d$art_map < 65
  d$stratum <- cut(d$nee, breaks=c(-Inf,.1,.25,.5,1,Inf), labels=labels, right=TRUE)
  k <- quantile(d$art_map, probs=c(.05,.275,.5,.725,.95), type=7, names=FALSE)
  # Individual-observation binomial likelihood and R's natural-spline basis.
  # This does not import the production sufficient-count fit or spline basis.
  fit <- glm(outcome ~ splines::ns(art_map, knots=k[2:4], Boundary.knots=k[c(1,5)]),
             data=d, family=binomial(), control=glm.control(epsilon=1e-12,maxit=100))
  stopifnot(fit$converged, fit$rank==5L, all(is.finite(coef(fit))))
  spline <- as.numeric(predict(fit, type='response'))
  bins <- ifelse(d$art_map<40,'below40',ifelse(d$art_map>120,'above120',
                 ifelse(d$art_map==120,'exact120',paste0('interval_',floor(d$art_map)))))
  pbin <- tapply(d$outcome,bins,mean)
  empirical <- as.numeric(pbin[bins])
  model_checks[[paste(db,replicate)]] <- data.frame(database=db,replicate=replicate,converged=fit$converged,rank=fit$rank,
    iterations=fit$iter,knots=paste(format(k,trim=TRUE),collapse=';'),
    predicted_cuff_ge65_total=sum(spline),observed_cuff_ge65_total=sum(d$outcome),
    empirical_global_low_expected=sum(empirical[d$low]),
    observed_joint_low_normal=sum(d$outcome[d$low]),R_version=R.version.string)
  for (method in c('primary','empirical_1mmHg')) {
    p <- if(method=='primary') spline else empirical
    temp <- list()
    for (s in labels) {
      take <- d$stratum==s
      numerator <- sum(d$outcome[take & d$low]); expected_num <- sum(p[take & d$low])
      observed <- c(numerator/sum(take & d$low),numerator/sum(d$outcome[take]))
      expected <- c(expected_num/sum(take & d$low),expected_num/sum(p[take]))
      temp[[s]] <- list(observed=observed,expected=expected)
      if(method=='primary' || s=='>1.0') {
        for(j in 1:2) results[[length(results)+1]] <- data.frame(database=db,replicate=replicate,analysis=method,
          metric=c('missed_hypotension','art_lt65_given_cuff_ge65')[j],stratum=s,
          observed=observed[j],expected=expected[j],difference=observed[j]-expected[j],
          observed_increase=NA_real_,expected_increase=NA_real_,explained_ratio=NA_real_,
          low_observed=NA_real_,low_expected=NA_real_)
      }
    }
    od <- temp[['>1.0']]$observed[2]-temp[['<=0.1']]$observed[2]
    ed <- temp[['>1.0']]$expected[2]-temp[['<=0.1']]$expected[2]
    results[[length(results)+1]] <- data.frame(database=db,replicate=replicate,analysis=method,
      metric='explained_ratio',stratum='>1.0_vs_<=0.1',observed=NA_real_,expected=NA_real_,difference=NA_real_,
      observed_increase=od,expected_increase=ed,explained_ratio=if(abs(od)<1e-12) NA_real_ else ed/od,
      low_observed=temp[['<=0.1']]$observed[2],low_expected=temp[['<=0.1']]$expected[2])
  }
 }
}
write.csv(do.call(rbind,results),file.path(out,'phase5_independent_R_points.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,model_checks),file.path(out,'phase5_independent_R_models.csv'),row.names=FALSE,na='')
cat('INDEPENDENT_R_POINT_FITS_COMPLETE\n')
