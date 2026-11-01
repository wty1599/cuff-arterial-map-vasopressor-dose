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

suppressPackageStartupMessages({library(data.table);library(lme4)})
root<-'@CUFF_WORKSPACE@/phase1'
d<-fread(file.path(root,'_restricted/ba_validation_mimic.csv'))
d[,patient:=factor(patient)]
m<-lmer(delta~1+(1|patient),data=d,REML=TRUE,control=lmerControl(optimizer='bobyqa',optCtrl=list(rhoend=1e-10)))
r<-fread(file.path(root,'phase1_secondary_bland_altman.csv'))[database=='mimic'&stratification=='overall']
actual<-c(fixef(m)[1],as.numeric(VarCorr(m)$patient[1]),sigma(m)^2)
expected<-c(r$bias,r$patient_variance,r$residual_variance)
out<-data.table(parameter=c('bias','patient_variance','residual_variance'),aggregate_REML=expected,independent_lmer_REML=actual,absolute_difference=abs(expected-actual),passed=abs(expected-actual)<1e-3)
fwrite(out,file.path(root,'phase1_bland_altman_validation.csv'))
stopifnot(all(out$passed))
cat('INDEPENDENT_BA_VALIDATION_COMPLETE\n')
