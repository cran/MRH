\name{summary.MRH}
\alias{summary.MRH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Summarizes the results of an MRH object
}
\description{
	This function summarizes the MRH object (either the fitted object or the MCMC chains) and returns estimates of the covariate effects, the hazard rate, the survival curve, the cumulative hazard, and the log ratio for non-propoortional hazards models.  If desired, the user can specify the alpha-level for the credible intervals.
}
\usage{
\method{summary}{MRH}(object, alpha.level, maxStudyTime, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	object is an MRHobject: Either 1) an MRH fitted object or 2) MCMC chains produced by the estimateMRH routine, converted to an MRH object using MRH() or as.MRH().
}
  \item{alpha.level}{
	The width of the credible intervals for the parameters, with the lower bound calculated as the alpha/2 percentile, and the upper bound calculated as the 1-alpha/2 percentile.
}
  \item{maxStudyTime}{
	The maximum study period (or the censoring time) used in the analysis. This is only required if the MCMC chains are used instead of the fitted MRH model.
}
\item{...}{
	Arguments to be passed to methods, such as graphical parameters (see par).
}
}
\value{
	\item{d}{The hazard increments within each bin (i.e. the cumulative hazard within each bin)}
	\item{hazardRate}{The hazard rate across the bins}
	\item{SurvivalCurve}{The survival curve (i.e. S(t)) over the course of the study.}
	\item{CumulativeHazard}{The cumulative hazard (i.e. H(t)) over the course of the study.}
	\item{beta}{The estimated covariate effects.  If the model contains both PH and NPH parameters, the covariate effects will be denoted with 'betaPH' and 'betaNPH'.}
	\item{H}{The parameter H, which is the cumulative hazard at the end of the study.}
	\item{Rmp}{The split parameters.} 
	\item{maxStudyTime}{The maximum study time (i.e. the censoring time).}
}
\author{
	Yolanda Hagar <yolanda.hagar@colorado.edu> and Vanja Dukic
}
\examples{
# These MRH fit portion of the examples are from the 
# estimateMRH() help page.  
# They do not need to be re-run if the objects 
# are already in the active workspace.
data(cancer)
cancer$censorvar = cancer$status - 1
\dontrun{
fit.lung = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung')}

data(tongue)
\dontrun{
fit.tongue = estimateMRH(formula = Surv(time, delta) ~ 
	nph(type), data = tongue, M = 4, 
	burnIn = 200, maxIter = 2000, thin = 1, outfolder = 'MRH_tongue_nph')}
	
# Summarize the models
\dontrun{
summary(fit.lung)
summary(fit.lung, alpha.level = .01)

mcmc.lung = as.MRH(read.table('MRH_lung/MCMCchains.txt', header = TRUE))
summary(mcmc.lung, maxStudyTime = 960)
summary(mcmc.lung, maxStudyTime = 960, alpha.level = .01)
}
}
