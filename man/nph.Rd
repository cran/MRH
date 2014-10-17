\name{nph}
\alias{nph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	INTERNAL FUNCTION NEEDED FOR NON-PROPORTIONAL HAZARDS MODELS.
}
\description{
	nph is used in the MRH survival model formula to denote that a variable should be modeled under the non-proportional hazards assumption.
}
\usage{
nph(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	x is the covariate that will be modeled under the non-proportional hazards assumption.  It must be a categorical variable, although it need not be entered using the "factor" class in R.  A separate hazard rate will be estimated for each group in x.
}
}
\references{
Dukic, V., Dignam, J. (2007), Bayesian hierarchical multiresolution hazard model for the study of time-dependent failure patterns in early stage breast cancer. \emph{Bayesian Analysis}. \bold{2,} 591--610.
}
\author{
	Yolanda Hagar <yolanda.hagar@colorado.edu> and Vanja Dukic
}
\examples{
# Examine the tongue data set (from the KMest package), and
# quantify how the rumor DNA profile (1=Aneuploid Tumor, 2=Diploid Tumor)
# affects survival time (in weeks).
data(tongue)

# Fit the MRH model, including tumor type using the 
# non-proportional hazards model.  
# With 16 bins (M = 4), each bin represents 25 weeks.
# Generally it is recommended to use a higher burn-in value, 
# thinning value, and maximum number
# of iterations, but for illustrative purposes these values 
# have been lowered.
\donttest{fit.tongue = estimateMRH(formula = Surv(time, delta) ~ 
	nph(type), data = tongue, M = 4, 
	burnIn = 200, maxIter = 2000, thin = 1, 
	outfolder = 'MRH_tongue_nph')}
}
