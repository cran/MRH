\name{MRH}
\alias{MRH}
\alias{MRH.default}
\alias{as.MRH}
\alias{is.MRH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Converts an MRH MCMC chains text file into an MRH object
}
\description{
	NOTE: This function converts the MRH MCMC chains file to an MRH object.  For help with the MRH package, please see \code{\link{MRH-package}} for more information.
	
	When the MRH model is run as a background routine, the MRH object is not available in the workspace.  The MRH() and as.MRH() functions convert the MCMC chains text file produced by the estimateMRH() routine in to an MRH object. This then allows the user to use the summary and plotting functions on the MRH object. 
}
\usage{
MRH(x)
as.MRH(x)
is.MRH(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	The MCMC chains text file produced by the MRH estimation routine.
}
}
\author{
	Yolanda Hagar <yolanda.hagar@colorado.edu>
}
\examples{
# The MRH fit example is from the estimateMRH() help page.  It does not need to 
# be re-run if the code has previously been run and the outfolder ('MRH_tongue_nph') is
# saved.
data(tongue)
\dontrun{
fit.tongue = estimateMRH(formula = Surv(time, delta) ~ nph(type), data = tongue, M = 4, 
	burnIn = 200, maxIter = 2000, thin = 1, outfolder = 'MRH_tongue_nph')}

\dontrun{
MRHchains = read.table('MRH_tongue_nph/MCMCchains.txt', header = TRUE)
MRHobject = MRH(MRHchains)
# When the chains are read in, the maximum study time must
# be entered in the summary function, and the total number
# of subjects (n) must be entered for BIC to be calculated.
summary(MRHobject, maxStudyTime = max(tongue$time))
DIC(MRHobject, n = nrow(tongue))}
}
