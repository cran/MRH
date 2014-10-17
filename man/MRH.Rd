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
\donttest{fit.tongue = estimateMRH(formula = Surv(time, delta) ~ nph(type), data = tongue, M = 4, 
	burnIn = 200, maxIter = 2000, thin = 1, outfolder = 'MRH_tongue_nph')}

\donttest{MRHchains = read.table('MRH_tongue_nph/MCMCchains.txt', header = TRUE)
MRHobject = MRH(MRHchains)
summary(MRHobject)
DIC(MRHobject)}
}
