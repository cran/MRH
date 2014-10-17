\name{DIC}
\alias{DIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	DIC CALCULATES INFORMATION CRITERION GIVEN THE MCMC CHAINS FROM THE estimateMRH routine.
}
\description{
	The DIC function calculates the Deviance Information Criterion given the MCMC chains from an estimateMRH routine, using the formula: DIC = .5*var(D)+mean(D), where D is the chain of -2*log(L), calculated at each retained iteration of the MCMC routine.  It also provides the Aikaike Information Criterion (AIC) = 2*k + D and the Bayesian Information Criterion (BIC) = k*ln(n) + D, where k is the number of parameters in the model.  Both AIC and BIC report the maximum (i.e. "worst") values in the chain.
}
\usage{
DIC(mrhobject, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mrhobject}{ The chains found in the MCMCchains.txt file, created using the estimateMRH routine, or the MRH results object.}
  \item{n}{The sample size of the original dataset. If n is not entered, the BIC calculation will not be returned to the user.}
}
\references{
Gelman, A., Carlin, J., Stern, H., Dunson, D., Vehtari, A., Rubin, D. (2004) \emph{Bayesian Data Analysis}. 

Spiegelhalter, D.J., Best, N.G., Carlin, B.P., van der Linde, A. (2002), Bayesian measures of model complexity and fit (with discussion). \emph{Journal of the Royal Statistical Society, Series B}. \bold{64, } 583--639.

}
\value{ DIC returns the DIC, AIC, and BIC values, as well as a summary of D (-2*loglike).
	\item{loglik.summ}{The summary of the chain of -2*loglike values.}
	\item{ICtable}{Table containing the DIC, AIC and BIC values.}
}
\author{
	Yolanda Hagar <yolanda.hagar@colorado.edu> and Vanja Dukic
}
\examples{
# These MRH fit portion of the examples 
# are from the estimateMRH() help page.  
# They do not need to be re-run if the 
# objects are already in the active workspace.
data(cancer)
cancer$censorvar = cancer$status - 1
\donttest{fit.lung = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung')

fit.lung.prune = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, prune = TRUE, outfolder = 'MRH_lung_prune')}
	
# Compare the DIC of the pruned and unpruned models:
\donttest{DIC(fit.lung)
DIC(fit.lung.prune)}

}
