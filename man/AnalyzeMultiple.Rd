\name{AnalyzeMultiple}
\alias{AnalyzeMultiple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
%%  ~~function to do ... ~~
ANALYZE MULTIPLE MCMC CHAINS
}
\description{
Analyzes multiple MCMC chains and returns the median and alpha-level quantiles of the marginal posterior distribution for each parameter. The Gelman-Rubin test can be employed to check for convergence if the initial values for the parameter chains vary across the sample space.  
}
\usage{
AnalyzeMultiple(datalist, fileNames, alpha.level, maxStudyTime, GR = TRUE)
}
\arguments{
  \item{datalist}{
A list object that contains one set of MCMC chains per list element.  
}
  \item{fileNames}{
If datalist is empty, a list of filenames entered as a string of character names can be entered here.
}
  \item{alpha.level}{
The alpha-level bounds for the credible intervals for the parameter estimates.
}
  \item{maxStudyTime}{
The maximum observed or censored failure time, or the end of the study period.  
}
  \item{GR}{
A TRUE or FALSE value denoting whether the Gelman-Rubin diagnostic convergence test should be performed on the chains.  Default is TRUE, but can be set to FALSE if the user only wishes to obtain estimates across multiple data sets. If the Gelman-Rubin test is employed, the chains must have the same burn-in, thinning value, and maximum number of iterations, and the initialized parameter values must cover the parameter space.  The most efficient way to ensure this is to set GR = TRUE when fitting the MRH object using the estimateMRH() function.
}
}
\details{
If the Gelman-Rubin test is used, the MCMC chains must satisfy the assumptions needed for the test and the initial parameter values used to estimate the posteriors must be dispersed across the parameter space.  In addition, the burn-in number, thinning value, and maximum number of iterations must be equal across all chains.  This can be enforced in estimateMRH() by setting GR equal to TRUE.  Because setting GR equal to TRUE fixes the burn-in number, thinning value, and maximum number of iterations, it is critical that the user enter values that are sufficient for convergence, while optimizing the run time for the fitted MRH model.
}
\value{
\item{hazardRate}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for the hazard rate within each bin.}
\item{beta}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for each covariate.}
\item{SurvivalCurve}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for the survival curve.}
\item{CumulativeHazard}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for the cumulative hazard over the course of the study.}
\item{d}{The estimated cumulative hazard rate increments for each of the $2^M$ bins, with the medians and (alpha/2, 1-alpha/2) quantiles of the posterior distribution for each d reported.}
\item{H}{The median and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for the cumulative hazard H at the end of the study.}
\item{Rmp}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for each of the split parameters (Rmp).}
\item{gamma}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for each of the gamma parameters (used in the posterior of the Rmp split parameters).  This is only returned to the user if the gamma parameters were sampled in the original MCMC sampling routine.}
\item{k}{The medians and alpha/2\% and 1-alpha/2\% quantiles of the posterior distribution for k (used in the posterior of the Rmp split parameters).  This is only returned to the user if k was sampled in the original MCMC sampling routine.}
\item{gelman.rubin}{If option is set to TRUE, the results of the Gelman-Rubin test.  The chains have converged of the scale reduction factor is "close enough" to 1.}
}
\author{
	Yolanda Hagar
}
\examples{
# Generate 3 chains for the same model.  
# Set GR = TRUE so that the burn-in number, 
# thinning value, 
# and maximum number of iterations are the same 
# across all MCMC chains, and so that initialized 
# parameter values cover the parameter space.

# Note that the routine may produce a warning 
# message that the algorithm has not converged,
# as typically more iterations are needed for convergence.
# However, for the purposes of this example, the number
# of iterations is sufficient.

data(cancer)
cancer$censorvar = cancer$status - 1
\dontrun{
fit.lung1 = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung1', GR = TRUE)
fit.lung2 = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung2', GR = TRUE)
fit.lung3 = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung2', GR = TRUE)}

# Calculate the results of all three chains and the 
# scale reduction factor entering the data sets in a list
\dontrun{
AnalyzeMultiple(datalist = 
	list(read.table('MRH_lung1/MCMCchains.txt', header = TRUE), 
	read.table('MRH_lung2/MCMCchains.txt', header = TRUE),
	read.table('MRH_lung2/MCMCchains.txt', header = TRUE)), 
	maxStudyTime = 960)}

# Calculate the results of all three chains and the 
# scale reduction factor entering the data file names
\dontrun{
AnalyzeMultiple(fileNames = c('MRH_lung1/MCMCchains.txt', 
				'MRH_lung2/MCMCchains.txt', 'MRH_lung2/MCMCchains.txt'), 
				maxStudyTime = 960)}
}
