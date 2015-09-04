\name{CalcFunction}
\alias{CalcFunction}
\title{CALCFUNCTION CALCULATES THE HAZARD RATE, CUMULATIVE HAZARD AND/OR SURVIVAL FUNCTION OF AN MRH OBJECT
	}
\description{
This function calculates the hazard rate, cumulative hazard and/or the survival function of an MRH object.  The alpha-level bounds are included in the calculation.}
\usage{
CalcFunction(mrhobject, function.type = c("h", "H", "S"), 
maxStudyTime, alpha.level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mrhobject}{
%%     ~~Describe \code{mrhobject} here~~
An MRH object: either the fitted model from the estimateMRH() routine or the text file of MCMC chains converted to an MRH object.
}
  \item{function.type}{
  	The function the user would like returned, with 'h' denoting the hazard rate, 'H' denoting the cumulative hazard, and 'S' denoting the survival function.  The user may request any subset of these functions, and multiple functions can be requested and returned.
}
  \item{maxStudyTime}{
	The maximum observed or censored failure time, or the end of the study period. This is only needed if the text file of the chains are used instead of the actived fitted MRH object.
}
  \item{alpha.level}{
	The alpha.level for the bounds of the credible intervals.
}
}
\value{
\item{hazardrate}{The hazard rate and alpha-level bounds. Returned to the user if 'h' is entered in  function.type.}
\item{cumulhaz}{The cumulative hazard functions and alpha-level bounds.  Returned to the user if 'H' is entered in function.type.}
\item{survfunction}{The survival function and alpha-level bounds.  Returned to the user if 'S' is entered in function.type}
}
\author{
	Yolanda Hagar <yolanda.hagar@colorado.edu> 
}
\examples{
# The MRH fit portion of the example is from 
# the estimateMRH() help page.  They do not need to 
# be re-run if the code has previously 
# been run and the outfolder ('MRH_lung') is saved.
data(cancer)
cancer$censorvar = cancer$status - 1
\dontrun{
fit.lung = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, 
	maxIter = 1000, thin = 1, outfolder = 'MRH_lung')}

# Get the cumulative hazard with 99% credible interval bounds
\dontrun{
CalcFunction(fit.lung, function.type = 'H', 
	alpha.level = .01)}

# Get the hazard rate, cumulative hazard, and 
# survival function using the MCMC chains text file
\dontrun{
lung.chains = MRH(read.table('MRH_lung/MCMCchains.txt', header = TRUE))}

# When the chains are read in from a file and used, 
# the maximum study time 
# must be entered to calculate the hazard rate.
# (i.e. CalcFunction(lung.chains) does not work)
\dontrun{
CalcFunction(lung.chains, maxStudyTime = 960)}
}
