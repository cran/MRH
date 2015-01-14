\name{plot.MRH}
\alias{plot.MRH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plots the hazard rate (with credible bands) of an MRH object.
}
\description{
	plot.mrh plots the hazard rate of an MRH object with the 95\% credible interval bounds included as dashed lines.
}
\usage{
\method{plot}{MRH}(x, maxStudyTime, main = "", xlab = "Time", 
ylab = "Hazard Rate", plot.type = 'h', 
interval = TRUE, alpha.level = 0.05, smooth.graph = FALSE, 
smooth.df = NULL, combine.graphs = TRUE, log.ratio = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	x is an MRHobject: Either 1) an MRH fitted object or 2) MCMC chains produced by the estimateMRH routine, converted to an MRH object using MRH() or as.MRH().
}
  \item{maxStudyTime}{
	The maximum observed or censored failure time needs to be entered when the text file of MCMC chains is used instead of the active MRH fitted object for calculation of the hazard rate.
}
  \item{main}{
	The main title of the graph.
}
  \item{xlab}{
	The label of the x-axis.
}
  \item{ylab}{
	The label of the y-axis.
}
\item{plot.type}{Denotes whether the hazard rate, cumulative hazard, survival function or hazard ratio (for the case of non-proportional hazards) should be plotted, specified with 'h', 'H', 'S', or 'r' respectively.  Default plot is the hazard rate.}
\item{interval}{Set to TRUE or FALSE if the credible interval bounds should also be included in the plot.  Default is TRUE.}
\item{alpha.level}{1 minus the credible level of the interval.  Default is 0.05.}
  \item{smooth.graph}{
	Should be set equal to TRUE if the user would like a graph of the smoothed hazard rate using smooth.spline().  Default is set to FALSE.
}
  \item{smooth.df}{
	The degrees of freedom to be used in the smoothed graph of the hazard rate.  If the user does not enter a value, it is set to 1/2 of the number of bins.
}
  \item{combine.graphs}{
	If the MCMC chains are from a non-proportional hazards model, combine.graphs can be set to FALSE if the user would like 
	to view the hazard rates or ratios separately.  The default is TRUE, combining the hazard rates or ratios in to one graph. 
}
\item{log.ratio}{In the non-proportional hazards models, this denotes whether the log-ratio of the hazard rates should be plotted over time (log.ratio = TRUE) or if the ratio of the hazard rates should be plotted instead (log.ratio = FALSE).  The default is a plot of the log-ratio of the hazard rates.}
\item{...}{
	Arguments to be passed to methods, such as graphical parameters (see par).
}
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
	
# Plot the hazard rate, cumulative hazard, 
# and survival function of the lung model side-by-side:
\dontrun{
par(mfrow = c(1,3))
plot(fit.lung, main = 'Hazard rate')
plot(fit.lung, plot.type = 'H', main = 'Cumulative hazard')
plot(fit.lung, plot.type = 'S', main = 'Survival function')}

# Plot the hazard rates for the tongue (i.e. non-proportional hazards) 
# model with smoothed estimates.
\dontrun{
plot(fit.tongue, smooth.graph = TRUE)}
# Separate the graphs for less crowding
\dontrun{
plot(fit.tongue, smooth.graph = TRUE, combine.graphs = FALSE)}
# Plot the hazard ratios
\dontrun{
plot(fit.tongue, plot.type = 'r')}

# Plot the hazard rate of the lung model using the chains. 
# This requires maxStudyTime to be entered.
\dontrun{
lung.chains = MRH(read.table('MRH_lung/MCMCchains.txt', header = TRUE))
plot(lung.chains, maxStudyTime = 960)}
}
