\name{estimateMRH}
\alias{estimateMRH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Joint estimation of the hazard rate and covariate effects using multi-resolution hazard methodology}
\description{
estimateMRH is a function used to jointly estimate the hazard rate for time-to-event data using the multi-resolution hazard method.  Covariates can be included with or without the proportional hazards assumption.
}
\usage{
estimateMRH(formula, data, M, numberBins, maxStudyTime, outfolder = "MRHresults", 
	prune = FALSE, prune.alpha = 0.05, prune.levels = NULL, burnIn = 50000, 
	maxIter = 5e+05, thin = 10, Rmp.init = NULL, a.init = 10, lambda.init, 
	beta.init = NULL, k.fixed, gamma.fixed, GR = FALSE, convGraphs = TRUE, 
	fix.burnIn = FALSE, fix.thin = FALSE, fix.max = FALSE, continue.chain = FALSE)
}
\arguments{
  \item{formula}{
A formula object, with the response on the left of a ~ operator, and the covariates on the right. The response must be a survival object as returned by the Surv function.  The non-proportional covariates can be included in the formula by using the function nph(). (See Example.)
  }
  \item{data}{A data frame in which to interpret the variables named in the formula.}
  \item{M}{
	The value of M that dictates the time "resolution" of the hazard estimate, approximated via  2^M hazard increments.
}
  \item{numberBins}{
	An optional value that can be used instead of 'M' to dictate the number of bins.  This value must be a power of 2.
}
\item{maxStudyTime}{
	Allows the user to set the end of the study time to something other than the maximum observed or censored failure time.  This is typically used if a specific bin width is desired.
}
  \item{outfolder}{
The name given to the folder (subdirectory of the working directory) that stores the output for the MRH estimation routine (output includes the MCMC chains text file, the convergence graphs, and hazard rate graphs).  Note that this folder will be automatically placed in the working directory.  A pathname for the folder may also be specified, but the path must also be accessible from the working directory.
}
  \item{prune}{
  	If set to TRUE, two adjacent bins that are constructed via the same split parameter (Rm,p, or 1-Rm,p) are merged if the estimated hazard rate in these two bins are statistically similar. Pruning hypothesis tests (with the null hypothesis that the two bins have equal hazard rates) are performed using a modified Fisher's exact test, based on the 2 by 2 table composed of the number of failures within the time interval and at-risk patients at the end of the  time interval, for each pair of adjacent bins sharing a split parameter. If the test fails to reject the null hypothesis at the alpha level, that split parameter (Rm,p) is set to 0.5.  If a non-proportional covariate is used, each hazard rate is pruned separately.  The user may also enter a pruning vector of their choice, which must contain a 0 or 1 for each Rmp in the model. For the non-proportional hazards model, the user-entered pruning indicator can either be a vector, which prunes all hazard rates the same, or it can be a matrix with each column containing a pruning vector for each hazard rate group, allowing the pruning for each hazard rate to be different. To make a user-defined pruning vector, please email Vanja Dukic <vanja.dukic@colorado.edu> or Yolanda Hagar <yolanda.hagar@colorado.edu>.  
}
  \item{prune.alpha}{
	The significance level for Fisher's exact test if the "prune" option is set to TRUE.  The default is alpha = 0.05.
}
  \item{prune.levels}{
	If pruning is used, the number of desired levels to be pruned (ranges from 1 to M).  The default to prune all M levels.  For multiple hazards, a single value or a vector with different pruning levels may be entered, with one value for each hazard.
}
  \item{burnIn}{
The number of iterations in the burn-in for the MCMC chain.  Default is 50,000.  See details.
}
  \item{maxIter}{
%%     ~~Describe \code{maxIter} here~~
The maximum number of iterations in the MCMC routine.  Default is 1,000,000. 
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
The thinning parameter that denotes the thinning of the MCMC chain to reduce autocorrelation.  Default value is 10.  See Details.   
}
  \item{Rmp.init}{
%%     ~~Describe \code{RmpInit} here~~
The initial  values for Rmp (split) parameters.  If no values are entered, the default initial values is 0.5 for all Rmp.
}
  \item{a.init}{
%%     ~~Describe \code{a.init} here~~
The initial value for "a", a parameter in the (Gamma) prior for the baseline cumulative hazard H.  If no value is entered, the default value is 10.

}
  \item{lambda.init}{
%%     ~~Describe \code{lambda.init} here~~
The initial value for "lambda", a parameter in the (Gamma) prior for the baseline cumulative hazard H.  If no value is entered, the default value is -log(mean(delta))/a, where delta is the vector of censoring indicators.  This initialization of lambda approximates the cumulative hazard rate divided by the value of a.
}
  \item{beta.init}{
%%     ~~Describe \code{beta.init} here~~
The initial value(s) for the beta (covariate effect) parameter(s).  If no value is entered, the default values are obtained from the Cox proportional hazards estimates.
}
\item{k.fixed}{
	Can either be a boolean indicator or a numeric value assigned to k, a parameter in the Beta prior for the split parameters Rmp. By default (or if set equal to TRUE), k is fixed at 0.5, implying zero prior correlation among the hazard increments.  The user may also specify the fixed value of k by entering a numeric value for k.fixed greater than 0.  When k is greater than 0.5, the increments are positively correlated a priori.  Correspondingly, for k less than 0.5, the hazard increments will have negative prior correlation. If k.fixed is set to FALSE, then k will be sampled in the MCMC routine.
}
\item{gamma.fixed}{ 
	Can either be a boolean indicator or a numeric vector of length 2^M-1 for values assigned to the gamma parameters (parameters in the Beta prior for the split parameters Rmp), with one gamma.mp associated with each Rmp, and E(Rmp) = gamma.mp.  The gamma parameter allows the user to a priori "center" the baseline hazard in each bin. If gamma.fixed is omitted or set equal to TRUE, the values are fixed at 0.5, assuming centering occurs in the middle of the bin.  The user may also specify values (between 0 and 1) for each gamma.  If gamma.fixed is set to FALSE, then each gamma parameter will be sampled in the MCMC routine.
}	
\item{GR}{
%%     ~~Describe \code{GR} here~~
Denotes whether the Gelman-Rubin test statistic will be applied to the results of the MCMC chain.  See Details.
}
  \item{convGraphs}{
%%     ~~Describe \code{convGraphs} here~~
Allows user to specify if graphs for convergence should be saved.  Includes trace plots, density plots, and running means.  Default is TRUE.
}
  \item{fix.burnIn}{
%%     ~~Describe \code{fix.burnIn} here~~
If set equal to TRUE, the specified burn-in number (i.e. the ""burnIn"" value) will be fixed throughout the routine, and will not be increased, regardless of the convergence diagnostics performed by the routine.  See details.
}
  \item{fix.thin}{
%%     ~~Describe \code{fix.thin} here~~
If set equal to TRUE, the specified thinning value (i.e. the "thin" value) will be fixed throughout the routine, and will not be increased, regardless of the autocorrelation diagnostics performed by the routine.  See details.
}
  \item{fix.max}{
%%     ~~Describe \code{fix.max} here~~
If set equal to TRUE, the specified maximum number of iterations (i.e. the "maxIter" value) will be fixed throughout the routine, and will not be decreased, regardless of convergence diagnostics.  See details.
}
\item{continue.chain}{
	If set to TRUE, the MCMC routine will continue running where the routine left off.  The same output folder (``outfolder"") must be used, and parameters will be initialized with the last recorded sample from the MCMC chain in the previous run.  Any thinning or burn-in specifications made in the model call will be ignored, and only values from the text file containing the chains will be used.  The maximum number of new iterations is specified by ``maxIter", and the new MCMC chains are appended to the existing text file in the output folder.
}
}
\details{
This function returns the estimate of the hazard rate using the Multi-Resolution Hazard (MRH) method.  The user must have the survival time and censoring variable for each subject.  Parameters are estimated using MCMC.

After the first 100,000 MCMC iterations, the results are checked to determine an appropriate thinning value.  The auto-correlation estimate, available via the coda package, is examined, and the first lag with an autocorrelation below 0.10 for any parameter is taken as the new thinning value.  Default or user entered burn-in and thin values are taken as the minimum possible thinning values, and the routine only checks if the values should be greater (not smaller) than those specified by the user.

Every 100,000 iterations, the results are checked to determine if there is evidence of convergence and to determine an appropriate burn-in number.  This is done through the Geweke diagnostic test using the geweke.diag() function available in the coda package.  The geweke.diag() function returns the resulting z-score from the Geweke test for each parameter.  If any of the z-scores are outside of the 0.005 range, then convergence is assumed not to have been reached.  If convergence is not reached, 20,000 more iterations are burned, and the Geweke diagnostic test is performed again.  This continues until there are fewer than 1,000 retained iterations left (i.e. retained after thinning).  If there is no evidence of convergence, and if there are fewer than 1,000 retained iterations, the MCMC sampler runs another 100,000 times.  This continues until the maximum number of iterations (maxIter) is reached, at which point the routine stops and the results are returned to the user.  If the maximum number of iterations is reached before convergence and the user would like to continue the routine, the last parameter estimates in the output file can be used as initial values in the another call of the routine.  The user has the option to fix the burn-in number, maximum number of iterations, and thinning value if fix.burnIn, fix.max, or fix.thin are set equal to TRUE.  The convergence routine will not change or check these values if they are fixed.

If the user would like to run the routine multiple times on the same data set and test for convergence using the Gelman-Rubin diagnostic test,  the user should specify GR = TRUE when calling the MRH estimation function.  This way the initial parameter values will be randomly sampled to explore the parameter space.  Additionally, this will fix the user-entered or default burn-in number, and thinning value, and the number of MCMC iterations will be fixed at the maximum number of iterations (maxIter) so that the diagnostic test can be performed.  NOTE: If the user desires less than the default 1,000,000 maximum number of iterations and has set GR = TRUE, maxIter should be changed.  Likewise, if fix.max is set equal to TRUE, the routine will run exactly maxIter number of times, regardless of convergence.

There may be instances in which the user may desire to combine some of the hazard rate bins, particularly if there is evidence that the hazard rate does not change from one bin to the next.  In these cases, a pruning indicator can be used to denote which bins can be combined, with a "0" indicating the bin should not be combined with another bin, and a "1" indicating it should.  While it is possible to specify this manually, the user can  use the built-in pruning function by setting prune = TRUE, which performs Fisher's exact test to determine if bins can be combined based on observed failures and censoring in each bin.  Based on the results of the hypothesis tests, certain bins may be pruned. If there are multiple hazards being calculated (i.e. the non-proportional hazards model is used), by default each hazard rate is pruned separately.  The user may create their own pruning vector or matrix, although care must be taken in accounting for the tree-like structure of the Rmp parameters.  For assistance with creation of a pruning vector or matrix, please contact Vanja Dukic <vanja.dukic@colorado.edu> or Yolanda Hagar <yolanda.hagar@colorado.edu>.
}
\value{
estimateMRH returns a list with class MRH that contains the results of the MCMC algorithm. The routine also writes a file of the thinned and burned chain of MCMC iterations and saves pdf graphs of the hazard function with credible bounds in an MRH results folder (with the default title "MRHresults").  The components returned in the MRH fitted object are:
\item{summary}{A summary table containing the rounded estimates and central credible intervals for the estimates of the within-bin hazard rate (labeled ``h.binj" for the jth bin) and the covariate effects.  Estimates are calculated as the median of the thinned and burned MCMC chain, and bounds on the credible intervals are calculated as the 2.5\% and 97.5\% of the thinned and burned MCMC chain.}
\item{hazardRate}{The estimate and lower and upper bounds of the 95\% central credible interval for the hazard rate in each bin (labeled ``h.binj" for the hazard rate of bin j).}
\item{beta}{The estimate and lower and upper bounds of the 95\% central credible interval for the covariate effects.  This can include the log-ratios of the non-proportional hazards if the non-proportional assumption is used, which are labeled with ``binj" following the covariate name to denote the log-ratio estimate within the jth bin.}
\item{SurvivalCurve}{The estimate and lower and upper bounds of the 95\% central credible interval for the survival curve. In the non-proportional hazards setting, separate survival curves and credible intervals are provided for each strata of the covariate.}
\item{CumulativeHazard}{The estimate and lower and upper bounds of the 95\% central credible interval for the cumulative hazard. In the non-proportional hazards setting, separate cumulative hazard curves and credible intervals are provided for each strata of the covariate.}
\item{d}{The estimate and lower and upper bounds of the 95\% central credible interval for the within-bin cumulative hazard (denoted as `d') for each bin.}
\item{H}{The estimate and lower and upper bounds of the 95\% central credible interval for the baseline cumulative hazard H at the end of the study (denoted as `H00' or 'H(tJ)' in manuscripts).}
\item{Rmp}{The estimate and lower and upper bounds of the 95\% central credible interval for the split parameters Rmp.}
\item{gamma}{The estimate and lower and upper bounds of the 95\% central credible interval for the gamma parameter (used in the prior for the Rmp split parameters). This is only returned to the user if k.fixed is set equal to FALSE when fitting the model.}
\item{k}{The estimate and lower and upper bounds of the 95\% central credible interval for the k parameters (used in the priors for the Rmp split parameters). This is only returned to the user if gamma.fixed is set equal to FALSE when fitting the model.}
\item{AIC}{Akaike's Information Criterion, calculated as 2p-2ln(L), where k is the number of parameters in the model and ln(L) is the minimum of the likelihood values calculated at each chain iteration.  The number of parameters `p' is calculated as 2^M (one for each split parameter Rmp, and one for the cumulative hazard at H), plus 2 for a and lambda (parameters in the Gamma prior for H), and one for each covariate included under the proportional hazards assumption.  If k and/or gamma (parameters in the prior for Rmp) are sampled, the number of estimated parameters is increased by 1 for k and 2^M-1 for gamma.  If a covariate is included under the non-proportional hazards assumption, the number of estimated parameters (excluding any covariates included under the proportional hazards assumption) is multiplied by the number of strata in the non-proportional covariate.}
\item{BIC}{Bayesian Information Criterion, calculated as -2ln(L) + p*log(n)}
\item{DIC}{Deviance Information Criterion, calculated as .5*var(-2*ln(L)) + mean(-2ln(L)), where ln(L) is calculated at each retained chain iteration.}
\item{burnIn}{The number of iterations burned during the MCMC routine.}
\item{thin}{The thinning parameter used to account for autocorrelation.}
\item{TotalIters}{The number of iterations the MCMC algorithm performed before convergence.}
\item{convergence}{A TRUE/FALSE indicator denoting whether there is evidence that the algorithm converged based on methods used in the MRH routine.}
\item{gelman.rubin.used}{Denotes whether the Gelman-Rubin option was used in the routine.}
\item{fix.thin}{A TRUE/FALSE indicator denoting whether the thinning value was fixed by the user.}
\item{fix.burnin}{A TRUE/FALSE indicator denoting whether the burn-in value was fixed by the user.}
\item{fix.max}{A TRUE/FALSE indicator denoting whether the maximum value was fixed by the user.}
\item{InitialValues}{The initial values of the parameters at the start of the MCMC chain.}
\item{gamma.fixed}{Returns the fixed value of gamma or FALSE if the gamma parameters are sampled.}
\item{k.fixed}{Returns the fixed value of k or FALSE if k is sampled.}
\item{runtime}{The amount of time the routine ran in hours}
\item{outfolder}{The pathname of the folder containing the MCMC chains and the graphical output produced by the routine.}
\item{maxStudyTime}{The maximum time in the study (either as a censored observation or an observed failure)}
}

\references{
P. Bouman, J. Dignam, V. Dukic, XL. Meng.  (2005) Bayesian multiresolution hazard model with application to an aids reporting delay study.  \emph{Statistica Sinica}, \bold{102}, 1145--1157.

Bouman, P., Dignam, J., Dukic, V. (2007), A multiresolution hazard model for multi-center survival studies: Application to Tamoxifen treatment in early stage breast cancer. \emph{JASA}. \bold{102,} 1145--1157.

Dukic, V., Dignam, J. (2007), Bayesian hierarchical multiresolution hazard model for the study of time-dependent failure patterns in early stage breast cancer. \emph{Bayesian Analysis}. \bold{2,} 591--610.

Chen, Y., Hagar, Y., Dignam, J., Dukic, V. (2014), Pruned Multiresolution Hazard (PMRH) models for time-to-event data. \emph{In review.} Available upon request via email to vanja.dukic@colorado.edu.

http://amath.colorado.edu/faculty/vdukic/software/MRH.html

}
\author{
Yolanda Hagar <yolanda.hagar@colorado.edu>
}
\examples{
#####################################################
# NOTE: Examples may take a few minutes, so please be
# patient. In addition, warning messages about
# convergence may appear as more iterations are
# typically needed for chain convergence.  
	
########### PROPORTIONAL HAZARDS EXAMPLE ############
# Examine the NCCTG lung cancer data set (from the survival package),  
# and quantify how age, gender, and physician rated 
# Karnofsky performance scores affect survival times (in days).  
# Assume the hazards are proportional for all covariates.

data(cancer)

# Adjust "status" so that it is a 0/1 
# variable (currently it is 1 = censored, 2 = observed death)
cancer$censorvar = cancer$status - 1

# Run the estimateMRH routine.  Set the maximum 
# study time to 960 days, which makes each bin
# 120 days long.  This censors 0 extra subjects 
# (see FindBinWidth() for an example).  Save
# the output in a folder titled 'MRH_lung' 
# (default is 'MRHresults').
# Generally it is recommended to use a higher burn-in value, 
# thinning value, and maximum number
# of iterations, but for illustrative purposes 
# these values have been lowered.
# Note that the routine may produce a warning 
# message that the algorithm has not converged,
# as typically more iterations are needed for convergence.
# However, for the purposes of this example, the number
# of iterations is sufficient.
\dontrun{
fit.lung = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, outfolder = 'MRH_lung')}
\dontshow{
fit.lung = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 2, maxStudyTime = 960, burnIn = 50, maxIter = 100, 
	thin = 1, outfolder = 'MRH_lung')}

# See all items returned in the model fit
\dontrun{
	fit.lung
}
# See the main summary
\dontrun{
fit.lung$summary
}
# NOTE: If estimateMRH is run as a background job, 
# or if the output folder has been saved for use
# at a later instance, then  the fit can be calculated 
# using the as.MRH() and summary.MRH() functions.  
# See the those help pages or the vignette for 
# more information.

# Run the same model as above, but with pruning. 
# Save the output in a folder titled 'MRH_lung_prune'
\dontrun{
fit.lung.prune = estimateMRH(formula = Surv(time, censorvar) ~ 
	age + as.factor(sex) + ph.karno, data = cancer,
	M = 3, maxStudyTime = 960, burnIn = 200, maxIter = 1000, 
	thin = 1, prune = TRUE, outfolder = 'MRH_lung_prune')}


########### NON-PROPORTIONAL HAZARDS EXAMPLE ############
# Examine the tongue data set (from the KMsurv package), and
# quantify how the rumor DNA profile # affects survival time (in weeks).
data(tongue)

# Fit the MRH model, including tumor type using 
# the non-proportional hazards model.  
# With 16 bins (M = 4), each bin represents 25 weeks.
# Generally it is recommended to use a higher burn-in value, 
# thinning value, and maximum number
# of iterations, but for illustrative purposes 
# these values have been lowered.
# Note that the routine may produce a warning 
# message that the algorithm has not converged,
# as typically more iterations are needed for convergence.
# However, for the purposes of this example, the number
# of iterations is sufficient.
\dontrun{
fit.tongue = estimateMRH(formula = Surv(time, delta) ~ 
	nph(type), data = tongue, M = 4, 
	burnIn = 200, maxIter = 2000, thin = 1, 
	outfolder = 'MRH_tongue_nph')}

# Get the time-varying hazard ratios 
\dontrun{
	fit.tongue$beta
}
}
