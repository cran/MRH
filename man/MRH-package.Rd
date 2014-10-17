\name{MRH-package}
\alias{MRH-package}
\docType{package}
\title{
	Multi-resolution estimation of the hazard rate.
}
\description{
	This package uses the multi-resolution hazard estimator for estimation and inference for the hazard rate.  The multi-resolution hazard estimator is a Bayesian semi-parametric method for estimating the hazard rate jointly with covariate effects.  This methodology splits the hazard rate into 2^M bins, and estimates a constant hazard rate within each bin using a tree-like structure, providing robust estimates of the hazard rate even through periods of sparse observations.  This package allows for covariates in the model with or without the proportional hazards assumption.
}
\details{
\tabular{ll}{
Package: \tab MRH\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-03\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Yolanda Hagar, Yuanting Chen, Vanja Dukic

Maintainer: Yolanda Hagar <yolanda.hagar@colorado.edu>
}
\references{
Bouman, P., Dukic, V., Meng, X.L. (2005), Bayesian multiresolution hazard model with application to an {AIDS} reporting delay study. \emph{Statistica Sinica}. \bold{15,} 325--357.

Bouman, P., Dignam, J., Dukic, V. (2007), A multiresolution hazard model for multi-center survival studies: Application to Tamoxifen treatment in early stage breast cancer. \emph{JASA}. \bold{102,} 1145--1157.

Dukic, V., Dignam, J. (2007), Bayesian hierarchical multiresolution hazard model for the study of time-dependent failure patterns in early stage breast cancer. \emph{Bayesian Analysis}. \bold{2,} 591--610.

Hagar, Y., Albers, D., Pivavarov, R., Chase, H., Dukic, V., Elhadad, N. (2014), Survival analysis with Electronic Health Record data: Experiments with Chronic Kidney Disease. \emph{Statistical Analysis and Data Mining}. \bold{To appear}.

Chen, Y., Hagar, Y., Dignam, J., Dukic, V. (2014), Pruned Multiresolution Hazard (PMRH) models for time-to-event data. \emph{In review.} Available upon request via email to vanja.dukic@colorado.edu.

This work utilized the Janus supercomputer, which is supported by the National Science Foundation (award number CNS-0821794) and the University of Colorado Boulder. The Janus supercomputer is a joint effort of the University of Colorado Boulder, the University of Colorado Denver and the National Center for Atmospheric Research.

Yolanda Hagar and Vanja Dukic were supported in part by grants NSF-DEB 1316334 and NSF-GEO 1211668.

}
\keyword{multi-resolution hazard, survival analysis, hazard rate, non-proportional hazards}