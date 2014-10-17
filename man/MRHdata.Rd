\name{MRHdata}
\alias{MRHdata}
\docType{data}
\title{
SIMULATED SURVIVAL DATA SET
}
\description{
This data set is a simulated survival data set, which is useful for code examples, particularly in the vignette.  In addition to the failure times and a censoring variable, the data set also includes gender, treatment, and age as covariates.  This treatment hazards are non-proportional, and there are periods of sparse observations and high rates of censoring, which is similar to studies with longer follow-up periods.
}
\usage{data(MRHdata)}
\format{
  A data frame with 3000 observations on the following 5 variables:
  \describe{
    \item{\code{time}}{A numeric vector containing the survival time from start of treatment to failure.}
    \item{\code{delta}}{A numeric vector of the censoring variable, which equals '1' if the failure is observed, and '0' otherwise.}
    \item{\code{gender}}{A numeric binary vector, with 0 denoting males and 1 denoting females.}
    \item{\code{treatment}}{A numeric vector with three treatment groups.  The treatment hazard rates were simulated under the non-proportional hazards assumption.}
    \item{\code{age}}{A numeric vector for age, with the measurements standardized.}
  }
}
