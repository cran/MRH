\name{FindBinWidth}
\alias{FindBinWidth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A pre-processing function that calculates the length of time per bin for different values of M.  
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Before fitting a MRH model, the optimal number of bins must be determined.  The MRH methodology divides the total study time in to $2^M$ bins, so the choice of M can be determined through biological rationale or can be based on the ideal length of time per bin.  In some instances, the number of bins may be relatively easy to determine.  However, there are many cases where it is not clear what the ideal bin length should be.  In these instances, the FindBinWidth() provides a table of lengths of time per bin for different units of time (ranging from seconds to years) for different M values (ranging from M = 2 to M = 10).
}
\usage{
FindBinWidth(time, delta, time.unit, maxStudyTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{
%%     ~~Describe \code{time} here~~
The vector of failure times for the subjects in the study.
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
The censoring vector, with a '1' denoting the an observed failure, and a '0' denoting a censored failure.
}
  \item{time.unit}{
%%     ~~Describe \code{time.unit} here~~
The unit of time for the failures, with options for seconds ('s'), minutes ('min'), days ('d'), weeks ('w'), months ('mon'), and years ('y').
}
  \item{maxStudyTime}{
%%     ~~Describe \code{censortime} here~~
The maximum failure time (observed or censored) or the length of the study.
}
}
\author{
%%  ~~who you are~~
Yolanda Hagar
}
\examples{
# Examine the options for the NCCTG lung cancer data set (from the survival package)
data(cancer)

# Code the censoring variable delta as 0/1 instead of 1/2
cancer$censorvar = cancer$status - 1

# The time unit in the cancer data set is in days, so specify time.unit as "d". 
FindBinWidth(cancer$time, cancer$censorvar, time.unit = 'd')

# None of the bin options show an optimal/rounded length.  
# Set maxStudyTime to 960.  This will show the results if we use 8 bins, with 120 days per bin.
FindBinWidth(cancer$time, cancer$censorvar, time.unit = 'd', maxStudyTime = 8*120)

# Results show rounded bin lengths, and we see that with the shortened maximum study time
# zero extra failures are censored.
}
