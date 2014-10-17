\name{Prune}
\alias{Prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Merges bins in the MRH model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates an indicator variable that denotes which bins can be merged in the MRH model based on a modified Fisher's exact test. 
}
\usage{
Prune(Ti, delta, M, NumPrune, maxStudyTime = NULL, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ti}{
%%     ~~Describe \code{Ti} here~~
The vector of failure times.
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
The censoring variable for each subject, which is equal to 1 if the event of interest was observed, and 0 if the event was not observed (i.e. the subject was censored).  
}
  \item{M}{
%%     ~~Describe \code{M} here~~
M denotes the depth of the multi-resolution tree, and is directly related to the maximum number of bins, 2^M.  
}
  \item{NumPrune}{
%%     ~~Describe \code{NumPrune} here~~
The number of levels in the MRH model subject to pruning, starting from the level with the most bins (the bottom level). For example, M=5 and NumPrune=3 means the 3, 4 and 5 level of  this five-level MRH model are subject to pruning.
}
  \item{maxStudyTime}{
%%     ~~Describe \code{censortime} here~~
The maximum time (either censored or observed) in the study, or the end of the study period.  
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
The significance level for Fisher's exact test.  The default value is 0.05.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Two adjacent bins that are constructed via the same split parameter (Rm,p, or 1-Rm,p) are merged if the estimated hazard rate in these two bins are statistically similar. Pruning hypothesis tests (with the null hypothesis that the two bins have equal hazard rates) are performed using a modified Fisher's exact test, based on the 2 by 2 table composed of the number of failures within the time interval and at-risk patients at the end of the  time interval, for each pair of adjacent bins sharing a split parameter. If the test fails to reject the null hypothesis at the alpha level, that split parameter (Rm,p) is set to 0.5.  An indicator vector is returned to the user, with  '1' for each  Rmp/split value  set to 0.5 (related to the bins that are fused), and a '0' otherwise.
}
\value{
An indicator vector denoting splits that are associated with merged bins with a '1' and splits that are not with a '0'.
}
\references{
Y. Chen, Y. Hagar, J. Dignam, V. Dukic.  (2013) Pruned Multiresolution Hazard (PMRH) Models for Time-to-Event Data.  \emph{Submitted; available upon request from vanja.dukic@colorado.edu}.}
\author{
%%  ~~who you are~~
 Yuanting Chen, Yolanda Hagar <yolanda.hagar@colorado.edu>, and Vanja Dukic 
}
