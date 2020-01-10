\name{randomMutGen}
\alias{randomMutGen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randomMutGen(mutationNum, simNum, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationNum}{
%%     ~~Describe \code{mutationNum} here~~
}
  \item{simNum}{
%%     ~~Describe \code{simNum} here~~
}
  \item{noise}{
%%     ~~Describe \code{noise} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mutationNum, simNum, noise) 
{
    min_weight <- 0.05
    min_weight <- 0.02
    simSig <- sample(seq(2, 31), simNum)
    simPercentPool <- sort(runif(simNum - 1, max = 1 - noise - 
        simNum * min_weight), decreasing = T)
    simSigProfile <- rep(0, 30)
    simSpec <- sig[, simSig[1]] * (1 - noise - simNum * min_weight - 
        simPercentPool[1] + min_weight)
    simSigProfile[simSig[1] - 1] <- 1 - noise - simNum * min_weight - 
        simPercentPool[1] + min_weight
    simSpec <- simSpec + sig[, simSig[simNum]] * (simPercentPool[simNum - 
        1] + min_weight)
    simSigProfile[simSig[simNum] - 1] <- simPercentPool[simNum - 
        1] + min_weight
    if (simNum > 2) {
        for (i in 2:(simNum - 1)) {
            simSpec <- simSpec + sig[, simSig[i]] * (simPercentPool[i - 
                1] - simPercentPool[i] + min_weight)
            simSigProfile[simSig[i] - 1] <- simPercentPool[i - 
                1] - simPercentPool[i] + min_weight
        }
    }
    return(list(mutDraw(mutationNum, simSpec, noise), simSigProfile, 
        simSpec))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
