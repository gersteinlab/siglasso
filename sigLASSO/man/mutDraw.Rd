\name{mutDraw}
\alias{mutDraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mutDraw(mutationNum, simSpec, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutationNum}{
%%     ~~Describe \code{mutationNum} here~~
}
  \item{simSpec}{
%%     ~~Describe \code{simSpec} here~~
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
function (mutationNum, simSpec, noise) 
{
    noise_permute <- sample(seq(from = 1, to = 96), floor(runif(1, 
        min = 1, max = 25)))
    noise_random <- abs(rnorm(length(noise_permute)))
    noise_weights <- noise_random/sum(noise_random) * noise
    for (i in 1:length(noise_permute)) {
        simSpec[noise_permute[i]] <- simSpec[noise_permute[i]] + 
            noise_weights[i]
    }
    simAccm <- rep(NA, 96)
    simAccm[1] <- simSpec[1]
    for (i in 2:96) {
        simAccm[i] <- simAccm[i - 1] + simSpec[i]
    }
    simRandom <- runif(mutationNum)
    simMut <- rep(0, 96)
    simMut[1] <- length(which(simRandom <= simAccm[1]))
    for (i in 2:96) {
        simMut[i] = length(which(simRandom <= simAccm[i] & simRandom > 
            simAccm[i - 1]))
    }
    return(simMut)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
