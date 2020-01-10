\name{siglasso}
\alias{siglasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
siglasso(sample_spectrum, signature, conf = 1, prior, adaptive = T, gamma = 1, alpha_min = 400, iter_max = Inf, sd_multiplier = 0.5, elastic_net = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_spectrum}{
%%     ~~Describe \code{sample_spectrum} here~~
}
  \item{signature}{
%%     ~~Describe \code{signature} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
  \item{prior}{
%%     ~~Describe \code{prior} here~~
}
  \item{adaptive}{
%%     ~~Describe \code{adaptive} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{alpha_min}{
%%     ~~Describe \code{alpha_min} here~~
}
  \item{iter_max}{
%%     ~~Describe \code{iter_max} here~~
}
  \item{sd_multiplier}{
%%     ~~Describe \code{sd_multiplier} here~~
}
  \item{elastic_net}{
%%     ~~Describe \code{elastic_net} here~~
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
function (sample_spectrum, signature, conf = 1, prior, adaptive = T, 
    gamma = 1, alpha_min = 400, iter_max = Inf, sd_multiplier = 0.5, 
    elastic_net = F) 
{
    if (missing(sample_spectrum)) {
        stop("siglasso(spectrum, signature, prior, adaptive=T, gamma=1, \n                    alpha_min=400, iter_max=20, sd_multiplier=0.5")
    }
    if (missing(signature)) {
        print("No signature supplied, will use the COSMIC signatures")
        signature = load_sig()
    }
    if (nrow(signature) != length(sample_spectrum)) {
        stop("The number of rows of signatures do not equal the number of \n            rows of the sample spectrum")
    }
    if (length(unique(sample_spectrum)) == 1) {
        if (unique(sample_spectrum) == 0) {
            stop("The spectrum is an empty vector")
        }
        random_idx <- ceiling(length(sample_spectrum) * runif(2))
        sample_spectrum[random_idx[1]] <- sample_spectrum[random_idx[1]] - 
            1
        sample_spectrum[random_idx[2]] <- sample_spectrum[random_idx[2]] + 
            1
    }
    if (missing(prior)) {
        prior = rep(1, ncol(signature))
    }
    else if (length(prior) != ncol(signature)) {
        stop("The length of prior does not equal the number of signatures!")
    }
    return(siglasso_internal(sample_spectrum, signature, prior, 
        adaptive, elastic_net, gamma, alpha_min, iter_max, sd_multiplier, 
        conf))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
