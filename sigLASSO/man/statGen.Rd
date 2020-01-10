\name{statGen}
\alias{statGen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
statGen(N, mut_number, sig_number, noise, priorPer = 1, priorWeight = 0.1, Nsubsample = 0, adaptive = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{mut_number}{
%%     ~~Describe \code{mut_number} here~~
}
  \item{sig_number}{
%%     ~~Describe \code{sig_number} here~~
}
  \item{noise}{
%%     ~~Describe \code{noise} here~~
}
  \item{priorPer}{
%%     ~~Describe \code{priorPer} here~~
}
  \item{priorWeight}{
%%     ~~Describe \code{priorWeight} here~~
}
  \item{Nsubsample}{
%%     ~~Describe \code{Nsubsample} here~~
}
  \item{adaptive}{
%%     ~~Describe \code{adaptive} here~~
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
function (N, mut_number, sig_number, noise, priorPer = 1, priorWeight = 0.1, 
    Nsubsample = 0, adaptive = T) 
{
    generate <- randomMutGen(mut_number, sig_number, noise)
    spectrum <- unlist(generate[1])
    sigs <- unlist(generate[2])
    for (i in 1:(N - 1)) {
        generate <- randomMutGen(mut_number, sig_number, noise)
        spectrum <- cbind(unlist(generate[1]), spectrum)
        sigs <- cbind(unlist(generate[2]), sigs)
    }
    spectrum.deconstruct <- data.frame(t(spectrum)/mut_number)
    colnames(spectrum.deconstruct) <- colnames(sigs.input)
    rownames(spectrum.deconstruct) <- seq(1, N)
    coef_deconstruct <- NULL
    coef_all <- NULL
    coef_all2 <- NULL
    coef_all3 <- NULL
    coef_all4 <- NULL
    for (i in c(1:ncol(spectrum))) {
        print(i)
        penalty <- rep(1, ncol(sig[, -1]))
        if (priorPer > 0) {
            priorSet = which(sigs[, i] > 0)
            penalty[priorSet] = priorWeight
        }
        coef_all = cbind(coef_all, siglasso(spectrum[, i], data.matrix(sig[, 
            -1])))
        coef_all2 = cbind(coef_all2, siglasso(spectrum[, i], 
            data.matrix(sig[, -1]), conf = 10))
        coef_all3 = cbind(coef_all3, siglasso(spectrum[, i], 
            data.matrix(sig[, -1]), prior = penalty))
        coef_all4 = cbind(coef_all4, siglasso(spectrum[, i], 
            data.matrix(sig[, -1]), conf = 10, prior = penalty))
        coef_deconstruct <- cbind(coef_deconstruct, t(whichSignatures(tumor.ref = spectrum.deconstruct, 
            signatures.ref = signatures.cosmic, sample.id = i, 
            contexts.needed = TRUE, tri.counts.method = "default")[1]$weights))
    }
    return(list((coef_all - sigs)^2, support_stat(coef_all, sigs), 
        (coef_deconstruct - sigs)^2, support_stat(coef_deconstruct, 
            sigs), (coef_all2 - sigs)^2, support_stat(coef_all2, 
            sigs), coef_all, coef_deconstruct, sigs, coef_all2, 
        spectrum, (coef_all3 - sigs)^2, support_stat(coef_all3, 
            sigs), (coef_all4 - sigs)^2, support_stat(coef_all4, 
            sigs), coef_all3, coef_all4))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
