\name{getmin}
\alias{getmin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getmin(lambda, cvm, cvsd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{cvm}{
%%     ~~Describe \code{cvm} here~~
}
  \item{cvsd}{
%%     ~~Describe \code{cvsd} here~~
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
function (lambda, cvm, cvsd) 
{
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + 3 * cvsd)[idmin]
    idmin = cvm <= semin
    lambda.1se = max(lambda[idmin], na.rm = TRUE)
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
