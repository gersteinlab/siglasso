\name{load_sig}
\alias{load_sig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load_sig()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    signature <- read.table("./data/signatures_probabilities.txt", 
        sep = "\t", header = T)
    sig <- signature[, 3:33]
    sig$Somatic.Mutation.Type <- gsub("\\[", "", sig$Somatic.Mutation.Type)
    sig$Somatic.Mutation.Type <- gsub("]", "", sig$Somatic.Mutation.Type)
    sig$Somatic.Mutation.Type <- gsub(">", "", sig$Somatic.Mutation.Type)
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 4, 
        4)), ]
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 1, 
        1)), ]
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 2, 
        3)), ]
    return(data.matrix(sig[, -1]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
