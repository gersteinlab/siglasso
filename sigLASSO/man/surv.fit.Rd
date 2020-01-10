\name{surv.fit}
\alias{surv.fit}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("surv.fit")}
\format{
  The format is:
List of 14
 $ n        : int [1:2] 28 6
 $ time     : int [1:34] 14 26 32 35 62 64 147 307 360 378 ...
 $ n.risk   : num [1:34] 28 27 26 25 24 23 22 21 20 19 ...
 $ n.event  : num [1:34] 0 0 0 0 0 0 0 0 0 0 ...
 $ n.censor : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
 $ surv     : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
 $ type     : chr "right"
 $ strata   : Named int [1:2] 28 6
  ..- attr(*, "names")= chr [1:2] "V4=0" "V4=1"
 $ std.err  : num [1:34] 0 0 0 0 0 0 0 0 0 0 ...
 $ upper    : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
 $ lower    : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
 $ conf.type: chr "log"
 $ conf.int : num 0.95
 $ call     : language survfit(formula = Surv(as.numeric(V3), as.numeric(V2)) ~ V4, data = mydata)
 - attr(*, "class")= chr "survfit"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(surv.fit)
## maybe str(surv.fit) ; plot(surv.fit) ...
}
\keyword{datasets}
