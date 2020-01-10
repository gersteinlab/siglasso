\name{signature}
\alias{signature}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("signature")}
\format{
  A data frame with 96 observations on the following 40 variables.
  \describe{
    \item{\code{Substitution.Type}}{a factor with levels \code{C>A} \code{C>G} \code{C>T} \code{T>A} \code{T>C} \code{T>G}}
    \item{\code{Trinucleotide}}{a factor with levels \code{ACA} \code{ACC} \code{ACG} \code{ACT} \code{ATA} \code{ATC} \code{ATG} \code{ATT} \code{CCA} \code{CCC} \code{CCG} \code{CCT} \code{CTA} \code{CTC} \code{CTG} \code{CTT} \code{GCA} \code{GCC} \code{GCG} \code{GCT} \code{GTA} \code{GTC} \code{GTG} \code{GTT} \code{TCA} \code{TCC} \code{TCG} \code{TCT} \code{TTA} \code{TTC} \code{TTG} \code{TTT}}
    \item{\code{Somatic.Mutation.Type}}{a factor with levels \code{A[C>A]A} \code{A[C>A]C} \code{A[C>A]G} \code{A[C>A]T} \code{A[C>G]A} \code{A[C>G]C} \code{A[C>G]G} \code{A[C>G]T} \code{A[C>T]A} \code{A[C>T]C} \code{A[C>T]G} \code{A[C>T]T} \code{A[T>A]A} \code{A[T>A]C} \code{A[T>A]G} \code{A[T>A]T} \code{A[T>C]A} \code{A[T>C]C} \code{A[T>C]G} \code{A[T>C]T} \code{A[T>G]A} \code{A[T>G]C} \code{A[T>G]G} \code{A[T>G]T} \code{C[C>A]A} \code{C[C>A]C} \code{C[C>A]G} \code{C[C>A]T} \code{C[C>G]A} \code{C[C>G]C} \code{C[C>G]G} \code{C[C>G]T} \code{C[C>T]A} \code{C[C>T]C} \code{C[C>T]G} \code{C[C>T]T} \code{C[T>A]A} \code{C[T>A]C} \code{C[T>A]G} \code{C[T>A]T} \code{C[T>C]A} \code{C[T>C]C} \code{C[T>C]G} \code{C[T>C]T} \code{C[T>G]A} \code{C[T>G]C} \code{C[T>G]G} \code{C[T>G]T} \code{G[C>A]A} \code{G[C>A]C} \code{G[C>A]G} \code{G[C>A]T} \code{G[C>G]A} \code{G[C>G]C} \code{G[C>G]G} \code{G[C>G]T} \code{G[C>T]A} \code{G[C>T]C} \code{G[C>T]G} \code{G[C>T]T} \code{G[T>A]A} \code{G[T>A]C} \code{G[T>A]G} \code{G[T>A]T} \code{G[T>C]A} \code{G[T>C]C} \code{G[T>C]G} \code{G[T>C]T} \code{G[T>G]A} \code{G[T>G]C} \code{G[T>G]G} \code{G[T>G]T} \code{T[C>A]A} \code{T[C>A]C} \code{T[C>A]G} \code{T[C>A]T} \code{T[C>G]A} \code{T[C>G]C} \code{T[C>G]G} \code{T[C>G]T} \code{T[C>T]A} \code{T[C>T]C} \code{T[C>T]G} \code{T[C>T]T} \code{T[T>A]A} \code{T[T>A]C} \code{T[T>A]G} \code{T[T>A]T} \code{T[T>C]A} \code{T[T>C]C} \code{T[T>C]G} \code{T[T>C]T} \code{T[T>G]A} \code{T[T>G]C} \code{T[T>G]G} \code{T[T>G]T}}
    \item{\code{Signature.1}}{a numeric vector}
    \item{\code{Signature.2}}{a numeric vector}
    \item{\code{Signature.3}}{a numeric vector}
    \item{\code{Signature.4}}{a numeric vector}
    \item{\code{Signature.5}}{a numeric vector}
    \item{\code{Signature.6}}{a numeric vector}
    \item{\code{Signature.7}}{a numeric vector}
    \item{\code{Signature.8}}{a numeric vector}
    \item{\code{Signature.9}}{a numeric vector}
    \item{\code{Signature.10}}{a numeric vector}
    \item{\code{Signature.11}}{a numeric vector}
    \item{\code{Signature.12}}{a numeric vector}
    \item{\code{Signature.13}}{a numeric vector}
    \item{\code{Signature.14}}{a numeric vector}
    \item{\code{Signature.15}}{a numeric vector}
    \item{\code{Signature.16}}{a numeric vector}
    \item{\code{Signature.17}}{a numeric vector}
    \item{\code{Signature.18}}{a numeric vector}
    \item{\code{Signature.19}}{a numeric vector}
    \item{\code{Signature.20}}{a numeric vector}
    \item{\code{Signature.21}}{a numeric vector}
    \item{\code{Signature.22}}{a numeric vector}
    \item{\code{Signature.23}}{a numeric vector}
    \item{\code{Signature.24}}{a numeric vector}
    \item{\code{Signature.25}}{a numeric vector}
    \item{\code{Signature.26}}{a numeric vector}
    \item{\code{Signature.27}}{a numeric vector}
    \item{\code{Signature.28}}{a numeric vector}
    \item{\code{Signature.29}}{a numeric vector}
    \item{\code{Signature.30}}{a numeric vector}
    \item{\code{X}}{a logical vector}
    \item{\code{X.1}}{a logical vector}
    \item{\code{X.2}}{a logical vector}
    \item{\code{X.3}}{a logical vector}
    \item{\code{X.4}}{a logical vector}
    \item{\code{X.5}}{a logical vector}
    \item{\code{X.6}}{a logical vector}
  }
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
data(signature)
## maybe str(signature) ; plot(signature) ...
}
\keyword{datasets}
