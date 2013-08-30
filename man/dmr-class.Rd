\name{dmr-class}
\docType{class}
\alias{dmr-class}
\alias{predict,dmr-method}
\title{Class \code{"dmr"}}
\description{
The extended \code{dgCMatrix} for output from \code{dmr}.
}
\section{Objects from the Class}{
Objects can be created only by a call to the \code{dmr} function.
}
\section{Slots}{
  \describe{
    \item{\code{lambda}:}{A vector of the selected \code{gamlr} penalty for each coefficient.} 
    \item{\code{i}:}{From \code{dgCMatrix}: the row indices. }
    \item{\code{p}:}{From \code{dgCMatrix}: the column pointers. }
    \item{\code{Dim}:}{From \code{dgCMatrix}: the dimensions. }
    \item{\code{Dimnames}:}{From \code{dgCMatrix}: the list of labels. }
    \item{\code{x}:}{From \code{dgCMatrix}: the nonzero entries. }
    \item{\code{factors}:}{From \code{dgCMatrix}.}
}
}
\section{Extends}{
  Class \code{\linkS4class{dgCMatrix}}, directly.
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "dmr")}: see \code{help(predict.dmr)} }
	 }
}
\seealso{dmr,predict.dmr}
\author{Matt Taddy \email{taddy@chicagobooth.edu} }
\examples{showClass("dmr")}
\keyword{classes}
