\name{dmrcoef-class}
\docType{class}
\alias{dmrcoef-class}
\alias{predict,dmrcoef-method}
\title{Class \code{"dmrcoef"}}
\description{
The extended \code{dgCMatrix} class for output from \code{coef.dmr}.
}
\details{
This is the class for a covariate matrix from \code{dmr} regression;  it inherits the \code{dgCMatrix} class as defined in the \code{Matrix} library.  
 In particular, this is the \code{ncol(covars)} by \code{ncol(counts)} matrix of logistic regression coefficients chosen in \code{coef.dmr} from the regularization paths for each category. }
\section{Objects from the Class}{
Objects can be created only by a call to the \code{coef.dmr} function.
}
\section{Slots}{
  \describe{
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
    \item{predict}{\code{signature(object = "dmrcoef")}: 
      Prediction for a given \code{dmrcoef} matrix.  Takes the same arguments as \code{predict.dmr}, but will be faster (since \code{coef.dmr} is called inside \code{predict.dmr}).  }
	 }
}
\seealso{dmr, coef.dmr, predict.dmr}
\author{Matt Taddy \email{taddy@chicagobooth.edu} }
\examples{showClass("dmrcoef")}
\keyword{classes}
