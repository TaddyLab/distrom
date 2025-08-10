\name{dmrcoef-class}
\docType{class}
\alias{dmrcoef-class}
\alias{predict,dmrcoef-method}
\title{Class \code{"dmrcoef"}}
\description{
The extended \code{\link[Matrix]{Matrix}} class for output from \code{coef.dmr}.
}
\details{
This is the class for a covariate matrix from \code{dmr} regression;  it inherits the \code{\link[Matrix]{Matrix}} class as defined in the \code{Matrix} library.  
 In particular, this is the \code{ncol(covars)} by \code{ncol(counts)} matrix of logistic regression coefficients chosen in \code{coef.dmr} from the regularization paths for each category. }
\section{Objects from the Class}{
Objects can be created only by a call to the \code{coef.dmr} function.
}
\section{Slots}{
  \describe{
    \item{\code{i}:}{From \code{\link[Matrix]{Matrix}}: the row indices. }
    \item{\code{p}:}{From \code{\link[Matrix]{Matrix}}: the column pointers. }
    \item{\code{Dim}:}{From \code{\link[Matrix]{Matrix}}: the dimensions. }
    \item{\code{Dimnames}:}{From \code{\link[Matrix]{Matrix}}: the list of labels. }
    \item{\code{x}:}{From \code{\link[Matrix]{Matrix}}: the nonzero entries. }
    \item{\code{factors}:}{From \code{\link[Matrix]{Matrix}}.}
}
}
\section{Extends}{
  Class \code{\link[Matrix]{Matrix}}, directly.
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "dmrcoef")}: 
      Prediction for a given \code{dmrcoef} matrix.  Takes the same arguments as \code{predict.dmr}, but will be faster (since \code{coef.dmr} is called inside \code{predict.dmr}).  }
	 }
}
\seealso{\link{dmr}, \link{coef.dmr}, \link{predict.dmr}}
\author{Matt Taddy \email{mataddy@gmail.com} }
\examples{showClass("dmrcoef")}
\keyword{classes}
