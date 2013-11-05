\name{collapse}
\alias{collapse}
\title{Data checking and binning}
\description{Collapses counts along equal levels of binned covariates.}
\usage{collapse(v,x,mu=NULL,bins=NULL,listx=TRUE)}
\arguments{
	\item{v}{Either matrix or Matrix of covariates (matches \code{covars} in \code{dmr}).}
	\item{x}{Either matrix or Matrix of multinomial counts, or a factor (matches \code{counts} in \code{dmr}).}
	\item{mu}{Possible pre-specified fixed effects for \code{dmr}; otherwise they are calculated here.}
	\item{bins}{The number of quantile bins into which we collapse \code{v}.  \code{bins=NULL} does no collapsing. }
	\item{listx}{Whether to split the columns of \code{x} out as a list for parallel
	processing.}
}
\details{
	For each column of \code{v}, aggregates
	the observations into \code{bins}  defined by their average value.  Both \code{v} and \code{x} are then collapsed according to levels of the interaction across implied bin-factors, and the number
	of observations in each bin is recorded as \code{n}.  Look at the code of the \code{dmr} function to see \code{collapse} used in practice.
}
\value{A list containing collapsed and formatted \code{v}, \code{x}, and \code{nbin}, along with \code{mu = log(rowMeans(x)+1)} (unless pre-specified, or zero for rows with zero total count).}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
