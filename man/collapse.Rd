\name{collapse}
\alias{collapse}
\title{Data checking and binning}
\description{Collapses counts along equal levels of binned covariates.}
\usage{collapse(x,v,nquant=NULL)}
\arguments{
	\item{x}{Either matrix or Matrix of multinomial counts.}
	\item{v}{Either matrix or Matrix of covariates.}
	\item{nquant}{The number of quantile bins into which we collapse \code{v}.  \code{nquant=NULL} does no collapsing. }
}
\details{
	For each column of \code{v}, aggregates
	the observations into \code{nquant} bins defined by their average value.  Both \code{v} and \code{x} are then collapsed according to levels of the interaction across implied bin-factors, and the number
	of observations in each bin is recorded as \code{n}.
}
\value{A list containing collapsed and formatted \code{x}, \code{v}, and \code{n}.}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
