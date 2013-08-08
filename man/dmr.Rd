\name{dmr}
\alias{dmr}
\alias{predict.dmr}
\title{Distributed Multinomial Regression}
\usage{
dmr(counts, covars, bins=NULL, 
                cores=1, k=2, ...)

\method{predict}{dmr}(object,newdata,
	type=c("link","response","reduction"))
}
\arguments{
\item{counts}{A dense \code{matrix} 
      or sparse \code{Matrix} of
      response counts. }
\item{covars}{A dense \code{matrix} 
      or sparse \code{Matrix} of covariates.
      This should not include the intercept.}
\item{bins}{Number of bins into which we will attempt to collapse each column of \code{covars}.  Since sums of multinomials 
with equal probabilities are also multinomial, the model is then fit to these collapsed `observations'. \code{bins=NULL}
 does no collapsing. }
\item{cores}{The number of shared memory processor cores to pass to \code{mclapply}.}
\item{k}{The information criteria penalty on degrees of freedom.  \code{k=2} is the AIC, \code{k=log(n)} is the BIC. }
\item{...}{Additional arguments to \code{gamlr}.}
\item{object}{A fitted \code{dmr} coefficient matrix.}
\item{newdata}{A Matrix with the same number of columns as \code{covars}, unless
\code{type="reduction"} in which case \code{newdata} is multinomial category count data with the same number of columns as \code{counts}.}
\item{type}{  Under "link", just the linear map \code{newdata} times \code{object}, under "response" the fitted multinomial probabilities, and under "reduction" the MNIR sufficient reduction \eqn{F\phi/m}.}
}
\details{
	\code{dmr} fits multinomial logistic regression by assuming that, unconditionally on the `size' (total count across categories) each individual category count has been generated as a Poisson
	\deqn{
	x_ij \sim Po(exp[\mu_i + \alpha_j + \beta v_i ]).
	}
	We plug-in estimate \eqn{\hat\mu_i = m_i/p}, where \eqn{m_i = \sum_j x_ij} and \eqn{p} is the dimension of \eqn{x_i}.  Then each individual is outsourced to Poisson regression in the \code{gamlr} package via the \code{mclapply} function of the \code{Parallel} library.  \code{gamlr} fits regularization paths of possible coefficients, and returns the best according to \code{AIC} with \code{k} as specificed here.  The combined coefficents across all dimensions are then returned from \code{dmr}.
}
\value{ A \code{dmr} covariate matrix, which inherits the \code{dgCMatrix} class as defined in the \code{Matrix} library.  This is the \code{ncol(covars)} by \code{ncol(counts)} matrix of logistic regression coefficients chosen by AIC from the regularization paths for each category. }
\references{Taddy (2013) Distributed Multinomial Regression}
\author{Matt Taddy \email{taddy@chicagobooth.edu}}
\seealso{Examples in \code{we8there}, and the \code{gamlr} package.}

