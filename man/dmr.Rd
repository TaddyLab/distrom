\name{dmr}
\alias{dmr}
\alias{logLik.dmr}
\alias{predict.dmr}
\alias{coef.dmr}
\title{Distributed Multinomial Regression}
\description{Gamma-lasso path estimation for a multinomial logistic regression factorized into independent  Poisson log regressions.}
\usage{
dmr(covars, counts, mu=NULL, bins=NULL, cl=NULL, ...)
\method{logLik}{dmr}(object, ...)
\method{coef}{dmr}(object, select=NULL, k=2, ...)
\method{predict}{dmr}(object, newdata,
	type=c("link","response","class","reduction"), ...)
}
\arguments{
\item{covars}{A dense \code{matrix} 
      or sparse \code{Matrix} of covariates.
      This should not include the intercept.}
\item{counts}{A dense \code{matrix} 
      or sparse \code{Matrix} of
      response counts. }
\item{mu}{
	Pre-specified fixed effects for each observation in the Poisson regression linear equation.  If \code{mu=NULL}, then this is fixed as \code{log(rowMeans(x)+1)}.  Note that if \code{bins} is non-null then this argument is ignored and \code{mu} is recalculated on the collapsed data.}
\item{bins}{Number of bins into which we will attempt to collapse each column of \code{covars}.  Since sums of multinomials 
with equal probabilities are also multinomial, the model is then fit to these collapsed `observations'. \code{bins=NULL}
 does no collapsing. }
\item{cl}{A \code{parallel} library socket cluster.  If \code{is.null(cl)},
 a cluster is initialized with 

 \code{makeCluster(min(detectCores(),ncol(x)),
	type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK"))} 
	and stopped via \code{stopCluster(cl)} upon exit.  See
	\code{help(parallel)} for details. }
\item{select}{For \code{coef.dmr}, this is the index of the 
regularization paths from which you want estimates.  Can either be a single value for all categories or a vector of values, one for each category.  If left \code{NULL} the coefficients are selected via an information criteria according to the argument \code{k}.}
\item{k}{For model selection when \code{select=NULL}, the information criteria penalty on degrees of freedom.  \code{k=2} is the AIC, \code{k=log(n)} is the BIC. }
\item{type}{
For \code{predict.dmr}, this is the scale upon which you want prediction. Under "link", just the linear map \code{newdata} times \code{object}, under "response" the fitted multinomial probabilities, under "class" the max-probability class label, and under "reduction" the MNIR sufficient reduction \eqn{X\phi/m}, plus an extra column \eqn{m} containing row totals of \eqn{X}.}
\item{newdata}{A Matrix with the same number of columns as \code{covars}, unless
\code{type="reduction"} in which case \code{newdata} is multinomial category count data with the same number of columns as \code{counts}.}
\item{...}{Additional arguments to \code{gamlr} from \code{dmr} or to \code{coef.dmr} from \code{predict.dmr}.}
\item{object}{A \code{dmr} list of fitted \code{gamlr} models for each response category. }
}
\details{
	\code{dmr} fits multinomial logistic regression by assuming that, unconditionally on the `size' (total count across categories) each individual category count has been generated as a Poisson
	\deqn{
	x_{ij} \sim Po(exp[\mu_i + \alpha_j + \beta v_i ]).
	}
	We [default] plug-in estimate \eqn{\hat\mu_i = log(m_i/p +1)}, where \eqn{m_i = \sum_j x_{ij}} and \eqn{p} is the dimension of \eqn{x_i}.  Then each individual is outsourced to Poisson regression in the \code{gamlr} package via the \code{parLapply} function of the \code{parallel} library.  The output from \code{dmr} is a list of \code{gamlr} fitted models.

	\code{coef.dmr} builds a matrix of multinomial logistic regression coefficients from the \code{length(object)} list of \code{gamlr} fits. If the \code{select} index is not null, \code{coef.dmr} returns the corresponding coefficients.  More likely, selection uses an information criteria via \code{AIC}
	on Poisson deviance for each individual response dimension. The complexity penalty is specified via argument \code{k}.  Combined coefficients across all dimensions are then returned as a \code{dmrcoef} s4-class object.

	\code{predict.dmr} takes either a \code{dmr} or \code{dmrcoef} object and returns predicted values for \code{newdata} on the scale defined by the \code{type} argument.  
}
\value{  \code{dmr} returns the \code{dmr} s3 object: an \code{ncol(counts)}-length list of fitted \code{gamlr} objects, with the added attributes \code{nlambda}, \code{mu}, and \code{nobs}. }
\references{
Taddy (2013) Distributed Multinomial Regression

Taddy (2013) The Gamma Lasso

Taddy (2013) Multinomial Inverse Regression for Text Analysis, with discussion and rejoinder, Journal of the American Statistical Association.
}
\author{Matt Taddy \email{taddy@chicagobooth.edu}}
\examples{

library(MASS)
data(fgl)

fits <- dmr(fgl[,1:9], fgl$type)

## Individual Poisson model fits and IC selection
par(mfrow=c(3,2))
for(j in 1:6){
	plot(fits[[j]])
	mtext(names(fits)[j],font=2,line=2) }

##  AIC model selection
B <- coef(fits)
log(B@lambda)

## Fitted probability by true response
par(mfrow=c(1,1))
P <- predict(B, fgl[,1:9], type="response")
boxplot(P[cbind(1:214,fgl$type)]~fgl$type, 
	ylab="fitted prob of true class")


}
\seealso{\code{dmrcoef-class}, \code{cv.dmr}, \code{AIC},  and the \code{gamlr} package.}

