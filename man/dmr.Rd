\name{dmr}
\alias{dmr}
\alias{logLik.dmr}
\alias{predict.dmr}
\alias{coef.dmr}
\title{Distributed Multinomial Regression}
\description{Gamma-lasso path estimation for a multinomial logistic regression factorized into independent  Poisson log regressions.}
\usage{
dmr(counts, covars, bins=NULL, 
    lambda.start=NULL, cores=1,  ...)
\method{logLik}{dmr}(object, ...)
\method{coef}{dmr}(object, select=NULL, 
		grouped=TRUE, k=2, cores=attributes(object)['cores'], ...)
\method{predict}{dmr}(object,newdata,
	type=c("link","response","reduction"), ...)
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
\item{lambda.start}{Where to start each regularization path.  If \code{NULL} it uses the maximum absolute gradient across all categories (i.e. the smallest \eqn{\lambda} such that all coefficients are set to zero). }
\item{cores}{The number of shared memory processor 
	cores to pass to \code{parallel} library functions.  
	Type \code{detectCores()} to get a guess at how many cores you have available.  We use \code{mclapply} in \code{dmr} and \code{mcmapply} in \code{coef.dmr}.  These functions rely on unix forking, and will thus unfortunately not work in windows (you'll be stuck with \code{cores=1}). }
\item{select}{For \code{coef.dmr}, this is the index of the 
regularization paths from which you want estimates.  Can either be a single value for all categories or a vector of values, one for each category.  If left \code{NULL} the coefficients are selected via an information criteria according to arguments \code{k} and \code{grouped}.}
\item{k}{The information criteria penalty on degrees of freedom.  \code{k=2} is the AIC, \code{k=log(n)} is the BIC. }
\item{grouped}{For model selection in \code{coef.dmr} with \code{select=NULL}. Under \code{grouped=TRUE} we use a single shared \code{gamlr} penalty with minimum sum IC (information criterion) across response classes. Otherwise coefficients are returned at IC-optimal penalization for each individual response dimension. }
\item{type}{
For \code{predict.dmr}, this is the scale upon which you want prediction. Under "link", just the linear map \code{newdata} times \code{object}, under "response" the fitted multinomial probabilities, and under "reduction" the MNIR sufficient reduction \eqn{F\phi/m}.}
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
	We plug-in estimate \eqn{\hat\mu_i = log(m_i/p +1)}, where \eqn{m_i = \sum_j x_{ij}} and \eqn{p} is the dimension of \eqn{x_i}.  Then each individual is outsourced to Poisson regression in the \code{gamlr} package via the \code{mclapply} function of the \code{parallel} library.  The output from \code{dmr} is a list of \code{gamlr} fitted models.

	\code{coef.dmr} builds a matrix of multinomial logistic regression coefficients from the \code{length(object)} list of \code{gamlr} fits. 
	This construction can also be outsourced to multiple cores (using \code{mapply}).  If the \code{select} index is not null, \code{coef.dmr} returns the corresponding coefficients.  More likely, selection is based on an information criteria via \code{AIC}, with arguments \code{k} (complexity penalty) and \code{grouped} (for shared or individual \eqn{\lambda} across response categories; see argument description).  The combined coefficients across all dimensions are then returned as a \code{dmrcoef} s4-class object.

	\code{predict.dmr} takes either a \code{dmr} or \code{dmrcoef} object and returns predicted values for \code{newdata} on the scale defined by the \code{type} argument.  For \code{type="reduction"} this is a \code{data.frame} of sufficient reduction factors as in the MNIR paper, and otherwise it
	is a matrix of predicted class probabilities (\code{type="response"})
	or pre-logit linear equations (\code{type="link"}).
}
\value{  The \code{dmr} s3 object: an \code{ncol(counts)}-length list of fitted \code{gamlr} objects, with the added attributes \code{nobs} and \code{cores}. }
\references{
Taddy (2013) Distributed Multinomial Regression

Taddy (2013) The Gamma Lasso

Taddy (2013) Multinomial Inverse Regression for Text Analysis, with discussion and rejoinder, Journal of the American Statistical Association.
}
\author{Matt Taddy \email{taddy@chicagobooth.edu}}
\examples{

library(MASS)
data(fgl)
fits <- dmr(fgl$type, fgl[,1:9])
B <- coef(fits)
log(B@lambda)
P <- predict(B, fgl[,1:9], type="response")

## plot fitted probability by true response
boxplot(P[cbind(1:214,fgl$type)]~fgl$type, 
	ylab="fitted prob of true class")

## plot the individual Poisson model fits and selection
par(mfrow=c(3,2))
for(j in 1:6){
	plot(fits[[j]])
	mtext(names(fits)[j],font=2,line=2)
	abline(v=log(attributes(B)$lambda[j]), col="darkorange") }

}
\seealso{\code{dmrcoef-class}, \code{AIC},  and the \code{gamlr} package.}

