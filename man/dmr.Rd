\name{dmr}
\alias{dmr}
\alias{distrom}
\alias{predict.dmr}
\alias{coef.dmr}
\title{Distributed Multinomial Regression}
\description{Gamma-lasso path estimation for a multinomial logistic regression factorized into independent  Poisson log regressions.}
\usage{
dmr(cl, covars, counts, mu=NULL, bins=NULL, verb=0, cv=FALSE, ...)
\method{coef}{dmr}(object, ...)
\method{predict}{dmr}(object, newdata,
	type=c("link","response","class"), ...)
}
\arguments{
\item{cl}{A \code{parallel} library socket cluster.  If \code{is.null(cl)},
 everything is done in serial.  See
	\code{help(parallel)}, \code{help(makeCluster)}, and our examples here for details. }
\item{covars}{A dense \code{matrix} 
      or sparse \code{Matrix} of covariates.
      This should not include the intercept.}
\item{counts}{A dense \code{matrix} 
      or sparse \code{Matrix} of
      response counts. }
\item{mu}{
	Pre-specified fixed effects for each observation in the Poisson regression linear equation.  If \code{mu=NULL}, then we use \code{log(rowSums(x))}.  Note that if \code{bins} is non-null then this argument is ignored and \code{mu} is recalculated on the collapsed data.}
\item{bins}{Number of bins into which we will attempt to collapse each column of \code{covars}.  Since sums of multinomials 
with equal probabilities are also multinomial, the model is then fit to these collapsed `observations'. \code{bins=NULL}
 does no collapsing. }
\item{verb}{Whether to print some info.  \code{max(0,verb-1)} is passed
on to gamlr and will print if you created an \code{outfile} when specifying \code{cl}.  }
\item{cv}{A flag for whether to use \code{cv.gamlr} instead of \code{gamlr} for each Poisson regression. }
\item{type}{
For \code{predict.dmr}, this is the scale upon which you want prediction. Under "link", just the linear map \code{newdata} times \code{object}, under "response" the fitted multinomial probabilities, under "class" the max-probability class label.  For sufficient reductions see the \code{srproj} function of the textir library.}
\item{newdata}{A Matrix with the same number of columns as \code{covars}.}
\item{...}{Additional arguments to \code{gamlr}, \code{cv.gamlr}, and their associated methods.  }
\item{object}{A \code{dmr} list of fitted \code{gamlr} models for each response category. }
}
\details{
	\code{dmr} fits multinomial logistic regression by assuming that, unconditionally on the `size' (total count across categories) each individual category count has been generated as a Poisson
	\deqn{
	x_{ij} \sim Po(exp[\mu_i + \alpha_j + \beta v_i ]).
	}
	We [default] plug-in estimate \eqn{\hat\mu_i = log(m_i)}, where \eqn{m_i = \sum_j x_{ij}} and \eqn{p} is the dimension of \eqn{x_i}.  Then each individual is outsourced to Poisson regression in the \code{gamlr} package via the \code{parLapply} function of the \code{parallel} library.  The output from \code{dmr} is a list of \code{gamlr} fitted models.

	\code{coef.dmr} builds a matrix of multinomial logistic regression
	coefficients from the \code{length(object)} list of \code{gamlr} fits.  Default selection under \code{cv=FALSE} uses an
	information criteria via \code{AICc} on Poisson deviance for each
	individual response dimension (see \code{gamlr}).  Combined coefficients
	across all dimensions are then returned as a \code{dmrcoef} s4-class
	object.

	\code{predict.dmr} takes either a \code{dmr} or \code{dmrcoef} object and returns predicted values for \code{newdata} on the scale defined by the \code{type} argument.  
}
\value{  \code{dmr} returns the \code{dmr} s3 object: an \code{ncol(counts)}-length list of fitted \code{gamlr} objects, with the added attributes \code{nlambda}, \code{mu}, and \code{nobs}. }
\references{
Taddy (2015 AoAS) Distributed Multinomial Regression

Taddy (2017 JCGS) One-step Estimator Paths for Concave Regularization, the Journal of Computational and Graphical Statistics

Taddy (2013 JASA) Multinomial Inverse Regression for Text Analysis
}
\author{Matt Taddy \email{mataddy@gmail.com}}
\examples{

library(MASS)
data(fgl)

## make your cluster 
## FORK is faster but memory heavy, and doesn't work on windows.
cl <- makeCluster(2,type=ifelse(.Platform$OS.type=="unix","FORK","PSOCK")) 
print(cl)

## fit in parallel
fits <- dmr(cl, fgl[,1:9], fgl$type, verb=1)

## its good practice stop the cluster once you're done
stopCluster(cl)

## Individual Poisson model fits and AICc selection
par(mfrow=c(3,2))
for(j in 1:6){
	plot(fits[[j]])
	mtext(names(fits)[j],font=2,line=2) }

##  AICc model selection
B <- coef(fits)

## Fitted probability by true response
par(mfrow=c(1,1))
P <- predict(B, fgl[,1:9], type="response")
boxplot(P[cbind(1:214,fgl$type)]~fgl$type, 
	ylab="fitted prob of true class")


}
\seealso{\code{dmrcoef-class}, \code{cv.dmr}, \code{AICc},  and the \code{gamlr} and \code{textir} packages.}

