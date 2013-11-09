\name{cv.dmr}
\alias{cv.dmr}
\alias{logLik.cv.dmr}
\alias{plot.cv.dmr}
\alias{coef.cv.dmr}
\title{OOS experiments for distributed multinomial regression}
\description{OOS experiments for a multinomial logistic regression factorized into independent penalized Poisson log regressions.}
\usage{
cv.dmr(cl, covars, counts, 
    mu=NULL, lambda.start=NULL, 
    nfold=5, foldid=NULL,
    verb=1, savek=FALSE, ...)
\method{logLik}{cv.dmr}(object, ...)
\method{coef}{cv.dmr}(object, select=c("1se","min"), ...)
\method{plot}{cv.dmr}(x, ...)
}
\arguments{
\item{cl}{ A \code{parallel} library socket cluster; see details of the same argument for the \code{dmr} function. }
\item{covars}{A dense \code{matrix} 
      or sparse \code{Matrix} of covariates.
      This should not include the intercept.}
\item{counts}{A dense \code{matrix} 
      or sparse \code{Matrix} of
      response counts. }
\item{mu}{Possible pre-specified fixed effects, as in the same argument to \code{dmr}.}
\item{lambda.start}{Where to start each regularization path.  If \code{NULL} it uses the maximum absolute gradient across all category (i.e. the smallest \eqn{\lambda} such that all coefficients are set to zero). }
\item{nfold}{ The number of cross validation folds. }
\item{foldid}{ An optional length-n vector of fold memberships for each observation.  If specified, this dictates \code{nfold}.}
\item{verb}{ Whether to print progress through folds.
 \code{max(0,verb-1)} is passed
on to dmr. }
\item{savek}{ Whether to save each fold 'k' dmr fit.}
\item{...}{ Arguments to \code{dmr}. }
\item{select}{ In prediction and coefficient extraction, 
  select which "best" model to return: 
  \code{select="min"} is that with minimum average OOS deviance,
  and  \code{select="1se"} is that whose average OOS deviance is
  no more than 1 standard error away from the minimum.}
\item{x}{A \code{cv.dmr} object.}
\item{object}{A \code{cv.dmr} object.}
}
\details{ Fits a \code{dmr} regression to the full dataset, and then performs \code{nfold} 
cross validation to evaluate out-of-sample (OOS)
multinomial deviance under different penalty weights.  Be aware: this requires
much extra communication across response categories; for big datasets use \code{dmr} with AIC or BIC.

Model selection
here is grouped: the same \eqn{\lambda} is assumed for each response category.  This
is in contrast to the core \code{dmr} function, where coefficients are selected
independently via AIC for each response category (based on Poisson deviance).  For this reason, \code{cv.dmr} can take a very long time for high-dimensional counts (e.g., the we8there data).

\code{plot.cv.dmr} can be used to plot the results: it 
shows mean OOS multinomial deviance with 1se error bars. Average degrees of
freedom \emph{per response category} are printed along the top.  This 
calls \code{plot.cv.gamlr} and takes arguments to that function.

\code{coef.cv.dmr} returns a \code{dmrcoef} class matrix of coefficients; these
can be used in, say, \code{predict.dmr}.

\code{logLik.cv.dmr} allows for \code{AIC} calculations
on the grouped multinomial deviance (as opposed to the
independent Poissons used in \code{logLik.dmr}).
} 
\value{
  \item{dmr}{ The full-data fitted \code{dmr} object.}
  \item{nobs}{ The number of observations.}
  \item{dev}{ The full-data fit multinomial deviance. }
  \item{df}{ The matrix of degrees of freedom by response category and lambda.}
  \item{lambda}{ The grid of lambda values (shared across response categories).}
  \item{nfold}{ The number of CV folds. }
  \item{foldid}{ The length-n vector of fold memberships. }
  \item{cvm}{ Mean OOS deviance by \code{lambda} }
  \item{cvs}{ The standard errors on \code{cvm}.}
  \item{seg.min}{ The index of minimum \code{cvm}. }
  \item{seg.1se}{ The index of \code{1se} \code{cvm} (see details). }
  \item{lambda.min}{ Penalty at minimum \code{cvm}. }
  \item{lambda.1se}{ Penalty at \code{1se} \code{cvm}. }
  \item{kfit}{If \code{savek}, the \code{nfold} trained fits.}
 }
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\references{Taddy (2013), Distributed Multinomial Regression}

\examples{
library(MASS)
data(fgl)
## just fit in serial this time
fits <- cv.dmr(cl=NULL, fgl[,1:9], fgl$type)

plot(fits)
abline(v=log(fits$lambda[which.min(AIC(fits))]), 
		col="darkorange", lty=3)
abline(v=log(fits$lambda[which.min(BIC(fits))]), 
		col="darkorange", lty=3)
}
\seealso{dmr}

