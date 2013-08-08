\name{we8there}
\alias{we8there}
\alias{we8thereRatings}
\alias{we8thereCounts}

\title{ On-Line Restaurant Reviews }

\description{
Counts for 2804 bigrams in 6175 restaurant reviews from the site www.we8there.com. }

\details{ The short user-submitted reviews are accompanied by a five-star rating on four specific aspects of restaurant quality - food, service, value, and atmosphere - as well as the overall experience.  The reviews originally appear in Maua and Cozman (2009), and the parsing details behind these specific counts are in Taddy (2011). }

\value{
   \item{we8thereCounts}{ A \code{dgCMatrix} of phrase counts indexed by review-rows and bigram-columns.}
   \item{we8thereRatings}{ A \code{matrix} containing the associated review ratings. }
}

\references{ 
Maua, D.D. and Cozman, F.G. (2009), \emph{Representing and classifying user reviews}. In ENIA '09: VIII Enconro Nacional de Inteligencia Artificial, Brazil.

Taddy (2013), \emph{Multinomial Inverse Regression for Text Analysis}.
\url{http://arxiv.org/abs/1012.2098}

Taddy (2013), \emph{Distributed Multinomial Regression}.

}

\author{ 
Matt Taddy, \email{taddy@chicagobooth.edu}
}

\seealso{ dmr }

\examples{
## some multinomial inverse regression
data(we8there)
B <- dmr(we8thereCounts, we8thereRatings$Overall, bins=5,cores=4)
Z <- predict(B,we8thereCounts,type="reduction")
summary(fwd <- lm(we8thereRatings$Overall ~ Z))
plot(fwd$fitted ~ factor(we8thereRatings$Overall))
}
