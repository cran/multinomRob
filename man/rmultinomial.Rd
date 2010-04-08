% $Id: rmultinomial.Rd,v 1.4 2005/09/23 03:12:42 wrm1 Exp $
%
\name{rmultinomial}
\alias{rmultinomial}
\title{Random Number Generator for the Multinomial Distribution}
\description{
  Generates a random count vector for one observation of a multinomial
  distribution for \code{n} trials with probability vector \code{pr}.
}
\usage{
rmultinomial(n = 5, pr = c(0.5, 0.5), long = FALSE)
}
\arguments{
  \item{n}{Number of trials.}
  \item{pr}{Probability vector.}
  \item{long}{\code{TRUE} to choose one generator, \code{FALSE} to choose another one.}
}
\details{
  Generates a random count vector for one observation of a multinomial
  distribution for \code{n} trials with probability vector \code{pr}.
}
\value{
  \item{x}{Vector of counts.}
}
\author{Walter R. Mebane, Jr., University of Michigan,
  \email{wmebane@umich.edu}, \url{http://www-personal.umich.edu/~wmebane} \cr

  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}
}
\note{
  This function is only used in the examples and not in the
  \code{multinomRob} code.
}
\examples{
 rmultinomial(10, c(.3, .3, .4));
}
\keyword{distribution}
\keyword{multivariate}
