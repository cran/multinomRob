% $Id: rmultinomial.Rd,v 1.3 2004/02/18 07:07:40 wrm1 Exp $
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
\author{Walter R. Mebane, Jr., Cornell University,
  \email{wrm1@cornell.edu}, \url{http://macht.arts.cornell.edu/wrm1/} \cr

  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
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