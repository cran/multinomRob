% $Id: multinomT.Rd,v 1.6 2005/09/23 03:12:42 wrm1 Exp $
%
\name{Multinomial Multivariate-T Regression}
\alias{multinomT}
\title{Multinomial Multivariate-T Estimation}
\description{
  \code{multinomT} fits the multinomial multivariate-t regression for grouped
  count data.  This function is not meant to be called directly by the
  user.  It is called by \code{multinomRob}, which constructs the
  various arguments.
}
\usage{multinomT(Yp, Xarray, xvec, jacstack, start = NA, nobsvec, fixed.df = NA)
}
\arguments{
  \item{Yp}{
    Matrix (observations by alternatives) of outcome proportions.
    Values must be between 0 and 1.  Missing data (\code{NA} values) are
    not allowed.}
  \item{Xarray}{
    Array of regressors.  dim(Xarray) = c(observations, parameters, alternatives).}
  \item{xvec}{
    Matrix (parameters by alternatives) that represents the model structure.
    It has a 1 for an estimated parameter, an integer greater than 1 for an
    estimated parameter constrained equal to another estimated parameter (all
    parameters constrained to be equal to one another have the same integer
    value in xvec) and a 0 otherwize.}
  \item{jacstack}{
    Array of regressors used to facilitate computing the gradient and the
    hessian matrix.
    dim(jacstack) = c(observations, unique parameters, alternatives).}
  \item{start}{
    A list of starting values of three kinds of parameters:
    \code{start$beta}, the values for the regression coefficients; \code{start$Omega}, the
    values for the variance-covariance matrix; \code{start$df}, the
    value for the multivariate-t degrees of freedom parameter.}
  \item{nobsvec}{
    Vector of the total number of counts for each observation.}
  \item{fixed.df}{
    The degrees of freedom to be used for the multivariate-t
    distribution.  When this is specified, the DF will not be estimated.}
}
\details{
  The function often provides good starting values for multinomRob's LQD
  estimator, but the standard errors it reports are not correct, in part
  because they ignore heteroscedasticity.
} 
\value{
  \item{call}{Names and values of all of the arguments which were passed
    to the function.  See \code{match.call} for further details.}
  \item{logL}{Log likelihood.}
  \item{deviance}{Deviance.}
  \item{par}{A list of three kinds of parameter estimates:
    \code{par$beta}, the estimates for the regression coefficients; \code{par$Omega}, the
    estimates for the variance-covariance matrix; \code{par$df}, the
    estimate of the multivariate-t degrees of freedom parameter.}
  \item{se}{Vector of standard errors for the regression coefficients.
    WARNING: these are not correct in part because the model ignores
    heteroscedasticity.}
  \item{optim}{Returned by \code{optim}.}
  \item{pred}{A matrix of predicted probabilities with the same
    dimentions as \code{Yp}.}
}
\references{
  Walter R. Mebane, Jr. and  Jasjeet Singh Sekhon. 2004.  ``Robust Estimation
  and Outlier Detection for Overdispersed Multinomial Models of Count Data.'' 
  \emph{American Journal of Political Science} 48 (April): 391--410.
  \url{http://sekhon.berkeley.edu/multinom.pdf} 
  
  For additional documentation please visit \url{http://sekhon.berkeley.edu/robust/}.
}
\author{Walter R. Mebane, Jr., University of Michigan,
  \email{wmebane@umich.edu}, \url{http://www-personal.umich.edu/~wmebane} \cr
  
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}
}
\seealso{
  \code{\link{match.call}}.  
  \code{\link{optim}}.
}
\keyword{models}
\keyword{regression}
