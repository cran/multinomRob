% $Id: multinomMLE.Rd,v 1.11 2005/09/23 03:12:42 wrm1 Exp $
%
\name{Multinomial Regression}
\alias{multinomMLE}
\title{Multinomial Regression Maximum Likelihood Estimator with Overdispersion}
\description{
  \code{multinomMLE} estimates the coefficients of the multinomial
  regression model for grouped count data by maximum likelihood, then
  computes a moment estimator for overdispersion and reports standard
  errors for the coefficients that take overdispersion into account.
  This function is not meant to be called directly by the user.  It is
  called by \code{multinomRob}, which constructs the various arguments.
}
\usage{multinomMLE(Y, Ypos, Xarray, xvec, jacstack, itmax=100, xvar.labels,
                   choice.labels, print.level)
}
\arguments{
  \item{Y}{
    Matrix (observations by alternatives) of outcome counts.
    Values must be nonnegative.  Missing data (\code{NA} values) are not allowed.}
  \item{Ypos}{
    Matrix indicating which elements of Y are counts to be analyzed (TRUE)
    and which are values to be skipped (FALSE).  This allows the set of
    outcome alternatives to vary over observations.}
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
  \item{itmax}{
    The maximum number of iterations to be done in the Gauss-Newton
    optimization.}
  \item{xvar.labels}{
    Vector of labels for observations.}
  \item{choice.labels}{
    Vector of labels for outcome alternatives.}
  \item{print.level}{
    Specify 0 for minimal printing (error messages only) or 3 to print details
    about the MLE computations.}
}
\details{
  Following the generalized linear models approach, the coefficient
  parameters in an overdispersed multinomial regression model may be
  estimated using the likelihood for a standard multinomial regression
  model.  A moment estimator may be used for the dispersion parameter, given
  the coefficient estimates, with little efficiency loss.
}
\value{
  multinomMLE returns a list of 20 objects.  The returned objects are:

   \item{coefficients}{
     The maximum likelihood coefficient estimates in matrix format.  The value 0 is
     used in the matrix to fill in for values that do not correspond to a
     regressor.}
   \item{coeffvec}{
     A vector containing the maximum likelihood coefficient estimates.}
   \item{dispersion}{
     Moment estimate of the dispersion:  mean sum of squared orthogonalized
     residuals (adjusted for degrees of freedom lost to estimated coefficients).}
   \item{se}{
     The tanh coefficient estimate standard errors in matrix format.  The
     format and labelling used for the matrix is the same as is used for the
     \code{coefficients}.  The standard errors are derived from the estimated
     asymptotic sandwich covariance estimate.}
   \item{se.opg}{
     The tanh coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the outer product of the gradient (expected
     information) divided by the moment estimate of the dispersion.}
   \item{se.hes}{
     The tanh coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the hessian matrix (observed information).}
   \item{se.vec}{
     A vector containing the tanh coefficient sandwich estimate standard errors.}
   \item{se.opg.vec}{
     A vector containing the tanh coefficient OPG estimate standard errors.}
   \item{A}{
     The outer product of the gradient (expected information) divided by the
     moment estimate of the dispersion.}
   \item{B}{
     The inverse of the hessian matrix (observed formation).}
   \item{covmat}{
     Sandwich estimate of the asymptotic covariance of the maximum likelihood
     coefficient estimates.}
   \item{iters}{
     Number of Gauss-Newton iterations.}
   \item{error}{
     Exit error code.}
   \item{GNlist}{
     List reporting final results of the Gauss-Newton optimization.  Elements:
     \code{coefficients}, vector of coefficient parameters (same as \code{coeffvec}
     value in list returned by multinomMLE);  \code{tvec}, matrix of coefficient
     parameters (same as \code{coefficients} value in list returned by multinomMLE);
     \code{formation}, inverse Hessian matrix;  \code{score}, score (or gradient
     element) matrix;  \code{LLvals}, list containing log-likelihood value;
     \code{convflag}, TRUE/FALSE convergence flag;  \code{iters}, number of
     iterations done in final Gauss-Newton stage;  \code{posdef}, TRUE if Hessian is
     positive definite.}
   \item{sigma2}{
     Moment estimate of the dispersion:  mean sum of squared orthogonalized
     residuals (adjusted for degrees of freedom lost to estimated coefficients).}
   \item{Y}{
     The same \code{Y} matrix that was supplied as input, except modified by 
     having done \code{Y[!Ypos] <- 0}.}
   \item{Ypos}{
     The same \code{Ypos} matrix that was supplied as input.}
   \item{fitted.prob}{
     The matrix of predicted probabilities for each category for each
     observation based on the coefficient estimates.}
   \item{jacstack}{
     The same \code{jacstack} that was supplied as an input argument.}
}
\references{
  Walter R. Mebane, Jr. and  Jasjeet Singh Sekhon. 2004.  ``Robust Estimation
  and Outlier Detection for Overdispersed Multinomial Models of Count Data.'' 
  \emph{American Journal of Political Science} 48 (April): 391--410 .
  \url{http://macht.arts.cornell.edu/wrm1/multinom.pdf}

  For the most current code and related material see
  \url{http://sekhon.polisci.berkeley.edu/robust/}
}
\author{Walter R. Mebane, Jr., Cornell University,
  \email{wrm1@cornell.edu}, \url{http://macht.arts.cornell.edu/wrm1/} \cr

  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.polisci.berkeley.edu/}}
}
\keyword{models}
\keyword{regression}
