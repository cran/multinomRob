% $Id: mGNtanh.Rd,v 1.7 2005/09/23 03:12:42 wrm1 Exp $
%
\name{Multinomial Regression Tanh Estimator Gauss-Newton Optimization}
\alias{mGNtanh}
\title{Multinomial Regression Hyperbolic Tangenet (Tanh) Estimator Gauss-Newton Optimization}
\description{
  \code{mGNtanh} uses Gauss-Newton optimization to compute the
  hyperbolic tangent (tanh) estimator for the overdispersed multinomial
  regression model for grouped count data.  This function is not meant
  to be called directly by the user.  It is called by
  \code{multinomRob}, which constructs the various arguments.
}
\usage{mGNtanh(bstart, sigma2, resstart, Y, Ypos, Xarray, xvec, tvec,
               jacstack, itmax = 100, print.level = 0)
}
\arguments{
  \item{bstart}{
    Vector of starting values for the coefficient parameters.}
  \item{sigma2}{
    Value of the dispersion parameter (variance).  The estimator does not update this
    value.}
  \item{resstart}{
    Array of initial orthogonalized (but not standardized) residuals.}
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
  \item{tvec}{
    Starting values for the regression coefficient parameters, as a matrix
    (parameters by alternatives).  Parameters that are involved in equality
    constraints are repeated in tvec.}
  \item{jacstack}{
    Array of regressors used to facilitate computing the gradient and the
    hessian matrix.
    dim(jacstack) = c(observations, unique parameters, alternatives).}
  \item{itmax}{
    Maximum number of Gauss-Newton stages.  Each stage does at most 100 Gauss-Newton
    steps.}
  \item{print.level}{
    Specify 0 for minimal printing (error messages only) or 2 to print details
    about the tanh computations.}
}
\details{
  The tanh estimator is a redescending M-estimator.  Given an estimate of the
  scale of the overdispersion, the tanh estimator estimates the coefficient
  parameters of the linear predictors of the multinomial regression model.
}
\value{
  mGNtanh returns a list of 16 objects.  The returned objects are:

  \item{coefficients}{
    The tanh coefficient estimates in matrix format.  The matrix has one
    column for each outcome alternative.  The label for each row of the matrix
    gives the names of the regressors to which the coefficient values in the row
    apply.  The regressor names in each label are separated by a forward
    slash (/), and \code{NA} is used to denote that no regressor is
    associated with the corresponding value in the matrix.  The value 0 is
    used in the matrix to fill in for values that do not correspond to a
    regressor.}
  \item{coeffvec}{
     A vector containing the tanh coefficient estimates.}
  \item{dispersion}{
    Value of the dispersion parameter (variance).  This is the value specified
    in the argument \code{sigma2} in the call to the function.}
  \item{w}{
    Vector of weights based on the tanh estimator's \code{psi} function for
    each observation.}
  \item{psi}{
    Vector of values of the tanh estimator's \code{psi} function for each observation.}
  \item{A}{
    The outer product of the gradient (expected information) divided by the
    moment estimate of the dispersion.}
  \item{B}{
    The inverse of the hessian matrix (observed formation).}
  \item{covmat}{
    Sandwich estimate of the asymptotic covariance of the tanh coefficient estimates.}
  \item{iters}{
    Number of Gauss-Newton iterations.}
  \item{error}{
    Error code:
    0,    no errors;
    2,   \code{sum(w) < nobs*(ncats-1)/2} (weights are too small);
    32,  Hessian not positive definite in the final Newton step.}
  \item{GNlist}{
    List reporting final results of the Gauss-Newton optimization.  Elements:
    \code{coefficients}, vector of coefficient parameters (same as \code{coeffvec}
    value in list returned by mGNtanh);  \code{tvec}, matrix of coefficient
    parameters (same as \code{coefficients} value in list returned by mGNtanh);
    \code{formation}, inverse Hessian matrix;  \code{score}, score (or gradient
    element) matrix;  \code{LLvals}, list containing weighted
    (\code{LLvals$LL}) and unweighted (\code{LLvals$LLu}) log-likelihood values;
    \code{convflag}, TRUE/FALSE convergence flag;  \code{iters}, number of
    iterations done in final Gauss-Newton stage;  \code{posdef}, TRUE if Hessian is
    positive definite.}
  \item{tanhsigma2}{
    The tanh overdispersion parameter estimate, which is a
    weighted moment estimate of the dispersion:  weighted mean sum of squared
    orthogonalized residuals (adjusted for effective sample size after
    weighting and degrees of freedom lost to estimated coefficients).}
  \item{Y}{
    The same \code{Y} matrix that was supplied as input, except modified by 
    having done \code{Y[!Ypos] <- 0}.}
  \item{Ypos}{
    The same \code{Ypos} matrix that was supplied as input.}
  \item{probmat}{
    The matrix of predicted probabilities for each category for each
    observation based on the coefficient estimates.}
  \item{jacstack}{
    The same \code{jacstack} that was supplied as an input argument.}
  \item{Xarray}{
    The same \code{Xarray} that was supplied as an input argument.}  
}
\references{
  Walter R. Mebane, Jr. and  Jasjeet Singh Sekhon. 2004.  ``Robust Estimation
  and Outlier Detection for Overdispersed Multinomial Models of Count Data.'' 
  \emph{American Journal of Political Science} 48 (April): 391--410.
  \url{http://sekhon.berkeley.edu/multinom.pdf}

  For additional documentation please visit \url{http://sekhon.berkeley.edu/robust/}.
}
\author{Walter R. Mebane, Jr., Cornell University,
  \email{wrm1@cornell.edu}, \url{http://macht.arts.cornell.edu/wrm1/} \cr

  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}}
\keyword{robust}
\keyword{models}
\keyword{regression}
\keyword{optimize}
