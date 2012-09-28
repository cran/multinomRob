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
                   choice.labels, MLEonly=FALSE, print.level=0)
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
  \item{MLEonly}{
    If \code{TRUE}, then only the standard maximum-likelihood MNL model
    is estimated---i.e., no overdispersion parameter is estimated.}  
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
  multinomMLE returns a list containing the following objects.  The
  returned objects are:
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
     The MLE coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the Hessian matrix (observed
     information).}
   \item{se.opg}{
     The MLE coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the outer product of the gradient (expected
     information) divided by the moment estimate of the dispersion.
     Not provided if \code{MLEonly==TRUE}.}
   \item{se.hes}{
     The MLE coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the Hessian matrix (observed
     information).  Same as \code{se}; included for backward
     compatibility.}
   \item{se.sw}{
     The MLE coefficient estimate standard errors derived from the asymptotic
     covariance estimated using the estimated asymptotic
     sandwich covariance estimate. Not provided if \code{MLEonly==TRUE}.}   
   \item{se.vec}{
     \code{se} in vector form.}   
   \item{se.opg.vec}{
     \code{se.opg} in vector form.}
   \item{se.hes.vec}{
     \code{se.hes} in vector form.}      
   \item{se.sw.vec}{
     \code{se.sw} in vector form.}   
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
  \emph{American Journal of Political Science} 48 (April): 391--410.
  \url{http://sekhon.berkeley.edu/multinom.pdf}

  For additional documentation please visit \url{http://sekhon.berkeley.edu/robust/}.
}
\author{Walter R. Mebane, Jr., University of Michigan,
  \email{wmebane@umich.edu}, \url{http://www-personal.umich.edu/~wmebane} \cr

  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}
}
\examples{
# make some multinomial data
x1 <- rnorm(50);
x2 <- rnorm(50);
p1 <- exp(x1)/(1+exp(x1)+exp(x2));
p2 <- exp(x2)/(1+exp(x1)+exp(x2));
p3 <- 1 - (p1 + p2);
y <- matrix(0, 50, 3);
for (i in 1:50) {
  y[i,] <- rmultinomial(1000, c(p1[i], p2[i], p3[i]));
}

# perturb the first 5 observations
y[1:5,c(1,2,3)] <- y[1:5,c(3,1,2)];
y1 <- y[,1];
y2 <- y[,2];
y3 <- y[,3];

# put data into a dataframe
dtf <- data.frame(x1, x2, y1, y2, y3);

#Do MLE estimation.  The following model is NOT identified if we
#try to estimate the overdispersed MNL.
dtf <- data.frame(y1=c(1,1),y2=c(2,1),y3=c(1,2),x=c(0,1))
summary(multinomRob(list(y1 ~ 0, y2 ~ x, y3 ~ x), data=dtf, MLEonly=TRUE))
}
\keyword{models}
\keyword{regression}
