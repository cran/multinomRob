% $Id: multinomTanh.Rd,v 1.6 2004/02/18 07:07:40 wrm1 Exp $
%
\name{Multinomial Regression Tanh Estimator}
\alias{multinomTanh}
\title{Multinomial Regression Hyperbolic Tangent (Tanh) Estimator}
\description{
  \code{multinomTanh} fits the overdispersed multinomial regression
  model for grouped count data using the hyperbolic tangent (tanh)
  estimator.  This function is not meant to be called directly by the
  user.  It is called by \code{multinomRob}, which constructs the
  various arguments.
}
\usage{multinomTanh(Y, Ypos, X, jacstack, xvec, tvec, pop, s2,
                    xvar.labels, choice.labels, print.level = 0)
}
\arguments{
  \item{Y}{
    Matrix (observations by alternatives) of outcome counts.
    Values must be nonnegative.  Missing data (\code{NA} values) are not allowed.}
  \item{Ypos}{
    Matrix indicating which elements of Y are counts to be analyzed (TRUE)
    and which are values to be skipped (FALSE).  This allows the set of
    outcome alternatives to vary over observations.}
  \item{X}{
    Array of regressors.  dim(X) = c(observations, parameters, alternatives).}
  \item{jacstack}{
    Array of regressors used to facilitate computing the gradient and the
    hessian matrix.
    dim(jacstack) = c(observations, unique parameters, alternatives).}
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
  \item{pop}{
    Vector giving the total number of counts for each observation.  In general,
    \code{pop <- apply(Y * ifelse(Ypos,1,0), 1, sum)}.}
  \item{s2}{
    Overdispersion value.  In multinomRob this is the square of the LQD
    scale estimate.}
  \item{xvar.labels}{
    Vector of labels for observations.}
  \item{choice.labels}{
    Vector of labels for outcome alternatives.}
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
  multinomTanh returns a list of 5 objects.  The returned objects are:

  \item{mtanh}{
     List of tanh estimation results from function \code{mGNtanh}.}
  \item{weights}{
    The matrix of tanh weights for the orthogonalized residuals.  The matrix
    has the same dimensions as the outcome count matrix \code{Y}.  The first
    column of the matrix has names for the observations, and the remaining
    columns contain the weights.  Each of the latter columns has a name
    derived from the \code{choice.labels} vector:  column \code{i+1} is named
    \code{paste("weights:",choice.labels[i],sep="")}.

    If \code{sum(Ypos[i,]==FALSE)>0}, then values of \code{NA} appear in
    \code{weights[i,]}, with \code{sum(is.na(weights[i,]))==sum(!Ypos[i,])}.
    The \code{NA} values will be the last values in the affected
    row of the \code{weights} matrix, regardless of which outcome alternatives
    were unavailable for the observation.}
  \item{Hdiag}{
    The matrix of weights used to fully studentize the orthogonalized
    residuals.  The matrix has the same dimensions as the outcome count matrix
    \code{Y}.  The first column of the matrix has names for the observations,
    and the remaining columns contain the weights.  Each of the latter columns
    has a name derived from the \code{choice.labels} vector:  column \code{i+1}
    is named \code{paste("Hdiag:",choice.labels[i],sep="")}.

    If \code{sum(Ypos[i,]==FALSE)>0}, then values of 0 appear in
    \code{Hdiag[i,]}, with \code{sum(is.na(Hdiag[i,]))==sum(!Ypos[i,])}.
    The \code{0} values created for this reason will be the last values in the
    affected row of the \code{Hdiag} matrix, regardless of which outcome
    alternatives were unavailable for the observation.}
  \item{cr}{
    List of predicted outcome counts, studentized residuals and
    standardized residuals.}
  \item{tvec}{
    The tanh coefficient estimates in matrix format.  The matrix has one
    column for each outcome alternative.  The label for each row of the matrix
    gives the names of the regressors to which the coefficient values in the row
    apply.  The regressor names in each label are separated by a forward
    slash (/), and \code{NA} is used to denote that no regressor is
    associated with the corresponding value in the matrix.  The value 0 is
    used in the matrix to fill in for values that do not correspond to a
    regressor.}
}
\references{
  Walter R. Mebane, Jr. and  Jasjeet Singh Sekhon. 2004.  ``Robust Estimation
  and Outlier Detection for Overdispersed Multinomial Models of Count Data.'' 
  \emph{American Journal of Political Science} 48 (April): 391--410 .
  \url{http://macht.arts.cornell.edu/wrm1/multinom.pdf}

  For the most current code and related material see
  \url{http://jsekhon.fas.harvard.edu/robest/}
}
\author{Walter R. Mebane, Jr., Cornell University,
  \email{wrm1@cornell.edu}, \url{http://macht.arts.cornell.edu/wrm1/} \cr

  Jasjeet S. Sekhon, Harvard University, \email{jasjeet\_sekhon@harvard.edu},
  \url{http://jsekhon.fas.harvard.edu/}
}
\keyword{robust}
\keyword{models}
\keyword{regression}
