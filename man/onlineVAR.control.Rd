\name{onlineVAR.control}
\alias{onlineVAR.control}

\title{Auxiliary Function for Controlling onlineVAR Fitting}

\description{
  Auxiliary function for \code{onlineVAR} fitting. 
}
\usage{
onlineVAR.control(lambda.ratio = NULL, nlambda = NULL, 
  lambda.min.ratio = NULL, abstol = 0.001, trace = FALSE, start = NULL, 
  parallel = FALSE, predall = FALSE)
}
\arguments{
  \item{lambda.ratio}{Vector of penalization parameters as fractions of the
  minimum lambda for which all coefficients are zero. If not specified
  a sequence of lambda values is generated based on \code{nlambda} and 
  \code{lambda.min.ratio}. If supplied, \code{nlambda} and 
  \code{lambda.min.ratio} are ignored.}
  \item{nlambda}{Number of lasso penalization parameters lambda. 
  Default is 10.}
  \item{lambda.min.ratio}{Smallest value of lambda.ratio. Default is 0.0001}
  \item{abstol}{Absolute tolerance for coordinate descent convergence.
  In each time step the algorithm stops when the sum of coefficient estimates
  does not change more than \code{abstol}. Default is 0.001.} 
  \item{trace}{If \code{TRUE} coefficient estimates are stored for all time 
  steps. If \code{FALSE} coefficient matrices are only stored for the last time
  step to save memory.}
  \item{start}{Object of class \code{"onlineVAR"}. Coefficient estimates from the 
  last time step are used as staring values. Can be used to continue updating
  the model with new data.}
  \item{parallel}{If \code{TRUE} the model fitting for the different lambda
  is parallelized.}
  \item{predall}{Logical whether predictions from all penalization parameters
  in the sequence are stored.}
}

\value{
  An list of components named as the arguments.
  \item{nlambda}{Number of lasso penalization parameters in the lambda
  sequence.}
  \item{lambda.min.ratio}{Smallest value for lambda.ratio.}
  \item{abs.tol}{Absolute tolerance for coordinate descent convergence.}
  \item{lambda.ratio}{Lambda sequence as fractions of the
  minimum lambda for which all coefficients are zero.}
  \item{trace}{Logical whether coefficients should be stored for all time 
  steps.}
  \item{start}{Starting values.}
  \item{parallel}{Logical whether the model fitting for the different lambda
  is parallelized.}
  \item{predall}{Logical whether prediction from all penalization parameters
  in the sequence are stored.}
}


\seealso{\code{\link{onlineVAR}}}


\keyword{vector time series regression}
