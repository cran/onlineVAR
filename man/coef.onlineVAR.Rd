\name{coef.onlineVAR}
\alias{coef.onlineVAR}
\alias{summary.onlineVAR}
\alias{print.summary.onlineVAR}
\alias{print.onlineVAR}

\title{Methods for onlineVAR Objects}

\description{
  Methods for extracting information from fitted \code{onlineVAR} objects.
}
\usage{
\method{coef}{onlineVAR}(object, time = "last", s = NULL, lags = NULL, \dots)
}
\arguments{
  \item{object}{An object of class \code{"onlineVAR"}.}
  \item{time}{Time step from which coefficient matrix estimates are taken.
  Default is "last" to take the estimates from the last update.
  Other time steps can only be specified if model has been fitted with
  \code{trace=TRUE}.}
  \item{s}{Optional index specifying which lambda in the lambda.ratio sequence
  to use. If not specified, optimum s is taken based on weighted 
  squared error.}
  \item{lags}{Optional vector specifying for which lags coefficients should
  be returned. If not specified, coefficient matrices for all lags are returned.}
  \item{\dots}{Required for S3 compatibility.}
}

\details{
  In addition to the \code{\link[stats]{coef}} method above, also the standard 
  extractor functions \code{\link[base]{print}} and
  \code{\link[base]{summary}} for \code{"onlineVAR"} objects are available. 
}


\seealso{\code{\link{onlineVAR}}}

\keyword{vector time series regression}
