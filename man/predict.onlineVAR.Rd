\name{predict.onlineVAR}
\alias{predict.onlineVAR}

\title{Predictions for onlineVAR Models}

\description{
  Obtains predictions for onlineVAR models. 
}
\usage{
\method{predict}{onlineVAR}(object, newdata, s = NULL, \dots)
}
\arguments{
  \item{object}{An object of class \code{"onlineVAR"}.}
  \item{newdata}{Optional time series data to predict. If supplied, predictions
  are derived with coefficients from the last time step 
  that was used for fitting \code{object}}
  \item{s}{Optional index specifying which lambda in the lambda.ratio sequence
  to use for prediction.}
  \item{\dots}{Required for S3 compatibility.}
}

\details{
  By default, predictions derived during the online fitting are returned.
  By specifying \code{newdata} predictions for a new time series
  are generated with coefficient estimates from the last time step of
  the online fitting.
}
\value{
  Time series matrix of predictions of the same size as the input data used for 
  model fitting or, if supplied, of the same size as \code{newdata}. 
}

\seealso{\code{\link{onlineVAR}}}

\keyword{vector time series regression}
