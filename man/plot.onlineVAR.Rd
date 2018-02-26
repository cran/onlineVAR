\name{plot.onlineVAR}
\alias{plot.onlineVAR}

\title{Coefficient matrix plots for onlineVAR objects.}

\description{
  Generates plots of estimated coefficient matrices for onlineVAR fits
  by using the \code{\link{levelplot}} function from \code{\link{lattice}}.
}
\usage{
\method{plot}{onlineVAR}(x, lag = 1, time = "last", s = NULL, col = NULL, 
  at = NULL, xlab = "", ylab = "", \dots)

}
\arguments{
  \item{x}{\code{"onlineVAR"} object, typically output from \code{onlineVAR}.}
  \item{lag}{Lag for which coefficient matrix is plotted.} 
  \item{time}{Time step from which coefficient matrix estimates are taken.
  Default is \code{"last"} to take the estimates from the last update.
  Other time steps can only be specified if model has been fitted with
  \code{trace=TRUE}.}
  \item{s}{Optional index specifying which lambda in the lambda.ratio sequence
  to use. If not specified, optimum s is taken based on weighted 
  squared error.}
  \item{col}{Color palette. Default is a color palette with reddish
  colors for negative and blueish colors for positive coefficients.}
  \item{at}{A numeric vector giving the breakpoints for the color palette.}
  \item{xlab}{A title for the x axis.}
  \item{ylab}{A title for the y axis.}
  \item{\dots}{Further arguments passed to \code{levelplot}.}
}



\seealso{\code{\link{onlineVAR}}, \code{\link{levelplot}}}

\examples{
data(aemo)

## use only subset of first 8000 time steps
y <- aemo[1:8000,]

## fit online lasso VAR
onlinefit <- onlineVAR(y, nu = 0.99, lags = 1, ahead = 1)

## plot coefficient matrix from last update
plot(onlinefit)

}
\keyword{vector time series regression}
