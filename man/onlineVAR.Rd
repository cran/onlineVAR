\name{onlineVAR}
\alias{onlineVAR}
\alias{onlineVAR.fit}

\title{Online Fitting of Time-Adaptive Lasso Vector Auto Regression.}

\description{
  Recursive fitting of time-adaptive lasso vector auto regressive models with 
  an online cyclic coordinate descent algorithm. 
}
\usage{
onlineVAR(y, nu = 1, lags = NULL, ahead = NULL, lambda.ratio = NULL, 
  burnin = NULL, control = onlineVAR.control(lambda.ratio, \dots), \dots)

onlineVAR.fit(x, y, nu, lags, ahead, lambda.ratio, fit, burnin,
  control = onlineVAR.control())
}
\arguments{
  \item{y}{Vector time series with columns for each variable and rows 
  for each time step.}
  \item{x}{Lagged time series data. Columns 1 to \code{ncol(y)} contain lag-1
  data, columns \code{ncol(y)+1} to \code{2*ncol(y)} lag-2 data and so on.}
  \item{nu}{Forgetting factor for exponential forgetting.} 
  \item{lags}{Number of considered lags. Order of VAR model.}
  \item{ahead}{Look ahead time for predictions. Most recent data that are
  available. The fitted model only considers lags larger than \code{ahead}.}
  \item{lambda.ratio}{Vector of penalization parameters as fractions of the
  minimum lambda for which all coefficients are zero. If not specified
  a sequence of lambda values is used. See \code{\link{onlineVAR.control}}
  for details.}
  \item{burnin}{Burn-in period, i.e., the number of time steps that are used to 
  stabilize the model fit. This is particularly important for selecting the
  optimum lambda. Default is the effective training data length \code{1/(1-nu)}
  for a vector of \code{lambda.ratio} and 0 for a single valued 
  \code{lambda.ratio}.}
  \item{control}{A list of control parameters. Default is 
  \code{\link{onlineVAR.control}(lambda.ratio)}.}
  \item{fit}{List of various starting values that are used for the online updating
  algorithm.}
  \item{\dots}{Arguments to be used to form the default \code{control} argument
  if it is not supplied directly.}
}

\details{
\code{onlineVAR} recursively fits time-adaptive lasso vector auto regressive 
models with an online coordinate descent algorithm that updates the model fit 
at each time step. By default these models are fitted for a sequence of lasso 
penalization parameters and predictions are generated for each time step with 
the parameter with the smallest weighted mean squared error. 


\code{onlineVAR.fit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of class \code{"onlineVAR"}, i.e., a list with the 
    following elements.
  \item{coef}{List of intercept vector and coefficient matrices for the
  different lags.}
  \item{fit}{List of various values that are used for the online updating
  algorithm. This is mainly used as starting values for further updates.}
  \item{nu}{Forgetting factor for exponential forgetting.}
  \item{pred}{Predictions. Matrix of the same size as the input data but with
  predictions from adaptive lasso VAR. The first \code{burnin+ahead+lags} 
  rows are filled with \code{NA}}
  \item{predall}{If \code{predall=TRUE} in \code{control} predictions for all 
  penalization parameters in \code{lambda.ratio}. Else \code{NULL}}
  \item{residuals}{Prediction errors.}
  \item{lags}{Number of considered lags. Order of VAR model.}
  \item{ahead}{Look ahead time for predictions.}
  \item{lambda.ratio}{Sequence of lasso penalization parameters}
  \item{trace}{If \code{trace=TRUE} a list of coefficient estimates for each
  time step. Else \code{FALSE}}
  \item{burnin}{Burn-in period, i.e., the number of data that are used to 
  stabilize the model fit.}
  \item{call}{Function call.}
}

\references{
  Messner JW, Pinson P (2016). Online Adaptive LASSO Estimation in Vector Auto Regressive Models for Wind Power Forecasting in High Dimension.
    \emph{International Journal of Forecasting}, in press.
    \url{http://pierrepinson.com/docs/MessnerPinson18.pdf}.
}

\seealso{\code{\link{predict.onlineVAR}}, \code{\link{plot.onlineVAR}}}

\examples{
data(aemo)

## use only subset of first 8000 time steps
y <- aemo[1:8000,]

## fit online lasso VAR
onlinefit <- onlineVAR(y, nu = 0.99, lags = 1, ahead = 1)

## plot coefficient matrix from last update
plot(onlinefit)

## compare mean root mean squared error to persistence 
c(onlinefit   = mean(sqrt(apply((predict(onlinefit)-y)^2, 2, 
    mean, na.rm = TRUE))),
  persistence = mean(sqrt(apply((aemo[1000:7999,]-y[1001:8000,])^2, 2, 
    mean))))
}
\keyword{vector time series regression}
