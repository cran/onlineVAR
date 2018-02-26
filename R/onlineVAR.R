
###############################################################################
## Online VAR control function ################################################
###############################################################################

onlineVAR.control <- function(lambda.ratio = NULL, nlambda = NULL, 
  lambda.min.ratio = NULL, abstol = 0.001, trace = FALSE, start = NULL, 
  parallel = FALSE, predall = FALSE) 
{
  if(is.null(start)) {
    if(!is.null(lambda.ratio)) {
      if(!is.null(nlambda) | !is.null(lambda.min.ratio)) warning("if lambda is supplied nlambda and lambda.min.ratio are ignored")
      nlambda <- length(lambda.ratio)
      lambda.min.ratio <- tail(lambda.ratio, 1)
    } else {
      if(is.null(nlambda)) nlambda <- 10
      if(is.null(lambda.min.ratio)) lambda.min.ratio <- 0.0001
      ## lambda sequence relative to lambda max (minimum lambda such that all coef are 0)
      lambda.ratio <- exp(-log(lambda.min.ratio)/(nlambda-1)*(nlambda-c(1:nlambda))) * 
        lambda.min.ratio
    }     
  } else {
    if(!is.null(nlambda) | !is.null(lambda.min.ratio) | !is.null(lambda.ratio)) warning("nlambda, lambda.min.ratio, lambda taken from start. Input is ignored")
    lambda.ratio <- start$lambda.ratio    
    nlambda <- length(lambda.ratio)
    lambda.min.ratio <- tail(lambda.ratio, 1)
    
  }
  rval <- list(nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
    abstol = abstol, lambda.ratio = lambda.ratio, trace = trace, start = start, 
    parallel = parallel, predall = predall)
  rval
}



###############################################################################
## Online VAR function ########################################################
###############################################################################

onlineVAR <- function(y, nu = 1, lags = NULL, ahead = NULL, lambda.ratio = NULL,
  burnin = NULL, control = onlineVAR.control(lambda.ratio, ...), ...) {
  
  cl <- match.call()

  N <- NROW(y)
  ## set default burnin period
  start <- control$start
  if(is.null(burnin)) {
    if(is.null(start)) {
      if(control$nlambda == 1) {
        burnin <- 0
      } else {
        burnin <- 1/(1-nu)
        if(burnin >= N) {
          warning("Effective training data length shorter than data length. burnin set to 0. Forecasts might not be optimal.")
          burnin <- 0
        }
      }
    } else {
      burnin <- 0
    }
  }

  ## set default lags and ahead
  if(is.null(start)) {
    if(is.null(lags)) lags <- 1
    if(is.null(ahead)) ahead <- 1
  } else {
    if(!is.null(lags) | !is.null(ahead)) warning("lags and ahead taken from start. Input is ignored")
    lags <- start$lags
    ahead <- start$ahead
  }

  Q <- NCOL(y)
  P <- Q*lags

  ## create starting values 
  if(is.null(start)) {
    coef2 <- matrix(0, ncol = P, nrow = P/lags)
    if(lags > 1) coef2 <- extendmatrix(coef2, lags)
    zeromatrix <- matrix(0, ncol = P, nrow = P)
    coef <- werror <- R <- r <- list()
    for(i in 1:control$nlambda) {
      coef[[i]] <- coef2
      R[[i]] <- zeromatrix
      r[[i]] <- zeromatrix
    }
    fit <- list(coef = coef, R = R, r = r, mu = rep(0, P), wN = 0, 
      werror = rep(0, control$nlambda))
    start <- list(fit = fit)
  } 

  ## create matrix y such that model can be fitted as VAR1
  names <- colnames(y)
  y <- duplicatedata(y, lags)
  names2 <- c(rownames(y), rep("NA", ahead))
  y <- na.omit(y)
  x <- head(y, -ahead)
  y <- tail(y, -ahead)

  ## call the actual workhorse lassovar.fit()
  rval <- onlineVAR.fit(x = x, y = y, nu = nu, lags = lags, ahead = ahead, 
    lambda.ratio = control$lambda.ratio, fit = start$fit, burnin = burnin,
    control = control)
    
  ## fill up prediction matrix
  ## NA predictions for the first lags+ahead rows
  rval$pred <- rbind(matrix(NA, nrow = lags-1 + ahead, ncol = Q), rval$pred)
  rval$residuals <- rbind(matrix(NA, nrow = lags-1 + ahead, ncol = Q), rval$residuals)
  ## row and column names
  rownames(rval$pred) <- names2[1:nrow(rval$pred)]
  colnames(rval$pred) <- colnames(y)[1:Q]
  rval$call <- cl

  return(rval)
}


###############################################################################
## Online VAR fitting function ################################################
###############################################################################

onlineVAR.fit <- function(x, y, nu, lags, ahead, lambda.ratio, fit, burnin, 
  control = onlineVAR.control()) {

  abstol <- control$abstol
  trace <- control$trace
  parallel <- control$parallel
  nlambda <- control$nlambda
  seq <- lambda.ratio

  N <- nrow(y)
  P <- ncol(y)
  Q <- P/lags

  pred <- coefall <- intercept <- pred2 <- list()
  
  ## Call C function for multiple penalization parameters
  fitfun <- function(s) {
    fit2 <- .Call("onlineVARfit", x = x, y = y, nu = nu, mui = fit$mu,  
      wNi=fit$wN, Ri = fit$R[[s]], ri = fit$r[[s]], seq = seq[s], coefi = fit$coef[[s]], 
      q = as.integer(Q), abstol = abstol, trace = as.integer(trace))
    fit2
  }

  fit2 <- if(parallel) mclapply(1:length(seq), fitfun) else 
    lapply(1:nlambda, fitfun)

  ## derive errors for different lambdas
  error <- lapply(1:length(seq), 
    function(s) c(fit$werror[s], rowSums((fit2[[s]][[6]] - y[,1:Q])^2)))

  werrorfun <- function(error) {
    werror0 <- error[1]
    error <- error[-1]
    rval <- rep(0, N)
    rval[burnin+1] <- nu*werror0 + error[burnin+1] 
    for(n in (burnin+2):N) rval[n] <- nu*rval[n-1] + error[n]
    rval
  }
  werror <- sapply(error, werrorfun)

  ## store predictions of best performing lambda
  pred <- matrix(NA, nrow = N, ncol = Q)
  for(n in (burnin+1):N) {
    pred[n,] <- fit2[[which.min(werror[n,])]][[6]][n,]
  }

  ## store list of predictions for all lambdas
  if(control$predall) {
    predall <- list()
    for(s in 1:length(seq)) predall[[s]] <- fit2[[s]][[6]]
  } else {
    predall <- NULL
  }


  ## if trace=TRUE: store coefficients for all time steps
  if(trace) {
    mu <- fit$mu
    wN <- fit$wN
    coefall <- list()
    for(n in 1:(N-ahead)+1) {
      mu <- nu*mu + y[n,]
      wN <- nu*wN + 1    
      coefn <- fit2[[which.min(werror[n-1,])]][[7]][((n-1)*P+1):(n*P),]
      interceptn <- (mu - coefn%*%mu)/wN
      coefall[[n+ahead]] <- list(intercept = interceptn[1:Q])
      for(l in 1:lags) coefall[[n+ahead]][[paste0("lag", l)]] <- coefn[1:Q, ((l-1)*Q+1):(l*Q)]
    }
    trace <- coefall
  }

  ## store additional informations required for further updates 
  coef <- R <- r <- list()
  for(i in 1:nlambda) {
    coef[[i]] <- fit2[[i]][[5]]
    R[[i]] <- fit2[[i]][[3]]
    r[[i]] <- fit2[[i]][[4]]
  }
  mu <- fit2[[1]][[1]]
  wN <- fit2[[1]][[2]]

  fit <- list(coef = coef, R = R, r = r, mu = mu, wN = wN, 
    werror = werror[N-ahead,])

  s <- which.min(werror[N,])

  ## coefficients from last time step for best performing lambda
  intercept <- (mu - coef[[s]]%*%mu)/wN
  coef2 <- list(intercept = intercept[1:Q])
  for(l in 1 : lags) coef2[[paste0("lag", l)]] <- coef[[s]][1:Q, ((l-1)*Q+1):(l*Q)]

  
  rval <- list(coef = coef2, 
    fit = fit, 
    nu = nu,
    pred = pred,
    predall = predall,
    residuals = y[,1:Q] - pred,
    lags = lags,
    ahead = ahead,
    lambda.ratio = seq,
    trace = trace,
    burnin = burnin)

  class(rval) <- "onlineVAR"
  return(rval)

}





###############################################################################
## Extractor functions ########################################################
###############################################################################
print.onlineVAR <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  cat(paste("\nCoefficients after", NROW(x$pred), "time steps:\n"))
  cat("\nIntercept:\n")
  print.default(format(x$coef[[1]], digits = digits), print.gap = 2, quote = FALSE)
  for(i in 2:length(x$coef)){
    cat(paste("\nLag-", i-1, " coefficients:\n", sep = ""))
    print.default(format(x$coef[[i]], digits = digits), print.gap = 2, quote = FALSE)
  }
  invisible(x)
}


summary.onlineVAR <- function(object, ...)
{
  ## return
  class(object) <- "summary.onlineVAR"
  object
}


print.summary.onlineVAR <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  cat("Residuals:\n")

  print(summary(x$residuals[-c(1:(x$burnin+x$lags+x$ahead-2)),], digits = digits))
  
  cat("\nCovariance matrix of residuals:\n")

  print.default(format(cov(x$residuals, use = "na.or.complete"), digits = digits), print.gap = 2, quote = FALSE)

  cat(paste("\nCoefficients after", NROW(x$pred), "time steps:\n"))
  cat("\nIntercept:\n")
  print.default(format(x$coef[[1]], digits = digits), print.gap = 2, quote = FALSE)
  for(i in 2:length(x$coef)){
    cat(paste("\nLag-", i-1, " coefficients:\n", sep = ""))
    print.default(format(x$coef[[i]], digits = digits), print.gap = 2, quote = FALSE)
  }  
  invisible(x)
}



predict.onlineVAR <- function(object, newdata=NULL, s = NULL, ...) {
  if(is.null(newdata)) {
    if(is.null(s)) {
      rval <- object$pred
    } else {
      if(is.null(object$predall)) stop("predictions for user defined s cannot be derived when predall = FALSE")
      rval <- object$predall[[s]]
    }
  } else {
    if(is.null(s)) s <- which.min(object$fit$werror)
    wN <- object$fit$wN
    std <- object$fit$sigma/wN
    mu <- object$fit$mu
    coef <- object$fit$coef
    coef2 <-  coef[[s]]
    intercept <- (mu - coef2%*%mu)/wN
    Q <- (ncol(coef2)/object$lags)
    y <- duplicatedata(newdata, object$lags, ahead = object$ahead)
    x <- head(y[, (Q + 1):ncol(y)], -(object$ahead+object$lags - 1))
  
    rval <- matrix(rep(intercept, nrow(x)), nrow = nrow(x), byrow = TRUE) + t(coef2 %*% t(x))
    rval <- rval[,1:Q]
    colnames(rval) <- colnames(y)[1:Q]
  }
  return(rval)
}


coef.onlineVAR <- function(object, time = "last", s = NULL, lags = NULL, ...) {
  if(is.null(s)) {
    s <- which.min(object$fit$werror)
    if(time == "last") {
      coef3 <- object$coef
    } else {
      if(length(object$trace)==1) stop("models fit with trace=FALSE only support time = 'last'")
      coef3 <- object$trace[[time]]
    }
  } else {
    if(time != "last") stop("if s not equal to NULL time has to be 'last'")
    Q <- ncol(object$pred)
    lags <- object$lags
    wN <- object$fit$wN
    std <- object$fit$sigma/wN
    mu <- object$fit$mu
    coef <- object$fit$coef
    coef2 <-  coef[[s]]
    intercept <- (mu - coef2%*%mu)/wN

    coef3 <- list(intercept = intercept[1:Q])
    for(l in 1 : lags) coef3[[paste0("lag", l)]] <- coef2[1:Q, ((l-1)*Q+1):(l*Q)]
  }

  if(lags == "all") {
    rval <- coef3
  } else {
    rval <- coef3[lags]
  }
  return(rval)
}




plot.onlineVAR <- function(x, lag = 1, time = "last", s = NULL, col = NULL, 
  at = NULL, xlab = "", ylab = "",  ...)  {
  if(is.null(col)) col <- c("#8E063B", "#941E44", "#9B2D4D", "#A13B56",     
    "#A84960", "#B0576B", "#B86577", "#C07585", "#CA8894", "#D49DA7",
    "#E2B9C0", "#FFFFFF", "#BDC2DF", "#A3AAD1", "#8E97C7", "#7D88BF", 
    "#6E7BB8", "#5F6FB3", "#5264AE", "#445AA9", "#3651A6", "#2448A5",     
    "#023FA5")
  ## color palette can also be created with diverge_hcl from colorspace
  ## rev(diverge_hcl(23, l = c(30, 100),power = 0.5))

  if(is.null(s)) {
    s <- which.min(x$fit$werror)
    if(time == "last") {
      coef3 <- x$coef
    } else {
      if(length(x$trace)==1) stop("time can only be 'last' for models fit with trace=FALSE")
      coef3 <- x$trace[[time]]
    }
  } else {
    if(time != "last") stop("time can only be 'last' when s not equal to NULL")
    Q <- ncol(x$pred)
    lags <- x$lags
    wN <- x$fit$wN
    std <- x$fit$sigma/wN
    mu <- x$fit$mu
    coef <- x$fit$coef
    coef2 <-  coef[[s]]
    intercept <- (mu - coef2%*%mu)/wN

    coef3 <- list(intercept = intercept[1:Q])
    for(l in 1 : lags) coef3[[paste0("lag", l)]] <- coef2[1:Q, ((l-1)*Q+1):(l*Q)]
  }
  coef <- coef3[[lag+1]]

  if(is.null(at)) {
    limits <- c(min(sapply(coef, min)[-1]), max(sapply(coef, max)[-1]))
    maxlim <- max(abs(limits))
    at <- seq(maxlim/11, maxlim, length.out= 11)
    at <- c(-rev(at), -maxlim/100, maxlim/100, at)
  }
  if(abs(limits[1])<abs(limits[2])) {
    col <- col[at[-1]>=limits[1]]
    at <- at[tail(which(at<limits[1]), 1):length(at)]
  } else {
    col <- col[at[-24]<=limits[1]]
    at <- at[1:which(at>limits[1])[1]]
  }

  

  if(is.null(rownames(coef))) rownames(coef) <- 1:nrow(coef)
  if(is.null(colnames(coef))) colnames(coef) <- 1:ncol(coef)
  z <- t(coef)[, ncol(coef):1]
  levelplot(z, col.regions=col, at = at, 
    xlab = xlab, ylab = ylab,...)  
}






###############################################################################
## additional functions required for other functions ##########################
###############################################################################


## function to get response matrix
duplicatedata <- function(y, lags, ahead = 0) {
  
  y2 <- y3 <- as.matrix(y)
  names <- rownames(y2)
  if(is.null(names)) names <- 1:nrow(y2)
  if(ahead) {
    for(i in 1:ahead) {
      y2 <- rbind(y2, NA)
      y3 <- rbind(NA, y3)
      names <- c(names, "NA")
    }
    y2 <- cbind(y2, y3)
  }

  if(lags>1) {
    for(i in 1:(lags-1)) {
      y3 <- rbind(NA, y3)
      y2 <- cbind(rbind(y2, NA), y3)
      names <- c(names, "NA")
    }
  } 
  rownames(y2) <- names
  
  return(y2)
}

## function to extend coefficientmatrix to transform VAR-p into VAR-1 model
extendmatrix <- function(coef, lags) {
  P <- ncol(coef)/lags
  coef <- rbind(coef, matrix(0, ncol = P*lags, nrow = P*(lags-1)))
  for(j in 1:(lags-1)) {
    coef[(j*P) + 1:P, ((j-1)*P) + 1:P] <- diag(P)
  }
  coef
}



 
