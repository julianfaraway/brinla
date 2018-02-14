#' linear regression model Baysian residual plot
#'
#' @author Xiaofeng Wang, \email{wangx6@ccf.org}
#' @param inla.obj an object of class "inla". 
#' @param covariate the covariate of interest for the residual plot.
#' @param m the size of posterior samples to be generated from the posterior distribution.
#' @param pmedian if pmedian = TRUE, posterior medians of Bayesian residuals will be calculated.
#' @param smooth if smooth = TRUE, a nonparametric smooth curve will be added.
#' @param xlab a title for the x axis.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param ... 
#' 
#' @details generate a Bayesian residual index plot or a plot of the residual versus the predictor for linear regression using INLA.
#' @return a list consisting of residuals (resid), the covariate used for plot (covariate), the response variable (y)
#' @export
bri.lmresid.plot <- function(inla.obj, covariate = NULL, m = 1000, pmedian = FALSE, smooth = FALSE, xlab = NULL, cex.lab = 1.25, cex.axis = 1.25, ylab = "Bayesian residual", main = "",...){	
  if(class(inla.obj) != "inla") stop("an 'inla' object is needed!")
  if(inla.obj$.args$family != "gaussian") stop("the function only supports Gaussian linear regression model!")
  # Get a single random draw of posterior distribution of parameters
  p1 <- length(names(inla.obj$marginals.fixed))
  beta <- matrix(0, nrow = m, ncol = p1)
  for (j in 1:p1){
    beta[, j] <- inla.rmarginal(m, marg = inla.obj$marginals.fixed[[j]])
  }
  yhat <- inla.obj$model.matrix %*% t(beta)
  
  org.formula <- inla.obj$.args$formula
  yname <- as.character(org.formula)[2]
  y <- inla.obj$.args$data[, yname]
  ymat <- matrix(rep(y, times = m), ncol = m)
  resid.sample <- ymat - yhat
  if (pmedian) {
    resid <- apply(resid.sample, 1, median)
  } else {
    resid <- apply(resid.sample, 1, mean)
  }
  
  if (is.null(covariate)) {
    covariate <- seq(1: length(resid))
    xlab <- "Index"		
  } else {
    if (length(resid) != length(covariate))
      stop("the 'covariate' variable does not match with the residual object!")		
    if (is.null(xlab)) xlab <- "Covariate"
  }
  if (smooth){
    res.dat <- data.frame(covariate = covariate, resid = resid)
    res.smooth.inla <- inla(resid ~ -1 + f(covariate, model = 'rw2', constr = FALSE), data = res.dat)
    bri.band.plot(res.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, cex.lab = cex.lab, cex.axis = cex.axis, type = 'random', ylim= c(min(resid), max(resid)))
    points(res.dat$covariate, res.dat$resid)		
  } else{
    plot(covariate, resid, cex.lab = cex.lab, cex.axis = cex.axis, xlab = xlab, main = main, ylab = ylab, ylim= c(min(resid), max(resid)), ...)
  }
  return(invisible(list(resid = resid, covariate = covariate, y = y))) 	
}


#' Baysian Pearson residuals for Poisson regression using INLA
#'
#' @author Xiaofeng Wang, \email{wangx6@ccf.org}
#' @param inla.obj an object of class "inla". 
#' @param covariate the covariate of interest for the residual plot.
#' @param m the size of posterior samples to be generated from the posterior distribution.
#' @param pmedian if pmedian = TRUE, posterior medians of Bayesian residuals will be calculated.
#' @param plot if plot = TRUE, a Bayesian residual plot will be generated
#' @param smooth if smooth = TRUE, a nonparametric smooth curve will be added.
#' @param xlab a title for the x axis.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param ... 
#' @details compute Bayesian Pearson residuals, and generate a index plot or a plot of the residual versus the predictor for Bayesian Possion regression using INLA.
#' @return a list consisting of residuals (resid), the covariate used for plot (covariate), the response variable (y)
#' @export
bri.Pois.resid <- function(inla.obj, covariate = NULL, m = 1000, pmedian = FALSE, plot = FALSE, smooth = FALSE, xlab = NULL, cex.lab = 1.25, cex.axis = 1.25, ylab = "Bayesian residual", ylim = NULL, main = "",...){	
  if(class(inla.obj) != "inla") stop("an 'inla' object is needed!")
  if(inla.obj$.args$family != "poisson") stop("the function only supports Poisson GLM!")
  # Get a single random draw of posterior distribution of parameters
  p1 <- length(names(inla.obj$marginals.fixed))
  beta <- matrix(0, nrow = m, ncol = p1)
  for (j in 1:p1){
    beta[, j] <- inla.rmarginal(m, marg = inla.obj$marginals.fixed[[j]])
  }
  yhat <- exp(inla.obj$model.matrix %*% t(beta))
  
  org.formula <- inla.obj$.args$formula
  yname <- as.character(org.formula)[2]
  y <- inla.obj$.args$data[, yname]
  ymat <- matrix(rep(y, times = m), ncol = m)
  if (is.null(inla.obj$.args$E)) {
    resid.sample <- (ymat - yhat)/sqrt(yhat)
  } else {
    Emat <- matrix(rep(inla.obj$.args$E, times = m), ncol = m)
    yhat.E <- yhat*Emat
    resid.sample <- (ymat - yhat.E)/sqrt(yhat.E)
  }
  
  if (pmedian) {
    resid <- apply(resid.sample, 1, median)
  } else {
    resid <- apply(resid.sample, 1, mean)
  }
  if (plot) {
    if (is.null(covariate)) {
      covariate <- seq(1: length(resid))
      xlab <- "Index"		
    } else {
      if (length(resid) != length(covariate))
        stop("the 'covariate' variable does not match with the residual object!")		
      if (is.null(xlab)) xlab <- "Covariate"
    }
    if (is.null(ylim)) {ylim <- c(min(resid), max(resid))}
    if (smooth){
      res.dat <- data.frame(covariate = covariate, resid = resid)
      res.smooth.inla <- inla(resid ~ -1 + f(covariate, model = 'rw2', constr = FALSE), data = res.dat)
      bri.band.plot(res.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, cex.lab = cex.lab, cex.axis = cex.axis, type = 'random', ylim = ylim)
      points(res.dat$covariate, res.dat$resid)		
    } else{
      plot(covariate, resid, cex.lab = cex.lab, cex.axis = cex.axis, xlab = xlab, main = main, ylab = ylab, ylim = ylim, ...)
    }
  } 	
  return(invisible(list(resid = resid, covariate = covariate, y = y))) 	
}

#' Baysian residuals for beta regression using INLA
#'
#' @author Xiaofeng Wang, \email{wangx6@ccf.org}
#' @param inla.obj an object of class "inla". 
#' @param covariate the covariate of interest for the residual plot.
#' @param m the size of posterior samples to be generated from the posterior distribution.
#' @param pmedian if pmedian = TRUE, posterior medians of Bayesian residuals will be calculated.
#' @param plot if plot = TRUE, a Bayesian residual plot will be generated
#' @param smooth if smooth = TRUE, a nonparametric smooth curve will be added.
#' @param xlab a title for the x axis.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param ... 
#' @details compute Bayesian residuals, and generate a index plot or a plot of the residual versus the predictor for Bayesian beta regression using INLA.
#' @return a list consisting of residuals (resid), the covariate used for plot (covariate), the response variable (y)
#' @export
bri.beta.resid <- function(inla.obj, covariate = NULL, m = 1000, pmedian = FALSE, plot = FALSE, smooth = FALSE, xlab = NULL, cex.lab = 1.25, cex.axis = 1.25, ylab = "Bayesian residual", ylim = NULL, main = "",...){	
  if(class(inla.obj) != "inla") stop("an 'inla' object is needed!")
  if(inla.obj$.args$family != "beta") stop("the function only supports beta GLM!")
  # Get a single random draw of posterior distribution of parameters
  p1 <- length(names(inla.obj$marginals.fixed))
  beta <- matrix(0, nrow = m, ncol = p1)
  for (j in 1:p1){
    beta[, j] <- inla.rmarginal(m, marg = inla.obj$marginals.fixed[[j]])
  }
  yhat <- exp(inla.obj$model.matrix %*% t(beta))/(1+exp(inla.obj$model.matrix %*% t(beta)))
  phihat <- inla.rmarginal(m, marg = inla.obj$marginals.hyperpar[[1]])
  phihatmat <- t(matrix(rep(phihat, times = nrow(yhat)), ncol = nrow(yhat)))
  varhat <- yhat*(1-yhat)/(1+phihatmat)
  
  org.formula <- inla.obj$.args$formula
  yname <- as.character(org.formula)[2]
  y <- inla.obj$.args$data[, yname]
  ymat <- matrix(rep(y, times = m), ncol = m)
  resid.sample <- (ymat - yhat)/sqrt(varhat)
  
  if (pmedian) {
    resid <- apply(resid.sample, 1, median)
  } else {
    resid <- apply(resid.sample, 1, mean)
  }
  
  if (plot) {
    if (is.null(covariate)) {
      covariate <- seq(1: length(resid))
      xlab <- "Index"		
    } else {
      if (length(resid) != length(covariate))
        stop("the 'covariate' variable does not match with the residual object!")		
      if (is.null(xlab)) xlab <- "Covariate"
    }
    if (is.null(ylim)) {ylim <- c(min(resid), max(resid))}
    if (smooth){
      res.dat <- data.frame(covariate = covariate, resid = resid)
      res.smooth.inla <- inla(resid ~ -1 + f(covariate, model = 'rw2', constr = FALSE), data = res.dat)
      bri.band.plot(res.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, cex.lab = cex.lab, cex.axis = cex.axis, type = 'random', ylim= ylim)
      points(res.dat$covariate, res.dat$resid)		
    } else{
      plot(covariate, resid, cex.lab = cex.lab, cex.axis = cex.axis, xlab = xlab, main = main, ylab = ylab, ylim = ylim, ...)
    }
  } 	
  return(invisible(list(resid = resid, covariate = covariate, y = y))) 	
}



#' Plot the baseline hazard functions for survival models using INLA
#' @param inla.obj an object of class "inla". 
#' @param plot if plot = TRUE, the plot of the baseline hazard functions will be generated.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ... 
#'
#' @details The function supports Weibull model, Expontenial model and Cox proportional hazards model for right censored data.
#' @return A list of evaluated time points (time), and the values of corresponding base line function (basehaz).
#' @export
bri.basehaz.plot <- function(inla.obj, plot = TRUE, cex.lab = 1.25, cex.axis = 1.25, ...){
  if(class(inla.obj) != "inla") stop("an 'inla' object is needed!")
  family <- inla.obj$.args$family
  switch(tolower(family),
         poisson = {
           if (length(inla.obj$.arg$data$baseline.hazard.values) == nrow(inla.obj$summary.random$baseline.hazard)) {
             eval.point <- inla.obj$.arg$data$baseline.hazard.values
             basehaz <- inla.obj$summary.random$baseline.hazard$mean + exp(inla.obj$summary.fixed[1,1])
             sfun <- stepfun(eval.point[-1], basehaz, f=0)
             if (plot) {plot(sfun, xval = eval.point, do.points = F, xaxs = "i", xlim = c(0, max(eval.point)), main="", xlab = "Time", ylab = "Baseline hazard function")}				
           } else {
             eval.point <- inla.obj$.arg$data$baseline.hazard.values
             basehaz <- matrix(inla.obj$summary.random$baseline.hazard$mean + exp(inla.obj$summary.fixed[1,1]), nrow = length(eval.point))
             maxbh <- max(basehaz)
             minbh <- min(basehaz)
             if (plot){
               if (ncol(basehaz) %in% 1:2) {a <- 1; b <- 2
               } else {if (ncol(basehaz) %in% 3:4) {a <- 2; b <- 2}
                 else {if (ncol(basehaz) %in% 5:6) {a <- 2; b <- 3}
                   else {if (ncol(basehaz) %in% 7:9) {a <- 3; b <- 3}
                     else {a <- b <- ceiling(sqrt(ncol(basehaz)))}
                   }}}
               par(mfrow=c(a,b))
               for (i in 1:ncol(basehaz)){
                 sfun <- stepfun(eval.point[-1], basehaz[,i], f=0)
                 plot(sfun, xval = eval.point, do.points = F, xaxs = "i", xlim = c(0, max(eval.point)), main = paste("Stratified group", i), ylim=c(minbh, maxbh), xlab = "Time", ylab = "Baseline hazard function", cex.lab = cex.lab, cex.axis = cex.axis, ...)	
               }									
             }
           }
         },
         weibull = {
           org.formula <- inla.obj$.args$formula
           org.x <- as.character(org.formula)[2]
           splitstr1 <- strsplit(org.x, "(", fixed = TRUE)[[1]]
           splitstr2 <- strsplit(splitstr1[2], ",", fixed = TRUE)[[1]]
           time <- inla.obj$.args$data[,paste(splitstr2[1])]
           
           alpha <- as.numeric(inla.obj$summary.hyperpar[1])
           sigma <- 1/alpha
           mu0 <- -1*inla.obj$summary.fixed[1,1]
           lambda <- exp(-mu0/sigma)
           eval.point <- seq(min(time), max(time), length.out = 101)
           h0 <- lambda*alpha*eval.point^(alpha-1)
           if (plot){
             plot(eval.point, h0, type = "l", xlab = "Time", ylab = "Baseline hazard function", cex.lab = cex.lab, cex.axis = cex.axis, ...)
           }
         },
         exponential = {	
           org.formula <- inla.obj$.args$formula
           org.x <- as.character(org.formula)[2]
           splitstr1 <- strsplit(org.x, "(", fixed = TRUE)[[1]]
           splitstr2 <- strsplit(splitstr1[2], ",", fixed = TRUE)[[1]]
           time <- inla.obj$.args$data[,paste(splitstr2[1])]
           
           mu0 <- -1*inla.obj$summary.fixed[1,1]
           lambda <- exp(-mu0)
           eval.point <- seq(min(time), max(time), length.out = 101)
           h0 <- rep(lambda, length.out = 101)
           if (plot){
             plot(eval.point, h0, type = "l", xlab = "Time", ylab = "Baseline hazard function", cex.lab = cex.lab, cex.axis = cex.axis, ...)
           }
         },
         stop("The function is only support Weibull, Exponential, and Cox proportional hazards models!!!")		
  )
  return(invisible(list(time = eval.point, basehaz = basehaz))) 
}


#' Compute different Bayesian residuals for survival models using INLA
#' @param inla.obj an object of class "inla". 
#' @param time the follow up time for right censored data, 
#' @param event the status indicator, 1=observed event, 0=right censored event
#'
#' @details The function supports Weibull model, Expontenial model and Cox proportional hazards model for right censored data. Three types of residuals are obtained: Cox-Snell residual; Martingale residuals; Deviance residuals. 
#' @return A list of different residuals, and time, event, family information.
#' @export 
bri.surv.resid <- function(inla.obj, time, event){
  if(class(inla.obj) != "inla") stop("an 'inla' object is needed!")
  family <- inla.obj$.args$family
  if (!is.numeric(time)) 
    stop("argument 'time' must be numeric")
  time <- as.vector(time)
  if (!is.numeric(event)) 
    stop("argument 'event' must be numeric")
  event <- as.vector(event)
  switch(tolower(family),
         weibullsurv = {
           if (nrow(inla.obj$summary.fitted.values) != length(time))
             stop("the 'time' variable does not match with the inla object!")
           if (nrow(inla.obj$summary.fitted.values) != length(event))
             stop("the 'event' variable does not match with the inla object!")
           alpha <- as.numeric(inla.obj$summary.hyperpar[1])
           mu0 <- inla.obj$summary.fixed[1,1]
           lambda <- exp(mu0*alpha)
           theta <- inla.obj$summary.fixed[-1,1]*alpha
           # Cox-Snell residual
           cs.resid <- as.vector(lambda*exp(inla.obj$model.matrix[,-1] %*% theta)*(time^alpha))
         },
         exponentialsurv = {
           if (nrow(inla.obj$summary.fitted.values) != length(time))
             stop("the 'time' variable does not match with the inla object!")
           if (nrow(inla.obj$summary.fitted.values) != length(event))
             stop("the 'event' variable does not match with the inla object!")
           mu0 <- inla.obj$summary.fixed[1,1]
           lambda <- exp(mu0)
           theta <- inla.obj$summary.fixed[-1,1]
           # Cox-Snell residual
           cs.resid <- as.vector(lambda*exp(inla.obj$model.matrix[,-1] %*% theta)*time)
         },
         poisson = {
           if (sum(inla.obj$.args$data$baseline.hazard.idx==1) != length(time))
             stop("The 'time' variable does not match with the inla object! The function does not support stratified proportional Hazard Model!")
           if (sum(inla.obj$.args$data$baseline.hazard.idx==1) != length(event))
             stop("The 'event' variable does not match with the inla object! The function does not support stratified proportional Hazard Model!")
           h0 <- approxfun(inla.obj$summary.random$baseline.hazard$ID, inla.obj$summary.random$baseline.hazard$mean + exp(inla.obj$summary.fixed[1,1]))
           H0.est <- unlist(sapply(time, function(x) integrate(h0, low=0, upper = x))[1,])
           idx1 <- which(inla.obj$.arg$data$baseline.hazard.idx == 1)
           mm <- inla.obj$model.matrix[idx1,]
           beta <- inla.obj$summary.fixed[-1,1]
           # Cox-Snell residual
           cs.resid <- as.vector(H0.est*exp(mm[,-1] %*% beta))			
         },
         stop("The function is only support Weibull, Exponential, and Cox proportional hazards models!!!")
  )
  martingale.resid <- event - cs.resid
  deviance.resid <- c(sign(martingale.resid)* sqrt(-2* (martingale.resid+ event * log(cs.resid))))
  
  return(list(cs = cs.resid, martingale = martingale.resid, deviance = deviance.resid, time = time, event = event, family = family))		
}




#' Cox-Snell residual plot
#' 
#' @param resid.obj object from bri.surv.resid function
#' @param lwd the line width, a positive number, defaulting to 2.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ... 
#'
#' @return
#' @export
bri.csresid.plot <- function(resid.obj, lwd = 2, xlab="Cox-Snell Residual", ylab="Estimated Cumulative Hazard Rates", main = "", cex.lab = 1.25, cex.axis = 1.25, ...){
  require(survival)
  cs.resid <- resid.obj$cs
  event <- resid.obj$event
  s.cs.res <- survfit(Surv(cs.resid, event) ~ 1, type="fleming-harrington")
  H.est <- cumsum(s.cs.res$n.event/s.cs.res$n.risk)
  xy.max <- max(c(s.cs.res$time, H.est))
  plot(s.cs.res$time, H.est,type='s',col='black', xlab=xlab, ylab=ylab, lwd=lwd, xlim=c(0,xy.max), ylim=c(0,xy.max), main = main, cex.lab = cex.lab, cex.axis = cex.axis, ...)
  abline(0, 1, lty=2)
}

#' Deviance residual plot
#' 
#' @param resid.obj object from bri.surv.resid function
#' @param covariate the covariate of interest for the residual plot.
#' @param smooth if smooth = TRUE, a nonparametric smooth curve will be added.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ... 
#'
#' @return
#' @export
bri.dresid.plot <- function(resid.obj, covariate = NULL, smooth = FALSE, xlab = NULL, ylab = "Deviance residual", main = "", cex.lab = 1.25, cex.axis = 1.25, ...){	
  d.resid <- resid.obj$deviance
  if (is.null(covariate)) {
    covariate <- seq(1: length(d.resid))
    xlab <- "Index"		
  } else {
    if (length(d.resid) != length(covariate))
      stop("the 'covariate' variable does not match with the residual object!")		
    if (is.null(xlab)) xlab <- "Covariate"
  }
  if (smooth){
    dev.dat <- data.frame(covariate = covariate, deviance = d.resid)
    dev.smooth.inla <- inla(deviance ~ -1 + f(covariate, model = 'rw2', constr = FALSE), data = dev.dat)
    bri.band.plot(dev.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, type = 'random', ylim= c(min(d.resid), max(d.resid)), cex.lab = cex.lab, cex.axis = cex.axis, ...)
    points(dev.dat$covariate, dev.dat$deviance)		
  } else{
    plot(covariate, d.resid, xlab = xlab, main = main, ylab = ylab, ylim= c(min(d.resid), max(d.resid)), cex.lab = cex.lab, cex.axis = cex.axis, ...)
  }	
}


#' Martingale Residual plot
#'
#' @param resid.obj object from bri.surv.resid function
#' @param covariate the covariate of interest for the residual plot.
#' @param smooth if smooth = TRUE, a nonparametric smooth curve will be added.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param main an overall title for the plot.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis The magnification to be used for x and y labels relative to the current setting of cex.
#' @param ... 
#'
#' @return
#' @export
bri.mresid.plot <- function(resid.obj, covariate = NULL, smooth = FALSE, xlab = NULL, ylab = "Martingale residual", main = "", cex.lab = 1.25, cex.axis = 1.25, ...){	
  m.resid <- resid.obj$martingale
  if (is.null(covariate)) {
    covariate <- seq(1: length(m.resid))
    xlab <- "Index"		
  } else {
    if (length(m.resid) != length(covariate))
      stop("the 'covariate' variable does not match with the residual object!")		
    if (is.null(xlab)) xlab <- "Covariate"
  }
  if (smooth){
    mart.dat <- data.frame(covariate = covariate, martingale = m.resid)
    mart.smooth.inla <- inla(martingale ~ -1 + f(covariate, model = 'rw2', constr = FALSE), data = mart.dat)
    bri.band.plot(mart.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, cex.lab = cex.lab, cex.axis = cex.axis, main = main, type = 'random', ylim= c(min(m.resid), max(m.resid)), ...)
    points(mart.dat$covariate, mart.dat$martingale)		
  } else{
    plot(covariate, m.resid, xlab = xlab, main = main, ylab = ylab, cex.lab = cex.lab, cex.axis = cex.axis, ylim= c(min(m.resid), max(m.resid)), ...)
  }	
}


#' Density estimation using INLA
#'
#' @param x the data from which the estimate is to be computed.
#' @param m the number of equally spaced points at which the density is to be estimated.
#' @param from,to the left and right-most points of the grid at which the density is to be estimated.
#' @param cut by default, the values of from and to are cut * diff(range(x)) beyond the extremes of the data.
#' @param diagonal An extra constant added to the diagonal of the precision matrix in INLA.
#' @param constr A boolean variable indicating whater to set a sum to 0 constraint on the term.
#' @param ... 
#'
#' @return a list containing the following components: x - the n coordinates of the points where the density is estimated; y - the estimated density values; y.lower - the estimated 2.5 percentile of density; y.lower - the estimated 97.5 percentile of density.
#' @export

bri.density <- function(x, m = 101, from, to, cut = 0.1, diagonal = 1e-03, constr = T, ...){
  if (any(is.na(x))) stop("'x' contains missing values!")
  if (missing(from)) from <- min(x) - cut * diff(range(x))
  if (missing(to)) to <- max(x) + cut * diff(range(x))
  
  Gmatrix <- function(x, sparse = TRUE){
    if (any(is.na(x))) stop("'x' contains missing values!")
    x <- sort(x)
    n <- length(x)
    d <- diff(x)
    d <- c(Inf, Inf, d, Inf, Inf)
    k <- 3:(n + 2)
    g <- 2/((d[k - 1]^2) * (d[k - 2] + d[k - 1])) + 2/(d[k - 1]*d[k]) * (1/d[k - 1] + 1/d[k]) + 2/((d[k]^2) * (d[k] + d[k + 1]))
    k <- 4:(n + 2)
    g1 <- -2/(d[k - 1]^2) * (1/d[k - 2] + 1/d[k])
    k <- 5:(n + 2)
    g2 <- 2/(d[k - 2] * d[k - 1] * (d[k - 2]+d[k - 1]))
    G <- diag(g)
    G[row(G) == col(G) + 1] <- g1
    G[col(G) == row(G) + 1] <- g1
    G[row(G) == col(G) + 2] <- g2
    G[col(G) == row(G) + 2] <- g2
    if (sparse == TRUE) G <- as(G, "dgTMatrix")
    return(G)
  }
  
  bins=seq(from, to, length.out = m)
  x.bins <- hist(x, breaks = bins, plot=FALSE)
  x.bins.root <- sqrt(x.bins$counts+1/4)
  idx <- 1:length(x.bins.root)
  Q <- Gmatrix(idx)
  inla.fit <- inla(x.bins.root ~ f(idx, model = "generic0", Cmatrix = Q, 
                                   diagonal = diagonal, constr = constr), 
                   data = as.data.frame(list(x.bins.root = x.bins.root, idx = idx)), 
                   control.predictor = list(compute = T))
  inla.est <- inla.fit$summary.linear.predictor[,1]
  inla.lower <- inla.fit$summary.linear.predictor[,3]
  inla.upper <- inla.fit$summary.linear.predictor[,5]
  SimpsonInt <- function (x, f, subdivisions = 256){
    ap <- approx(x, f, n = 2 * subdivisions + 1)
    integral <- diff(ap$x)[1] * (ap$y[2 * (1:subdivisions) - 1] 
                                 + 4 * ap$y[2 * (1:subdivisions)] + ap$y[2 * (1:subdivisions) + 1])/3
    return(sum(integral))
  }
  normalized <- SimpsonInt(x.bins$mids, inla.est^2)
  f <- inla.est^2/normalized
  f.lower <- inla.lower^2/normalized
  f.upper <- inla.upper^2/normalized
  return(structure(list(x = x.bins$mids, y = f, y.lower= f.lower, y.upper=f.upper)))	  
}







