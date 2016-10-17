#' compute the baseline hazard function in INLA survival models
#' 
#' the function supports Weibull model, Expontenial model and Cox proportional hazards model for right censored data.
#'
#' @param inla.obj 
#' @param plot 
#' @param ... 
#'
#' @return
#' @export
bri.basehaz.plot <- function(inla.obj, plot = TRUE, ...){
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
						plot(sfun, xval = eval.point, do.points = F, xaxs = "i", xlim = c(0, max(eval.point)), main = paste("Stratified group", i), ylim=c(minbh, maxbh), xlab = "Time", ylab = "Baseline hazard function")	
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
				plot(eval.point, h0, type = "l", xlab = "Time", ylab = "Baseline hazard function")
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
				plot(eval.point, h0, type = "l", xlab = "Time", ylab = "Baseline hazard function")
			}
		},
		stop("The function is only support Weibull, Exponential, and Cox proportional hazards models!!!")		
	)
	return(invisible(list(time = eval.point, basehaz = basehaz))) 
}

#' Cox-Snell residual plot
#'
#' @param resid.obj 
#' @param lwd 
#' @param xlab 
#' @param ylab 
#' @param main 
#' @param ... 
#'
#' @return
#' @export
bri.csresid.plot <- function(resid.obj, lwd = 2, xlab="Cox-Snell Residual", ylab="Estimated Cumulative Hazard Rates", main = "", ...){
	require(survival)
	cs.resid <- resid.obj$cs
	event <- resid.obj$event
	s.cs.res <- survfit(Surv(cs.resid, event) ~ 1, type="fleming-harrington")
	H.est <- cumsum(s.cs.res$n.event/s.cs.res$n.risk)
	xy.max <- max(c(s.cs.res$time, H.est))
	plot(s.cs.res$time, H.est,type='s',col='blue', xlab=xlab, ylab=ylab, lwd=lwd, xlim=c(0,xy.max), ylim=c(0,xy.max), main = main, ...)
	abline(0,1,col='red',lty=2)
}

#' Deviance residual plots
#'
#' @param resid.obj 
#' @param covariate 
#' @param smooth 
#' @param xlab 
#' @param ylab 
#' @param main 
#' @param ... 
#'
#' @return
#' @export
bri.dresid.plot <- function(resid.obj, covariate = NULL, smooth = FALSE, xlab = NULL, ylab = "Deviance residual", main = "", ...){	
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
		bri.band.plot(dev.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, type = 'random', ylim= c(min(d.resid), max(d.resid)))
		points(dev.dat$covariate, dev.dat$deviance)		
	} else{
		plot(covariate, d.resid, xlab = xlab, main = main, ylab = ylab, ylim= c(min(d.resid), max(d.resid)), ...)
	}	
}
#' Martingale Residual plot
#' @param resid.obj 
#' @param covariate 
#' @param smooth 
#' @param xlab 
#' @param ylab 
#' @param main 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
bri.mresid.plot <- function(resid.obj, covariate = NULL, smooth = FALSE, xlab = NULL, ylab = "Martingale residual", main = "",...){	
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
		bri.band.plot(mart.smooth.inla, name = 'covariate', alpha = 0.05, xlab = xlab, ylab = ylab, main = main, type = 'random', ylim= c(min(m.resid), max(m.resid)))
		points(mart.dat$covariate, mart.dat$martingale)		
	} else{
		plot(covariate, m.resid, xlab = xlab, main = main, ylab = ylab, ylim= c(min(m.resid), max(m.resid)), ...)
	}	
}

#' Residuals for inla survival models
#' 
#' the function supports Weibull model, Expontenial model and Cox proportional hazards model for right censored data.
#' Three types of residuals could be obtained: Cox-Snell residual; Martingale residuals; Deviance residuals. 

#' @param inla.obj 
#' @param time 
#' @param event 
#'
#' @return
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
		weibull = {
			if (nrow(inla.obj$summary.fitted.values) != length(time))
				stop("the 'time' variable does not match with the inla object!")
			if (nrow(inla.obj$summary.fitted.values) != length(event))
				stop("the 'event' variable does not match with the inla object!")
			alpha <- as.numeric(inla.obj$summary.hyperpar[1])
			mu0 <- inla.obj$summary.fixed[1,1]
			lambda <- exp(mu0*alpha)
			theta <- inla.obj$summary.fixed[-1,1]*alpha
			# Cox-Snell residual
			cs.resid <- c(lambda*exp(inla.obj$model.matrix[,-1] %*% theta)*(time^alpha))						
		},
		exponential = {
			if (nrow(inla.obj$summary.fitted.values) != length(time))
				stop("the 'time' variable does not match with the inla object!")
			if (nrow(inla.obj$summary.fitted.values) != length(event))
				stop("the 'event' variable does not match with the inla object!")
			mu0 <- inla.obj$summary.fixed[1,1]
			lambda <- exp(mu0)
			theta <- inla.obj$summary.fixed[-1,1]
			# Cox-Snell residual
			cs.resid <- c(lambda*exp(inla.obj$model.matrix[,-1] %*% theta)*time)			
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
			cs.resid <- c(H0.est*exp(mm[,-1] %*% beta))			
		},
		stop("The function is only support Weibull, Exponential, and Cox proportional hazards models!!!")
		)
	martingale.resid <- event - cs.resid
	deviance.resid <- c(sign(martingale.resid)* sqrt(-2* (martingale.resid+ event * log(cs.resid))))
	
	return(list(cs = cs.resid, martingale = martingale.resid, deviance = deviance.resid, time = time, event = event, family = family))		
}

