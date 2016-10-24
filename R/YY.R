#' Adaptive prior
#'
#' @param x
#' @param degree
#' @param nknot
#' @param theta.prec
#' @param type
#'
#' @return spde object
#' @export
bri.adapt.prior <- function(x, degree=3, nknot=5, theta.prec=0.01, type=c("indpt", "spde", "stat")){
	# xx <- sort(unique(x))
	require('splines')
	x.mesh = inla.mesh.1d(x)
	x.ind = inla.mesh.1d.bary(x.mesh,x,"nearest")$index[,1]

	n <- x.mesh$n
	xx <- x.mesh$loc
	x.fem = inla.mesh.1d.fem(x.mesh)
	H = -x.fem$g1;
	H[1,] = 0; H[n,] = 0; ## Free boundaries.
    G2 = t(H)%*%Diagonal(n, 1/diag(x.fem$c0))%*%H

	epsi <- 1e-8
	M0 <- Diagonal(n, epsi)
	M1 <- Diagonal(n, 0)
	M2 <- G2


	if(type=='indpt'){
		##B-spline basis
		prob <- seq(0, 1,, nknot)[2:(nknot-1)]
		xk <- as.vector(quantile(x, prob=prob)) ## Knots
		basis <- bs(xx, knots=xk, degree=degree, intercept=T)

		## mean and precision for theta
		theta.mu <- rep(0, dim(basis)[2])
		theta.Q <- diag(rep(theta.prec, dim(basis)[2]))
	} else if(type=='spde') {
		# SPDE precision
		prob.th <- seq(0, 1,, nknot)
		xk.th <- as.vector(quantile(x, prob=prob.th))
		xk.mesh = inla.mesh.1d(xk.th)
		basis <- as.matrix(inla.mesh.1d.A(xk.mesh, xx, method='linear'))

		## mean and precision for theta
		theta.mu <- rep(0, dim(basis)[2])
		xk.fem = inla.mesh.1d.fem(xk.mesh)
		nk <- xk.mesh$n
		kappa <- 1
		Hk = -xk.fem$g1;
		Bk = Diagonal(nk, diag(xk.fem$c0))
		Bk.inv = Diagonal(nk, 1/diag(xk.fem$c0))
		theta.Q <- theta.prec*(kappa^2*Bk - kappa*(t(Hk) + Hk) +
				t(Hk)%*%Bk.inv%*%Hk)
	} else{
		basis <- 1
		theta.mu <- 0
		theta.Q <- theta.prec
	}

	B0 <- cBind(0, basis)
	B1 <- cBind(0)
	B2 <- cBind(0)


	spde <- inla.spde2.generic(M0, M1, M2, B0, B1, B2,
			theta.mu, theta.Q, transform="identity", BLC = B0)

	spde$x.ind <- x.ind
	return(spde)
}

#' Plot credible bands for a nonlinear function
#'
#' @param result
#' @param name
#' @param alpha
#' @param ind
#' @param type
#' @param xlab
#' @param ylab
#' @param main
#' @param hpd
#' @param plot.data
#'
#' @return
#' @export

## Function to plot credible band for nonlinear function

bri.band.ggplot <- function(result, name = NULL, alpha = 0.05, ind = NULL, type = c('random', 'fitted', 'linear'), xlab = NULL, ylab = NULL, main = NULL,  hpd = FALSE)
{
  result
  
  if(is.null(ind) == TRUE){
    if(type == 'random'){
      post.summary <- result$summary.random[[name]]
      marg <- result$marginals.random[[name]]
    }
    if(type == 'fitted'){
      post.summary <- result$summary.fitted.values
      marg <- result$marginals.fitted.values
    }
    if(type == 'linear'){
      post.summary <- result$summary.linear.predictor
      marg <- result$marginals.linear.predictor
    }	
  }else{
    if(type == 'random'){
      post.summary <- result$summary.random[[name]][ind,]
      marg <- result$marginals.random[[name]][ind]
    }
    if(type == 'fitted'){
      post.summary <- result$summary.fitted.values[ind,]
      marg <- result$marginals.fitted.values[ind]
    }
    if(type == 'linear'){
      post.summary <- result$summary.linear.predictor[ind,]
      marg <- result$marginals.linear.predictor[ind]
    }
  }	
  
  
  if(hpd == TRUE){
    pp <- 1 - alpha
    tmp <- sapply(marg, function(x) inla.hpdmarginal(pp, x))
    fhat.lb <- tmp[1,]
    fhat.ub <- tmp[2,]
  } else{
    p.min <- alpha/2
    p.max <- 1 - alpha/2
    fhat.lb <- sapply(marg, function(x) inla.qmarginal(p.min, x))
    fhat.ub <- sapply(marg, function(x) inla.qmarginal(p.max, x))
    fhat.lb <- as.vector(fhat.lb)
    fhat.ub <- as.vector(fhat.ub)
  }
  
  require(ggplot2)
  
  fhat <- post.summary$mean	
  
  if(is.null(name) == FALSE){
    xx <- result$summary.random[[name]]$ID
  }else{
    xx <- 1:length(fhat)	
  }
  data.plot <- data.frame(x = xx, fhat = fhat, f.lb = fhat.lb, f.ub = fhat.ub)
  
  
  ggplot(data.plot, aes(x = x)) +
    geom_line(aes(y = fhat)) +
    geom_ribbon(aes(ymin = f.lb, ymax = f.ub), alpha = 0.2) +
    theme_bw(base_size = 20) + labs(x = xlab, y = ylab) + 
    ggtitle(main)
}
#' Plot credible band for nonlinear function
#'
#' @param result
#' @param name
#' @param alpha
#' @param ind
#' @param xlab
#' @param ylab
#' @param main
#' @param sub
#' @param xlim
#' @param ylim
#' @param type
#' @param hpd
#'
#' @return
#' @export
bri.band.plot <- 
  function (result, name = NULL, alpha = 0.05, ind = NULL, xlab = NULL, 
            ylab = NULL, main = NULL, sub = NULL, xlim = NULL, ylim = NULL, cex.lab = 1.25,
            cex.axis = 1.25, type = c("random", "fitted", "linear"), hpd = FALSE, ...) 
  {
    result
    if (is.null(ind) == TRUE) {
      if (type == "random") {
        post.summary <- result$summary.random[[name]]
        marg <- result$marginals.random[[name]]
      }
      if (type == "fitted") {
        post.summary <- result$summary.fitted.values
        marg <- result$marginals.fitted.values
      }
      if (type == "linear") {
        post.summary <- result$summary.linear.predictor
        marg <- result$marginals.linear.predictor
      }
    }
    else {
      if (type == "random") {
        post.summary <- result$summary.random[[name]][ind, 
                                                      ]
        marg <- result$marginals.random[[name]][ind]
      }
      if (type == "fitted") {
        post.summary <- result$summary.fitted.values[ind, 
                                                     ]
        marg <- result$marginals.fitted.values[ind]
      }
      if (type == "linear") {
        post.summary <- result$summary.linear.predictor[ind, 
                                                        ]
        marg <- result$marginals.linear.predictor[ind]
      }
    }
    if (hpd == TRUE) {
      pp <- 1 - alpha
      tmp <- sapply(marg, function(x) inla.hpdmarginal(pp, 
                                                       x))
      fhat.lb <- tmp[1, ]
      fhat.ub <- tmp[2, ]
    }
    else {
      p.min <- alpha/2
      p.max <- 1 - alpha/2
      fhat.lb <- sapply(marg, function(x) inla.qmarginal(p.min, 
                                                         x))
      fhat.ub <- sapply(marg, function(x) inla.qmarginal(p.max, 
                                                         x))
      fhat.lb <- as.vector(fhat.lb)
      fhat.ub <- as.vector(fhat.ub)
    }
    fhat <- post.summary$mean
    if (is.null(name) == FALSE) {
      xx <- result$summary.random[[name]]$ID
    }
    else {
      xx <- 1:length(fhat)
    }
    if (is.null(ylim) == TRUE) {
      ylim <- c(min(fhat.lb), max(fhat.ub))
    }
    plot(xx, fhat, ylab = ylab, xlab = xlab, main = main, sub = sub, 
         cex.lab = cex.lab, cex.axis = cex.axis, xlim = xlim, ylim = ylim, 
         type = "n", ...)
    lines(xx, fhat, lwd = 1.2)
    lines(xx, fhat.lb, lwd = 1.2, lty = 2)
    lines(xx, fhat.ub, lwd = 1.2, lty = 2)
  }

#' TPS Prior
#'
#' @param mesh
#' @param constr
#' @param extraconstr.int
#' @param extraconstr
#' @param theta.mean
#' @param theta.prec
#'
#' @return
#' @export
bri.tps.prior <- function(
    mesh, constr = FALSE,
	extraconstr.int = NULL,
	extraconstr = NULL,
	theta.mean = 0,
	theta.prec = 1e-3){
	B.tau <- cbind(0, 1)
    B.kappa <- cbind(log(1e-6), 0)
    output <- inla.spde2.matern(mesh,
    B.tau=B.tau, B.kappa=B.kappa,
    theta.prior.mean = theta.mean,
    theta.prior.prec = theta.prec)
	output
}

#' Excursions
#'
#' @param result.inla
#' @param name
#' @param ind
#' @param method
#' @param u
#' @param type
#' @param alpha
#'
#' @return
#' @export
excursions.brinla <- function(result.inla, name=NULL,
	ind=NULL, method, u, type, alpha=0.05){
	require(excursions)
	res.exc <- excursions.inla(result.inla, name=name,
		ind=ind, method=method, u=u, type=type)
	F.out <- res.exc$F
	if(name == 'Predictor'){
		x <- 1:dim(result.inla$summary.linear.predictor)[1]
	}else{
		x <- result.inla$summary.random[[name]]$ID
	}
	E.out <- x[res.exc$F >= 1-alpha]
	G.out <- res.exc$F >= 1-alpha
	rho.out <- res.exc$rho
	mean.out <- res.exc$mean
	vars.out <- res.exc$vars
	output <- list(E = E.out, F = F.out, G = G.out,
		rho = rho.out, mean = mean.out, vars = vars.out)
	output
}

