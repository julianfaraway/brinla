#' Construct a random walk prior for adaptive smoothing 
#' 
#' @param x the predictor vector
#' @param degree the degree of B-spline basis used for the smoothing function (default is 3)
#' @param nknot the number of knots for the B-spline basis (default is 5)
#' @param theta.prec the fixed precision of the Gaussian prior on the smoothing function (default value is 0.01)
#' @param type the type of the Gaussian prior on the smoothing function: 'indpt' independent Gaussian priors; 'spde' SPDE prior
#'
#' @return an object for the 'model' option used in inla()
#' @export
bri.adapt.prior <- function(x, degree=3, nknot=5, theta.prec=0.01, type=c("indpt", "spde")){
  require('splines')
  x.mesh = inla.mesh.1d(x)
  x.ind = inla.mesh.1d.bary(x.mesh,x,"nearest")$index[,1] 
  
  n <- x.mesh$n
  xx <- x.mesh$loc
  x.fem = inla.mesh.1d.fem(x.mesh)
  H = -x.fem$g1; 
  H[1,] = 0; H[n,] = 0; ## Free boundaries.
  G2 = Matrix::t(H)%*%Diagonal(n, 1/Matrix::diag(x.fem$c0))%*%H
  
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
    theta.Q <- Matrix::diag(rep(theta.prec, dim(basis)[2]))
  } else {
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
    Bk = Diagonal(nk, Matrix::diag(xk.fem$c0))
    Bk.inv = Diagonal(nk, 1/Matrix::diag(xk.fem$c0))
    theta.Q <- theta.prec*(kappa^2*Bk - kappa*(Matrix::t(Hk) + Hk) + Matrix::t(Hk)%*%Bk.inv%*%Hk)
  } 
  
  B0 <- cBind(0, basis)
  B1 <- cBind(0)
  B2 <- cBind(0)
  
  spde <- inla.spde2.generic(M0, M1, M2, B0, B1, B2, theta.mu, theta.Q, transform="identity", BLC = B0)
  
  spde$x.ind <- x.ind
  return(spde)
}


#' Plot credible band for a nonlinear function using ggplot()
#'
#' @param result the result object from INLA call
#' @param name the name of the component for which to do the plot
#' @param alpha specifies the credibility level as 100(1-alpha)\%
#' @param ind the indices for the part of a component that should be plotted
#' @param type the type of nonlinear function to be plotted: 'random' random effects; 'fitted' fitted values; 'linear' linear predictors; 'lincomb' linear combinations
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main an overall title for the plot
#' @param hpd TRUE if highest posterior density (HPD) interval is plotted (FALSE by default) 
#'
#' @return A plot of posterior mean and credible band for a nonlinear function estimated by INLA
#' @export
bri.band.ggplot <- function(result, name = NULL, x = NULL, alpha = 0.05, ind = NULL, type = c('random', 'fitted', 'linear', 'lincomb'), xlab = NULL, ylab = NULL, main = NULL,  hpd = FALSE)
{
  require(ggplot2)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }
  
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
    if(type == 'lincomb'){
      post.summary <- result$summary.lincomb.derived
      marg <- result$marginals.lincomb.derived
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
    if(type == 'lincomb'){
      post.summary <- result$summary.lincomb.derived[ind,]
      marg <- result$marginals.lincomb.derived[ind]
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
  fhat <- post.summary$mean	
  
  if(is.null(name) == FALSE){
    xx <- result$summary.random[[name]]$ID
  }else{
    xx <- 1:length(fhat)	
  }
  if(is.null(x) == FALSE){
    xx <- x
  }
  data.plot <- data.frame(x = xx, fhat = fhat, f.lb = fhat.lb, f.ub = fhat.ub)
  
  ggplot(data.plot, aes(x = x)) +
    geom_line(aes(y = fhat)) +
    geom_ribbon(aes(ymin = f.lb, ymax = f.ub), alpha = 0.2) +
    theme_bw(base_size = 20) + labs(x = xlab, y = ylab) + 
    ggtitle(main)
  # print(p)
  # invisible(data.plot)
}
#' Plot credible band for a nonlinear function using plot()
#'
#' @param result the result object from INLA call
#' @param name the name of the component for which to do the plot
#' @param alpha specifies the credibility level as 100(1-alpha)\%
#' @param ind the indices for the part of a component that should be plotted
#' @param x the predictor vector
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @param xlim sets the limits of x axis
#' @param ylim sets the limits of y axis
#' @param type the type of nonlinear function to be plotted: 'random' random effects; 'fitted' fitted values; 'linear' linear predictors
#' @param hpd TRUE if highest posterior density (HPD) interval is plotted (FALSE by default) 
#' @param gray.band TRUE (default) if the credible band is filled with gray color
#'
#' @return A plot of posterior mean and credible band for a nonlinear function estimated by INLA
#' @export
bri.band.plot <- function(result, name = NULL, alpha = 0.05, ind = NULL, x = NULL, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, xlim = NULL, ylim = NULL, cex.lab = 1.5, cex.axis = 1.5, type = c('random', 'fitted', 'linear'), hpd = FALSE, gray.band = TRUE)
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
  
  fhat <- post.summary$mean	
  
  if(is.null(name) == FALSE){
    xx <- result$summary.random[[name]]$ID
  }else{
    xx <- 1:length(fhat)	
  }
  if(is.null(ylim) == TRUE) {
    ylim <- c(min(fhat.lb), max(fhat.ub))
  }
  if(is.null(x) == FALSE){
    xx <- x
  }
  plot(xx, fhat, ylab = ylab, xlab = xlab, main = main, sub = sub, cex.lab = cex.lab, cex.axis = cex.axis, xlim = xlim, ylim = ylim, type = 'n')	
  if(gray.band == TRUE){
    xp <- c(xx, rev(xx))
    yp <- c(fhat.lb, rev(fhat.ub))
    polygon(xp, yp, col='light gray', border='light gray')
    lines(xx, fhat, lwd = 1.2)
  } else{
    lines(xx, fhat, lwd = 1.2)
    lines(xx, fhat.lb, lwd = 1.2, lty = 2)
    lines(xx, fhat.ub, lwd = 1.2, lty = 2)
  }
  
  data.plot <- data.frame(x = xx, fhat = fhat, f.lb = fhat.lb, f.ub = fhat.ub)
  invisible(data.plot)
}

#' Construct a thin-plate spline prior
#'
#' @param mesh the mesh to build the model on
#' @param constr TRUE if apply an integrate-to-zero constraint. Default FALSE.
#' @param extraconstr.int extra field integral constraints
#' @param extraconstr direct linear combination constraints on the basis weights
#' @param theta.mean prior mean for theta
#' @param theta.prec prior precision for theta
#'
#' @return an object for the 'model' option used in inla()
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
                B.tau = B.tau, B.kappa = B.kappa, 
                theta.prior.mean = theta.mean, 
                theta.prior.prec = theta.prec)
  output    
}

#' Excursion sets and contour credible regions for latent Gaussian models
#' 
#' @description A wrapper function for excursions.inla() in excursions package
#' @param result.inla result object from INLA call
#' @param name The name of the component for which to do the calculation
#' @param ind the indices for the part of a component that should be used in the calculations
#' @param method the method for handeling the latent Gaussian structure: 'EB' empirical Bayes; 'QC' quantile correction; 'NI' numerical integration; 'NIQC' numerical integration with quantile correction; 'iNIQC' improved integration with quantle correction
#' @param u the excursion or contour level
#' @param type the type of region: '>' positive excursions; '<' negative excursions; '!=' contour avoiding function; '=' contour credibility function 
#' @param alpha specifies the joint probability of an excursion set or contour credible region as 1-alpha
#'
#' @return a list that contains the following arguments: 'E' excursion set, contour credible region, or contour avoiding set with 1-alpha probability; 'F' the excursion function corresponding to the set E; 'G' the excursion function F thresholded at 1-alpha; 'rho' marginal excursion probabilities; 'mean' posterior mean; 'vars' marginal variances  
#' @export
excursions.brinla <- function(result.inla, name = NULL, ind = NULL, method, u, type, alpha = 0.05){
  require(excursions)
  if (!requireNamespace("excursions", quietly = TRUE)) {
    stop("Package excursions needed for this function to work. Please install it.", call. = FALSE)
  }
  res.exc <- excursions.inla(result.inla, name=name, ind=ind, method=method, u=u, type=type)
  
  if(is.null(ind) == TRUE){
    F.out <- res.exc$F
    if(name == 'Predictor'){
      x <- 1:dim(result.inla$summary.linear.predictor)[1]
    }else{
      x <- result.inla$summary.random[[name]]$ID
    }
    E.out <- x[F.out >= 1-alpha]
    G.out <- F.out >= 1-alpha
    rho.out <- res.exc$rho
    mean.out <- res.exc$mean
    vars.out <- res.exc$vars		
  }else{
    F.out <- res.exc$F[ind]
    if(name == 'Predictor'){
      x <- 1:dim(result.inla$summary.linear.predictor[ind,])[1]
    }else{
      x <- result.inla$summary.random[[name]]$ID[ind]
    }
    E.out <- x[F.out >= 1-alpha]
    G.out <- F.out >= 1-alpha
    rho.out <- res.exc$rho[ind]
    mean.out <- res.exc$mean[ind]
    vars.out <- res.exc$vars[ind]
  }
  output <- list(E = E.out, F = F.out, G = G.out, rho = rho.out, mean = mean.out, vars = vars.out)	
  output
}


#' Plot excursion sets and contour credible regions 
#'
#' @param res.exc result object from excursions.brinla()
#' @param xlab a title for x axis
#' @param ylab a title for y axis
#' @param main an overall title for the plot 
#'
#' @return a plot that shows the excursion set, the corresponding excursion function and the marginal excursion probabilities
#' @export
bri.excursions.ggplot <- function(res.exc, xlab = NULL, ylab = NULL, main = NULL)
{
  require(ggplot2)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }
  
  res.exc
  
  y <- res.exc$F
  x <- 1:length(y)
  yy <- rep(NA, length(y))
  yy[res.exc$G] <- y[res.exc$G]
  z <- res.exc$rho
  
  data.plot <- data.frame(x = x, y = y, z = z, yy = yy)
  
  ggplot(data.plot, aes(x = x, y = y)) + labs(x = xlab, y = ylab) + ggtitle(main) +
    geom_ribbon(aes(x = x, ymax = yy, ymin = 0), alpha = 0.5) + 
    geom_line(aes(x = x, y = y)) + 
    geom_line(aes(x = x, y = z), lty = 2) + 
    theme_bw(base_size = 20) 	
}

#' Draw Munich map
#'
#' @param results the vector of values to be shown on the map
#'
#' @return A map of Munich with color scale
#' @export				  

map.munich = function(results, ...)
{
  x = as.data.frame(cbind(mapping$V2, results))
  my.drawmap(x, munich.bnd, regionvar=1, plotvar=2, ...)
}

#' Draw map (internal function)

my.drawmap = function (data, map, regionvar = 2, plotvar = 3, limits, cols = "grey", nrcolors = 100, swapcolors = TRUE, pcat = FALSE, hcl.par = list(h = c(130, 25), c = 100, l = c(90, 70)), hsv.par = list(s = 1, v = 1), legend = TRUE, drawnames = FALSE, cex.names = 0.7, cex.legend = 1.5, mar.min = 2, ...) 
{
  if (!inherits(map, "bnd")) 
    stop("Argument 'map' is not an object of class 'bnd'!")
  regions <- names(map)
  S <- length(regions)
  is.in <- attr(map, "is.in")
  height2width <- attr(map, "height2width")
  height2width <- height2width * 1.1
  if (length(is.in)) {
    ind <- match(is.in, regions)
    regions <- c(regions[-ind], is.in)
    map <- c(map[-ind], map[ind])
  }
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (!is.null(mar.min)) {
    if (height2width > 1) {
      side <- 17.5 * (1 - 1/height2width) + mar.min/height2width
      par(mar = c(mar.min, side, mar.min, side))
    }
    else {
      top <- 17.5 * (1 - height2width) + mar.min * height2width
      par(mar = c(top, mar.min, top, mar.min))
    }
  }
  black <- grey(0)
  white <- grey(1)
  xmin <- 1:S
  xmax <- 1:S
  ymin <- 1:S
  ymax <- 1:S
  for (i in 1:S) {
    xmin[i] <- min(map[[i]][, 1], na.rm = TRUE)
    xmax[i] <- max(map[[i]][, 1], na.rm = TRUE)
    ymin[i] <- min(map[[i]][, 2], na.rm = TRUE)
    ymax[i] <- max(map[[i]][, 2], na.rm = TRUE)
  }
  xlimits <- c(min(xmin), max(xmax))
  ylimits <- c(min(ymin) - (max(ymax) - min(ymin)) * 0.1, max(ymax))
  if (missing(data)) {
    plot(xlimits, ylimits, type = "n", axes = FALSE, xlab = "", 
         ylab = "", ...)
    for (k in 1:S) polygon(map[[k]][, 1], map[[k]][, 2], 
                           lwd = 0.3, border = black)
  }
  else {
    if (!is.data.frame(data)) 
      data <- read.table(data, header = TRUE)
    ord <- order(data[, regionvar])
    plotvar <- data[, plotvar]
    plotvar <- plotvar[ord]
    regionvar <- data[, regionvar]
    regionvar <- regionvar[ord]
    if (cols != "hcl" && cols != "hsv" && cols != "grey") {
      nrcolors <- length(cols)
      if (swapcolors == TRUE) 
        cols <- rev(cols)
    }
    else {
      if (cols == "hcl") {
        h <- hcl.par$h
        c <- hcl.par$c
        l <- hcl.par$l
      }
      if (cols == "hsv") {
        s <- hsv.par$s
        v <- hsv.par$v
      }
    }
    maxim <- max(plotvar, na.rm = TRUE)
    minim <- min(plotvar, na.rm = TRUE)
    # maxim <- 1.4
    # minim <- -2.1
    if (cols != "hcl" && cols != "hsv" && cols != "grey") {
      upperlimit <- 1
      lowerlimit <- -1
    }
    else {
      if (missing(limits)) {
        lowerlimit <- minim
        upperlimit <- maxim
      }
      else {
        lowerlimit <- limits[1]
        upperlimit <- limits[2]
        if (lowerlimit > minim) {
          plotvar[plotvar < lowerlimit] <- lowerlimit
          cat(paste("Note: lowerlimit is above minimum value (", 
                    lowerlimit, " > ", minim, ")\n"))
        }
        if (upperlimit < maxim) {
          plotvar[plotvar > upperlimit] <- upperlimit
          cat(paste("Note: upperlimit is below maximum value (", 
                    upperlimit, " < ", maxim, ")\n"))
        }
      }
    }
    if (pcat) {
      nrcolors <- 3
      upperlimit <- 1
      lowerlimit <- -1
      if (cols != "hcl" && cols != "hsv" && cols != "grey") 
        cols <- c(cols[1], cols[round(length(cols)/2 + 
                                        0.5)], cols[length(cols)])
    }
    fill.colors <- cut(c(lowerlimit, plotvar, upperlimit), 
                       nrcolors)
    fill.colors <- fill.colors[c(-1, -length(fill.colors))]
    fill.colors <- as.vector(fill.colors, mode = "numeric")
    if (cols != "hcl" && cols != "hsv" && cols != "grey") {
      fill.colors <- cols[fill.colors]
      legend.colors <- cols
    }
    else {
      if (cols == "hcl") {
        if (swapcolors == TRUE) 
          h <- rev(h)
        fill.colors <- colorspace::diverge_hcl(nrcolors, 
                                               h = h, c = c, l = l)[fill.colors]
        legend.colors <- colorspace::diverge_hcl(nrcolors, 
                                                 h = h, c = c, l = l)
      }
      if (cols == "hsv") {
        fill.colors <- (fill.colors - 1)/(3 * (nrcolors - 
                                                 1))
        if (swapcolors == FALSE) 
          fill.colors <- 1/3 - fill.colors
        fill.colors <- hsv(h = fill.colors, s = s, v = v)
        legend.colors <- hsv(h = (0:(nrcolors - 1))/(3 * 
                                                       (nrcolors - 1)), s = s, v = v)
        if (swapcolors == FALSE) 
          legend.colors <- rev(legend.colors)
      }
      if (cols == "grey") {
        fill.colors <- (fill.colors - 1)/(nrcolors - 
                                            1)
        if (swapcolors == TRUE) 
          fill.colors <- 1 - fill.colors
        fill.colors <- grey(fill.colors)
        legend.colors <- grey((0:(nrcolors - 1))/(nrcolors - 
                                                    1))
        if (swapcolors == TRUE) 
          legend.colors <- rev(legend.colors)
      }
    }
    plot(xlimits, ylimits, type = "n", axes = FALSE, col = white, 
         xlab = "", ylab = "", ...)
    if (sum(!is.na(match(regions, regionvar))) == 0) 
      warning("map probably doesn't match datafile")
    block1 <- c()
    block2 <- c()
    for (k in 1:S) {
      if (is.na(map[[k]][1, 1]) && is.na(map[[k]][1, 2])) 
        block2 <- c(block2, k)
      else block1 <- c(block1, k)
    }
    m <- match(regions, regionvar)
    for (k in block1) {
      if (is.na(m[k])) {
        polygon(map[[k]][, 1], map[[k]][, 2], col = white, 
                border = FALSE)
        polygon(map[[k]][, 1], map[[k]][, 2], density = 35, 
                lwd = 0.3, col = black)
      }
      else polygon(map[[k]][, 1], map[[k]][, 2], col = fill.colors[m[k]], 
                   border = black)
    }
    for (k in block2) {
      if (is.na(m[k])) {
        polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = white, 
                border = FALSE)
        polygon(map[[k]][-1, 1], map[[k]][-1, 2], density = 35, 
                lwd = 0.3, col = black)
      }
      else polygon(map[[k]][-1, 1], map[[k]][-1, 2], col = fill.colors[m[k]], 
                   border = black)
    }
    if (legend == TRUE) {
      ylo <- yro <- ylimits[1] + (0.7 * (ylimits[2] - ylimits[1]))/11
      ylu <- yru <- ylimits[1] + (0.3 * (ylimits[2] - ylimits[1]))/11
      tylu <- tyru <- ylimits[1]
      xlu <- xlo <- xlimits[1] + 0.1 * (xlimits[2] - xlimits[1])
      xru <- xro <- xlimits[1] + 0.4 * (xlimits[2] - xlimits[1])
      step <- (xru - xlu)/nrcolors
      for (i in 0:(nrcolors - 1)) {
        polygon(c(xlo + step * i, xlo + step * (i + 1), 
                  xlu + step * (i + 1), xlu + step * i), c(ylo, 
                                                           yro, yru, ylu), col = legend.colors[i + 1], 
                border = legend.colors[i + 1])
      }
      lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, 
                                          ylu, ylo), col = black)
      par(cex = cex.legend)
      text(xlu + 0.5 * step, tylu, round(lowerlimit, 4), 
           col = black)
      text(xru - 0.5 * step, tyru, round(upperlimit, 4), 
           col = black)
      if (lowerlimit + (upperlimit - lowerlimit)/3 < 0 && 
          0 < upperlimit - (upperlimit - lowerlimit)/3) {
        help <- cut(c(0, lowerlimit, upperlimit), nrcolors)
        help <- as.vector(help, mode = "numeric")
        if (nrcolors%%2 == 0) 
          text(xlu + step * (help[1]), tylu, "0", col = black)
        else text(xlu + step * (help[1] - 0.5), tylu, 
                  "0", col = black)
      }
    }
  }
  par(cex = cex.names)
  if (drawnames == TRUE) {
    xpos <- (xmin + xmax)/2
    ypos <- (ymin + ymax)/2
    text(xpos, ypos, labels = regions, col = black)
  }
  return(invisible())
}



#' Choose the number of knots
#' 
#' @description An internal function for spline.mixed()
#' @return number of knots

num.knots <- function(x){
  x.uniq <- unique(x)
  num.knots <- max(5, min(floor(length(unique(x))/4), 35))
  return(num.knots)
}

#' Creat the vector of default knots
#' 
#' @description An internal function for spline.mixed()
#' @return quantiles

default.knots <- function(x,num.knots)
{
  if (missing(num.knots))
    num.knots <- max(5, min(floor(length(unique(x))/4), 35))
  
  return(quantile(unique(x), seq(0, 1, length = (num.knots + 2))[-c(1, (num.knots + 2))]))
}

#' Make matrices in the linear mixed models for truncated power basis splines and O'Sullivan splines
#' 
#' @param z the predictor vector
#' @param degree the degree of B-spline basis
#' @param Nknots the number of knots used in B-spline basis; the default number will be used if not specified
#' @param type the type of spline models: 'TPB' truncated power basis splines; 'OSS' O'Sullivan splines
#' 
#' @return X: matrix for fixed effects; Z: matrix for random effects; Q: penalty matrix in O'Sullivan splines; B: matrix of B-spline basis functions 
#' @export

spline.mixed <- function(z, degree = 3, Nknots = NULL, type = c("TPB", "OSS")){
  require('splines')
  if(is.null(Nknots) == TRUE){
    Nknots <- num.knots(z)
  }
  if(type == 'TPB'){
    N <- length(z)
    X  = matrix(0, nrow = N, ncol = (1 + degree)) 
    Z	= matrix(0, nrow = N, ncol = Nknots)
    knots <- quantile(unique(z), seq(0, 1, length = (Nknots + 2))[-c(1, (Nknots + 2))])
    for(i  in 1:N)
    {
      for(j  in 1:(degree+1))
      {
        X[i,j] = z[i]^(j-1)
      }
      for(j in 1:Nknots)
      {
        Z[i,j] = sign(z[i] - knots[j])*((z[i] - knots[j])^degree)
      }
    }
    colnames(X) <- c("Intercept", paste("x", 1:degree, sep = ""))
  }
  
  if(type == 'OSS'){
    a <- min(z); b <- max(z)
    intKnots <- quantile(unique(z), seq(0, 1, length = (Nknots + 2))[-c(1, (Nknots + 2))])
    
    formOmega <- function(a, b, intKnots)
    {
      allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
      K <- length(intKnots) ; L <- 3*(K+8)
      xtilde <- (rep(allKnots,each = 3)[-c(1,(L-1),L)] + rep(allKnots, each = 3)[-c(1, 2, L)])/2
      wts <- rep(diff(allKnots), each = 3)*rep(c(1, 4, 1)/6, K + 7)
      Bdd <- spline.des(allKnots, xtilde,derivs = rep(2,length(xtilde)),
                        outer.ok = TRUE)$design  
      Omega <- t(Bdd*wts)%*%Bdd     
      return(Omega)
    }
    Omega <- formOmega(a,b,intKnots)
    eigOmega <- eigen(Omega)
    indsZ <- 1:(Nknots + 2)
    UZ <- eigOmega$vectors[,indsZ]
    LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
    B <- bs(z,knots = intKnots, degree = degree, Boundary.knots = c(a, b), intercept = TRUE)
    ZS <- B%*%LZ   
    Z <- ZS
    X <- cbind(rep(1, length(z)), z)
    colnames(X) <- c("Intercept", "x")
  }
  if(type == 'OSS'){
    return(list(X = X, Z = Z, Q = Omega, B = B))
  }else{
    return(list(X = X, Z = Z))
  }
}

