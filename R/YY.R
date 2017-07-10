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
    Bk = Diagonal(nk, Matrix::diag(xk.fem$c0))
    Bk.inv = Diagonal(nk, 1/Matrix::diag(xk.fem$c0))
    theta.Q <- theta.prec*(kappa^2*Bk - kappa*(t(Hk) + Hk) + 
                             Matrix::t(Hk)%*%Bk.inv%*%Hk)
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


## Function to plot credible band for nonlinear function

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
                              B.tau = B.tau, B.kappa = B.kappa, 
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
excursions.brinla <- function(result.inla, name = NULL, 
                              ind = NULL, method, u, type, alpha = 0.05){
  require(excursions)
  if (!requireNamespace("excursions", quietly = TRUE)) {
    stop("Package excursions needed for this function to work. Please install it.",
         call. = FALSE)
  }
  res.exc <- excursions.inla(result.inla, name=name, 
                             ind=ind, method=method, u=u, type=type)
  
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
  output <- list(E = E.out, F = F.out, G = G.out, 
                 rho = rho.out, mean = mean.out, vars = vars.out)	
  output
}


#' Excursions
#'
#' @param res.exc 
#' @param xlab 
#' @param ylab 
#' @param main 
#'
#' @return
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
#' @param results
#'
#' @return
#' @export				  

map.munich = function(results, ...)
{
  x = as.data.frame(cbind(mapping$V2, results))
  my.drawmap(x, munich.bnd, regionvar=1, plotvar=2, ...)
}

#' Draw map
#'
#' @param data 
#' @param map 
#' @param regionvar 
#' @param plotvar 
#' @param limits 
#' @param cols 
#' @param nrcolors 
#' @param swapcolors 
#' @param pcat 
#' @param hcl.par 
#' @param hsv.par 
#' @param legend 
#' @param drawnames 
#' @param cex.names 
#' @param cex.legend 
#' @param mar.min 
#' @param ... 
#'
#' @return
#' @export

my.drawmap = function (data, map, regionvar = 2, plotvar = 3, limits, cols = "grey", 
                       nrcolors = 100, swapcolors = TRUE, pcat = FALSE, hcl.par = list(h = c(130, 
                                                                                             25), c = 100, l = c(90, 70)), hsv.par = list(s = 1, v = 1), 
                       legend = TRUE, drawnames = FALSE, cex.names = 0.7, cex.legend = 1.5, 
                       mar.min = 2, ...) 
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

#' Map testing
#'
#' @return
#' @export

map.munich.testing = function()
{
  a1 = c(0, 91, 135, 136, 137, 297)
  a2 = c(18, 13, 17 ,19, 78, 90, 332, 370) 
  a3 = c(378, 116, 117, 214, 247, 250, 251, 252, 253, 254, 274, 379)
  
  par(mfrow=c(2,2))
  x = rep(0, 380)
  x[a1+1] = 1
  map.munich(x, cols="grey")
  x = rep(0, 380)
  x[a2+1] = 1
  map.munich(x, cols="grey")
  x = rep(0, 380)
  x[a3+1] = 1
  map.munich(x, cols="grey")
  
}


num.knots <- function(x){
  x.uniq <- unique(x)
  num.knots <- max(5, min(floor(length(unique(x))/4), 35))
  return(num.knots)
}

default.knots <- function(x,num.knots)
{
  if (missing(num.knots))
    num.knots <- max(5, min(floor(length(unique(x))/4), 35))
  
  return(quantile(unique(x), seq(0, 1, length = (num.knots + 2))[-c(1, (num.knots + 2))]))
}

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


