#' Convert precision to SD
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param prec a precision density
#' @param internal logical indicating whether this is an internal representation
#'
#' @return an SD density
#' @export
bri.hyper.sd = function(prec,internal=FALSE){
  if(internal){
    inla.tmarginal(function(x) 1/sqrt(exp(x)),prec)
  }else{
    inla.tmarginal(function(x) 1/sqrt(x), prec)
  }
}

#' Compute a summary from a density
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param dens a density
#'
#' @return numerical summary
#' @export
bri.density.summary = function(dens){
  m = inla.emarginal(function(xx) c(xx, xx^2), dens)
  q = inla.qmarginal(c(0.025, 0.5, 0.975), dens)
  s = sqrt(max(0, m[2] - m[1]^2))
  md = inla.mmarginal(dens)
  c(mean = m[1], sd = s, q0.025 = q[1], q0.5 = q[2], q0.975 = q[3],mode=md)
}

#' Convert precisions to SD in INLA hyperparameter summary
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param inla model object
#'
#' @return summary of hyperparameters on SD scale (where appropriate)
#' @export
bri.hyperpar.summary = function(r){
  irp = r$internal.marginals.hyperpar
  hrp = r$marginals.hyperpar
  hypnames = names(irp)
  iip = grep("precision",hypnames)
  for(i in 1:length(irp)){
    if(i %in% iip){
      irp[[i]] = bri.hyper.sd(irp[[i]],internal=TRUE)
    }else{
      irp[[i]] = hrp[[i]]
      hypnames[i] = names(hrp)[i]
    }
  }
  ts = t(sapply(irp,bri.density.summary))
  hypnames = sub("Log precision","SD",hypnames)
  row.names(ts) = hypnames
  ts
}

#' Plot the hyperparameter posteriors
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param r an INLA model object
#' @param together TRUE if densities to be plotted on a single panel
#'
#' @return data frame containing the densities
#' @export
bri.hyperpar.plot = function(r,together=TRUE){
  if (!require("ggplot2")) stop("Function requires ggplot2 package. Please install this first.")
  irp = r$internal.marginals.hyperpar
  hrp = r$marginals.hyperpar
  hypnames = names(irp)
  iip = grep("precision",hypnames)
  for(i in 1:length(irp)){
    if(i %in% iip){
      irp[[i]] = bri.hyper.sd(irp[[i]],internal=TRUE)
    }else{
      irp[[i]] = hrp[[i]]
      hypnames[i] = names(hrp)[i]
    }
  }
  hypnames = sub("Log precision","SD",hypnames)
  hypnames = sub("the Gaussian observations","error",hypnames)
  names(irp) = hypnames
  cf = data.frame(do.call(rbind,irp))
  cf$parameter = rep(hypnames,times=sapply(irp,nrow))
  if(together){
    p=ggplot(cf,aes(x=x,y=y,linetype=parameter))+geom_line()+ylab("density")+xlab("")
    print(p)
  }else{
    p=ggplot(cf,aes(x=x,y=y))+geom_line()+facet_wrap(~parameter,scales="free")+ylab("density")+xlab("")
    print(p)
  }
  invisible(cf)
}

#' Plot the posterior densities of the random effects
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param r inla model object
#'
#' @return a data frame with the densities and group labels
#' @export
bri.random.plot = function(r){
  if (!require("ggplot2")) stop("Function requires ggplot2 package. Please install this first.")
  reff <- r$marginals.random
  irp = reff[[1]]
  cf = data.frame(do.call(rbind,irp))
  cf$group = rep(as.character(1:length(irp)),times=sapply(irp,nrow))
  p=ggplot(cf,aes(x=x,y=y,linetype=group))+geom_line()+ylab("density")+xlab("")
  print(p)
  invisible(cf)
}

#' Plot posterior densities of the fixed effects
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param r an inla model object
#'
#' @return a data frame containing the densities and parameter labels (invisible)
#' @export
bri.fixed.plot = function(r, together=FALSE){
  if (!require("ggplot2")) stop("Function requires ggplot2 package. Please install this first.")
  rmf = r$marginals.fixed
  cf = data.frame(do.call(rbind, rmf))
  cf$parameter = rep(names(rmf),times=sapply(rmf,nrow))
  if(together){
    p=ggplot(cf,aes(x=x,y=y,linetype=parameter))+geom_line()+geom_vline(xintercept=0)+ylab("density")
    print(p)
  }else{
    p = ggplot(cf,aes(x=x,y=y))+geom_line()+
      facet_wrap(~ parameter, scales="free")+geom_vline(xintercept=0)+ylab("density")
    print(p)
  }
  invisible(cf)
}

#' Gaussian Process Regression in 1D
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param x the predictor vector
#' @param y the response vector
#' @param pcprior limites for the penalised complexity prior (optional). If specified should be a vector 
#' of the form c(r,s) where P(range < r = 0.05) and P(SD(y) > s = 0.05)
#' @param nbasis - number of basis functions for the spline (default is 25)
#' @param degree - degree for splines (default is 2) - allowable possibilities are 0, 1 or 2.
#' @param alpha - controls shape of the GP kernel (default is 2) - 0 < alpha <=2 is possible
#' @param xout - grid on which posterior will be calculated (default is x)
#' @param sigma0 - prior mean for the signal SD (default is SD(y))
#' @param rho0 - prior mean for the range
#'
#' @return list consisting of xout, the posterior mean, the lower 95% credibility band, 
#' the upper 95% credibility band and the INLA object containing the fit
#' @export
bri.gpr <- function(x, y, pcprior, nbasis=25, degree=2, alpha=2, xout=x,
                    sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y)))) 
    stop("missing or infinite values in inputs are not allowed")
  mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis),degree = degree)
  nu <-  alpha - 1/2
  kappa0 <- sqrt(8 * nu)/rho0
  tau0 <- 1 / (4 * kappa0^3 * sigma0^2)^0.5
  if(missing(pcprior)){
    spde <- inla.spde2.matern(mesh, alpha=alpha, constr = FALSE,
                              B.tau = cbind(log(tau0), 1, 0),
                              B.kappa = cbind(log(kappa0), 0, 1),
                              theta.prior.prec = 1e-4)
  }else{
    spde <-  inla.spde2.pcmatern(mesh,alpha=alpha,prior.range=c(pcprior[1],0.05),prior.sigma=c(pcprior[2],0.05))
  }
  A <-  inla.spde.make.A(mesh, loc=x)
  Ap <- inla.spde.make.A(mesh, loc=xout)
  index <-  inla.spde.make.index("sinc", n.spde = spde$n.spde)
  st.est <- inla.stack(data=list(y=y), A=list(A),  effects=list(index),  tag="est")
  st.pred <- inla.stack(data=list(y=NA), A=list(Ap),  effects=list(index),  tag="pred")
  sestpred <- inla.stack(st.est,st.pred)
  formula <-  y ~ -1 + f(sinc, model=spde)
  data <-  inla.stack.data(sestpred)
  result <-  inla(formula, data=data,  family="normal",
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE))
  ii <- inla.stack.index(sestpred, tag='pred')$data
  list(xout=xout,
       mean=result$summary.fitted.values$mean[ii],
       lcb=result$summary.fitted.values$"0.025quant"[ii],
       ucb=result$summary.fitted.values$"0.975quant"[ii],
       inlaobj=result)
}

#' Smoothness bands for Gaussian Process Regression
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param x the predictor vector
#' @param y the response vector
#' @param nbasis - number of basis functions for the spline (default is 25)
#' @param degree - degree for splines (default is 2) - allowable possibilities are 0, 1 or 2.
#' @param alpha - controls shape of the GP kernel (default is 2) - 0 < alpha <=2 is possible
#' @param xout - grid on which posterior will be calculated (default is x)
#' @param sigma0 - prior mean for the signal SD (default is SD(y))
#' @param rho0 - prior mean for the range
#'
#' @return list consisting of xout, the posterior mean, the smoother 95% credibility band, 
#' the rougher 95% credibility band 
#' @export
bri.smoothband <- function(x, y, nbasis=25, degree=2, alpha=2, xout=x,
                           sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y)))) 
    stop("missing or infinite values in inputs are not allowed")
  mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis),degree = degree)
  nu <-  alpha - 1/2
  kappa0 <- sqrt(8 * nu)/rho0
  tau0 <- 1 / (4 * kappa0^3 * sigma0^2)^0.5
  spde <- inla.spde2.matern(mesh, alpha=alpha, constr = FALSE,
                            B.tau = cbind(log(tau0), 1, 0),
                            B.kappa = cbind(log(kappa0), 0, 1),
                            theta.prior.prec = 1e-4)
  A <-  inla.spde.make.A(mesh, loc=x)
  Ap <- inla.spde.make.A(mesh, loc=xout)
  index <-  inla.spde.make.index("sinc", n.spde = spde$n.spde)
  st.est <- inla.stack(data=list(y=y), A=list(A),  effects=list(index),  tag="est")
  st.pred <- inla.stack(data=list(y=NA), A=list(Ap),  effects=list(index),  tag="pred")
  sestpred <- inla.stack(st.est,st.pred)
  formula <-  y ~ -1 + f(sinc, model=spde)
  data <-  inla.stack.data(sestpred)
  result <-  inla(formula, data=data,  family="normal",
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE))
  mres <- inla.spde.result(result,"sinc",spde)
  kappa0 <- exp(mres$summary.log.kappa['0.025quant'])[,]
  sigma02 <- exp(mres$summary.log.variance.nominal['0.5quant'])[,]
  tau0 <- 1 / (4 * kappa0^3 * sigma02)^0.5
  spde <- inla.spde2.matern(mesh, alpha=alpha, constr = FALSE,
                            B.tau = cbind(log(tau0)),
                            B.kappa = cbind(log(kappa0)))
  formula <- y ~ -1 + f(sinc, model=spde)
  resulta <- inla(formula, data=data,  family="normal",
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE))
  kappa0 <- exp(mres$summary.log.kappa['0.975quant'])[,]
  sigma02 <- exp(mres$summary.log.variance.nominal['0.5quant'])[,]
  tau0 <- 1 / (4 * kappa0^3 * sigma02)^0.5
  spde <- inla.spde2.matern(mesh, alpha=alpha, constr = FALSE,
                            B.tau = cbind(log(tau0)),
                            B.kappa = cbind(log(kappa0)))
  formula <- y ~ -1 + f(sinc, model=spde)
  resultb <- inla(formula, data=data,  family="normal",
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE))
  ii <- inla.stack.index(sestpred, tag='pred')$data
  list(xout=xout, 
       mean=result$summary.fitted.values$mean[ii],
       scb=resulta$summary.fitted.values$mean[ii],
       rcb=resultb$summary.fitted.values$mean[ii])
}

#' Non-stationary smoothing for Gaussian Process Regression in 1D
#'
#' @author Julian Faraway, \email{jjf23@bath.ac.uk}
#' @param x the predictor vector
#' @param y the response vector
#' @param nbasis - number of basis functions for the spline (default is 25)
#' @param sbasis - number of basis functions for the smoothing of sigma and rho
#' @param degree - degree for splines (default is 2) - allowable possibilities are 0, 1 or 2.
#' @param alpha - controls shape of the GP kernel (default is 2) - 0 < alpha <=2 is possible
#' @param xout - grid on which posterior will be calculated (default is x)
#'
#' @return list consisting of xout, the posterior mean, the lower 95% credibility band, 
#' the upper 95% credibility band and the INLA object containing the fit
#' @export

bri.nonstat <- function(x, y, nbasis=25, sbasis=5, degree=2, alpha=2, xout=x,
                        sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y)))) 
    stop("missing or infinite values in inputs are not allowed")
  mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis),degree = degree)
  basis.T <-as.matrix(inla.mesh.basis(mesh, type="b.spline",
                                      n=sbasis, degree=2))
  basis.K <-as.matrix(inla.mesh.basis(mesh, type="b.spline",
                                      n=sbasis, degree=2))
  spde <- inla.spde2.matern(mesh, alpha=alpha,
                            B.tau = cbind(basis.T[-1,],0),
                            B.kappa = cbind(0,basis.K[-1,]/2),
                            theta.prior.prec = 1e-4)
  A <-  inla.spde.make.A(mesh, loc=x)
  Ap <- inla.spde.make.A(mesh, loc=xout)
  index <-  inla.spde.make.index("sinc", n.spde = spde$n.spde)
  st.est <- inla.stack(data=list(y=y), A=list(A),  effects=list(index),  tag="est")
  st.pred <- inla.stack(data=list(y=NA), A=list(Ap),  effects=list(index),  tag="pred")
  sestpred <- inla.stack(st.est,st.pred)
  formula <-  y ~ -1 + f(sinc, model=spde)
  data <-  inla.stack.data(sestpred)
  result <-  inla(formula, data=data,  family="normal",
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE))
  ii <- inla.stack.index(sestpred, tag='pred')$data
  list(xout=xout,
       mean=result$summary.fitted.values$mean[ii],
       lcb=result$summary.fitted.values$"0.025quant"[ii],
       ucb=result$summary.fitted.values$"0.975quant"[ii],
       inlaobj=result)
}

