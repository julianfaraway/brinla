#' Convert precision to SD
#'
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
#' @param r an INLA model object
#' @param together TRUE if densities to be plotted on a single panel
#'
#' @return data frame containing the densities
#' @export
bri.hyperpar.plot = function(r,together=TRUE){
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
#' @param r inla model object
#'
#' @return a data frame with the densities and group labels
#' @export
bri.random.plot = function(r){
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
#' @param r an inla model object
#'
#' @return a data frame containing the densities and parameter labels (invisible)
#' @export
bri.fixed.plot = function(r, together=FALSE){
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
