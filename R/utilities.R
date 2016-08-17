#' Convert precision to SD
#'
#' @param prec a precision density
#' @param internal logical indicating whether this is an internal representation
#'
#' @return an SD density
#' @export
inla.hyper.sd = function(prec,internal=FALSE){
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
inla.density.summary = function(dens){
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
inla.hyperpar.summary = function(r){
  irp = r$internal.marginals.hyperpar
  hrp = r$marginals.hyperpar
  hypnames = names(irp)
  iip = grep("precision",hypnames)
  for(i in 1:length(irp)){
    if(i %in% iip){
      irp[[i]] = inla.hyper.sd(irp[[i]],internal=TRUE)
    }else{
      irp[[i]] = hrp[[i]]
      hypnames[i] = names(hrp)[i]
    }
  }
  ts = t(sapply(irp,inla.density.summary))
  hypnames = sub("Log precision","SD",hypnames)
  row.names(ts) = hypnames
  ts
}
