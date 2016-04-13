#' Bayes Linear Model
#'
#' Fits a Bayes Linear Model using direct sampling
#'
#' @param lmfit a linear model object returned by \code{"lm"}
#' @param B number of samples to generate
#'
#' @return A list with the elements
#' \itemize{
#' \item{post.dist}{data frame containing samples from the posterior of the parameters.}
#' \item{summary.stat}{summary statistics for the posterior samples}
#' }
#'
#' @examples
#' cars.lm <- lm(dist ~ speed, data=cars)
#' set.seed(1)
#' cars.blm <- BayesLM.nprior(cars.lm,10000)
#'
#' @export
BayesLM.nprior <- function(lmfit, B){
  require(MASS)
  QR <- lmfit$qr
  df.res <- lmfit$df.residual
  V <- qr.R(QR)
  coef <- lmfit$coef
  Vb<- chol2inv(V)
  s2 <- (t(lmfit$residuals)%*%lmfit$residuals)
  s2<- s2[1,1]/df.res
  sigma2 <- df.res * s2/rchisq(B,df.res)
  coef.post <- data.frame(t(sapply(sigma2,function(x) mvrnorm(1,coef,Vb*x))))
  post.dist <- data.frame(coef.post, sigma = sqrt(sigma2))
  names(post.dist) <- c(names(lmfit$coef), "sigma")
  summary.post <- t(apply(post.dist, 2, function(x)
  {
    c("mean"=mean(x),
      "se"=sd(x),
      "0.025quant"=quantile(x,prob=0.025),
      "median"=median(x),
      "0.975quant"=quantile(x,prob=0.975))
  }))
  print(summary.post)
  return(list(post.dist=post.dist, summary.stat = summary.post))
}
