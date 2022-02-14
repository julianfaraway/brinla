#' ---
#' title: "Bayesian Regression Modeling with INLA"
#' author: "Chapter Five: Linear Mixed Models"
#' output:
#'   github_document:
#'       toc: true
#' ---
#+ echo=FALSE
knitr::opts_chunk$set(comment=NA, 
                      echo = TRUE,
                      fig.path="figs/", 
                      dev = 'svglite',  
                      fig.ext = ".svg",
                      warning=FALSE, 
                      message=FALSE)
library(svglite)
ggplot2::theme_set(ggplot2::theme_bw())
par(mgp=c(1.5,0.5,0), mar=c(3.1,3.1,3.1,0), pch=20)
#' Code from [Bayesian Regression Modeling with INLA](http://julianfaraway.github.io/brinla/)
#' 
#' Throughout you will see small differences from the output in the published book. Some
#' of this is due to random number generation and some due to minor changes in algorithms.
#' 
#' # Single Random Effect
data(reeds, package="brinla")
summary(reeds)
library(lme4)
mmod <- lmer(nitrogen ~ 1+(1|site), reeds)
summary(mmod)
library(INLA); library(brinla)
formula <- nitrogen ~ 1 + f(site, model="iid")
imod <- inla(formula,family="gaussian", data = reeds)
summary(imod)
#' fixed effects sd and quantiles differ noticeably from the printed book
invsqrt <- function(x) 1/sqrt(x)
sdt <- invsqrt(imod$summary.hyperpar[,-2])
row.names(sdt) <- c("SD of epsilson","SD of site")
sdt
prec.site <- imod$marginals.hyperpar$"Precision for site"
prec.epsilon <- imod$marginals.hyperpar$"Precision for the Gaussian observations"
c(epsilon=inla.emarginal(invsqrt,prec.epsilon),
  site=inla.emarginal(invsqrt,prec.site))
sigma.site <- inla.tmarginal(invsqrt, prec.site)
sigma.epsilon <- inla.tmarginal(invsqrt, prec.epsilon)
c(epsilon=inla.mmarginal(sigma.epsilon),
  site=inla.mmarginal(sigma.site))
#' This function removed from INLA `inla.contrib.sd()`
sampvars <- 1/inla.hyperpar.sample(1000,imod)
sampicc <- sampvars[,2]/(rowSums(sampvars))
quantile(sampicc, c(0.025,0.5,0.975))
bri.hyperpar.summary(imod)
alpha <- data.frame(imod$marginals.fixed[[1]])
library(ggplot2)
#+reedsmarg
ggplot(alpha, aes(x,y)) + geom_line() + geom_vline(xintercept = c(2.66, 3.40)) +
       xlim(2,4)+xlab("nitrogen")+ylab("density")
x <- seq(0,1,len=100)
d1 <- inla.dmarginal(x,sigma.site)
d2 <- inla.dmarginal(x,sigma.epsilon)
rdf <- data.frame(nitrogen=c(x,x),sigma=gl(2,100,
       labels=c("site","epsilon")),density=c(d1,d2))
#+reedssig
ggplot(rdf, aes(x=nitrogen, y=density, linetype=sigma))+geom_line()
bri.hyperpar.plot(imod)
sdres <- sd(reeds$nitrogen)
pcprior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01)))
formula <- nitrogen ~ f(site, model="iid", hyper = pcprior)
pmod <- inla(formula, family="gaussian", data=reeds)
pmod$summary.fixed
bri.hyperpar.summary(pmod)
#+reedshyper
bri.hyperpar.plot(pmod)
(lambda <- 3*sdres/tan(pi*0.99/2))
halfcauchy <- "expression:
              lambda = 0.022;
              precision = exp(log_precision);
              logdens = -1.5*log_precision-log(pi*lambda)-
                        log(1+1/(precision*lambda^2));
              log_jacobian = log_precision;
              return(logdens+log_jacobian);"
hcprior <- list(prec = list(prior = halfcauchy))
formula <- nitrogen ~ f(site, model="iid", hyper = hcprior)
hmod <- inla(formula, family="gaussian", data=reeds)
bri.hyperpar.summary(hmod)
reff <- pmod$marginals.random
x <- seq(-1.5,1.5,len=100)
d1 <- inla.dmarginal(x, reff$site[[1]])
d2 <- inla.dmarginal(x, reff$site[[2]])
d3 <- inla.dmarginal(x, reff$site[[3]])
rdf <- data.frame(nitrogen=x,density=c(d1,d2,d3),site=gl(3,100,
                  labels=LETTERS[1:4]))
#+reedsden
ggplot(rdf, aes(x=nitrogen, y=density, linetype=site))+geom_line()
bri.random.plot(pmod)
sdres <- sd(reeds$nitrogen)
pcprior <- list(prec = list(prior="pc.prec", param = c(3*sdres,0.01)))
formula <- nitrogen ~ f(site, model="iid", hyper = pcprior)
pmod <- inla(formula, family="gaussian", data=reeds, 
        control.compute=list(config = TRUE))
psamp <- inla.posterior.sample(n=1000, pmod)
psamp[[1]]
lvsamp <- t(sapply(psamp, function(x) x$latent))
colnames(lvsamp) <- row.names(psamp[[1]]$latent)
mean(lvsamp[,'site:3'] > lvsamp[,'site:1'])

#' # Longitudinal Data

data(reading, package="brinla")
#+piaddat
ggplot(reading, aes(agegrp, piat, group=id)) + geom_line()
library(lme4)
lmod <- lmer(piat ~ agegrp + (1|id), reading)
summary(lmod)
formula <- piat ~ agegrp + f(id, model="iid")
imod <- inla(formula, family="gaussian", data=reading)
imod$summary.fixed
bri.hyperpar.summary(imod)
#+piathyper
bri.hyperpar.plot(imod)
summary(imod$summary.random$id$mean)
#+piatran
bri.random.plot(imod)
reading$cagegrp <- reading$agegrp - 8.5
lmod <- lmer(piat ~ cagegrp + (cagegrp|id), reading)
summary(lmod)
nid <- length(levels(reading$id))
reading$numid <- as.numeric(reading$id)
reading$slopeid <- reading$numid + nid
formula <- piat ~ cagegrp + f(numid, model="iid2d", n = 2*nid) + 
           f(slopeid, cagegrp, copy="numid")
imod <- inla(formula, family="gaussian", data=reading)
imod$summary.fixed
bri.hyperpar.summary(imod)
postmean <- matrix(imod$summary.random$numid[,2],nid,2)
postmean <- sweep(postmean,2,imod$summary.fixed$mean,"+")
#+piatpost
p <- ggplot(reading, aes(cagegrp, piat, group=id)) + 
     geom_line(color="grey85") + xlab("centered age")
p+geom_abline(intercept=postmean[,1],slope=postmean[,2])
library(gridExtra)
sd.epsilon <- bri.hyper.sd(imod$internal.marginals.hyperpar[[1]],
              internal=TRUE)
sd.intercept <- bri.hyper.sd(imod$internal.marginals.hyperpar[[2]],
                internal=TRUE)
sd.slope <- bri.hyper.sd(imod$internal.marginals.hyperpar[[3]],
            internal=TRUE)
p1 <- ggplot(data.frame(sd.epsilon),aes(x,y))+geom_line()+
      ggtitle("Epsilon")+xlab("piat")+ylab("density")
p2 <- ggplot(data.frame(sd.intercept),aes(x,y))+geom_line()+
      ggtitle("Intercept")+xlab("piat")+ylab("density")
p3 <- ggplot(data.frame(sd.slope),aes(x,y))+geom_line()+
      ggtitle("Slope")+xlab("piat")+ylab("density")
p4 <- ggplot(data.frame(imod$marginals.hyperpar[[4]]),aes(x,y))+
      geom_line()+ggtitle("Rho")+ylab("density")
#+piatfour
grid.arrange(p1,p2,p3,p4,ncol=2)
#+piattog
bri.hyperpar.plot(imod, together=FALSE)
formula <- log(piat) ~ cagegrp + f(numid, model="iid2d", n = 2*nid) + 
           f(slopeid, cagegrp, copy="numid")
imod <- inla(formula, family="gaussian", data=reading)
bri.density.summary(imod$marginals.hyperpar[[4]])
#'
#' # Prediction
#' 
data(reading, package="brinla")
reading$id <- as.numeric(reading$id)
newsub <- data.frame(id=90, agegrp = c(6.5,8.5,10.5), 
          piat=c(18, 25, NA))
nreading <- rbind(reading, newsub)
formula <- piat ~ agegrp + f(id, model="iid")
#' `control.compute` option added to ensure marginals are computed
imod <- inla(formula, family="gaussian", data=nreading, 
        control.predictor = list(compute=TRUE),
        control.compute=list(return.marginals.predictor=TRUE))
pm90 <- imod$marginals.fitted.values[[270]]
p1 <- ggplot(data.frame(pm90),aes(x,y))+geom_line()+xlim(c(20,60))
newsub=data.frame(id=90, agegrp = c(6.5,8.5,10.5), piat=c(NA, NA, NA))
nreading <- rbind(reading, newsub)
formula <- piat ~ agegrp + f(id, model="iid")
#' `control.compute` option added to ensure marginals are computed
imodq <- inla(formula, family="gaussian", data=nreading, 
              control.predictor = list(compute=TRUE),
              control.compute=list(return.marginals.predictor=TRUE))
qm90 <- imodq$marginals.fitted.values[[270]]
# error here
#+readingden
p1+geom_line(data=data.frame(qm90),aes(x,y),linetype=2)+
   xlab("PIAT")+ylab("density")
nsamp <- 10000
randprec <- inla.hyperpar.sample(nsamp, imod)[,1]
neweps <- rnorm(nsamp, mean=0, sd=1/sqrt(randprec))
newobs <- inla.rmarginal(nsamp, pm90) + neweps
dens <- density(newobs)
#+readingden2
p1 + geom_line(data=data.frame(x=dens$x,y=dens$y),aes(x,y),linetype=2)

#' # Classical Z-matrix Model

library(lme4)
mmod <- lmer(nitrogen ~ 1+(1|site), reeds)
Z <- getME(mmod, "Z")
X <- getME(mmod, "X")
n <- nrow(reeds)
formula <- y ~ -1 + X +  f(id.z, model="z",  Z=Z)
imodZ <- inla(formula, data = list(y=reeds$nitrogen, id.z = 1:n, X=X))
bri.hyperpar.summary(imodZ)
imodZ$summary.random
data(meatspec, package="brinla")
trainmeat <- meatspec[1:172,]
testmeat <- meatspec[173:215,]
wavelengths <- seq(850, 1050, length=100)
modlm <- lm(fat ~ ., trainmeat)
rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(predict(modlm,testmeat), testmeat$fat)
#+meatspecwave
plot(wavelengths,coef(modlm)[-1], type="l",ylab="LM Coefficients")
n <- nrow(meatspec)
X <- matrix(1,nrow = n, ncol= 1)
Z <- as.matrix(meatspec[,-101])
y <- meatspec$fat
y[173:215] <- NA
scaley <- 100
formula <- y ~ -1 + X +  f(idx.Z, model="z", Z=Z)
zmod <- inla(formula, data = list(y=y/scaley, idx.Z = 1:n, X=X), control.predictor = list(compute=TRUE))
predb <- zmod$summary.fitted.values[173:215,1]*scaley
rmse(predb, testmeat$fat)
rcoef <- zmod$summary.random$idx.Z[216:315,2]
#+meatspec
plot(wavelengths, rcoef, type="l", ylab="Ridge Coefficients")
int.fixed <- list(prec = list(initial = log(1.0e-9), fixed=TRUE))
u.prec <- list(prec = list(param = c(1.0e-3, 1.0e-3)))
epsilon.prec <- list(prec = list(param = c(1.0e-3, 1.0e-3)))
idx.X <- c(1, rep(NA,100))
idx.Z <- c(NA, 1:100)
scaley <- 100
formula <- y ~ -1 + f(idx.X,  model="iid", hyper = int.fixed) + f(idx.Z,  model="iid", hyper = u.prec)
amod <- inla(formula, data = list(y=y/scaley, idx.X=idx.X, idx.Z=idx.Z), 
        control.predictor = list(A=cbind(X, Z),compute=TRUE), 
        control.family = list(hyper = epsilon.prec))
predb <- amod$summary.fitted.values[173:215,1]
rmse(predb, testmeat$fat/scaley)*scaley

#' # Generalized Linear Mixed Models


#' # Poisson GLMM

data(nitrofen, package="boot")
head(nitrofen)
library(dplyr)
library(tidyr)
lnitrofen <- select(nitrofen, -total) %>% 
             mutate(id=1:nrow(nitrofen)) %>% 
             gather(brood,live,-conc,-id) %>% 
             arrange(id)
lnitrofen$brood <- factor(lnitrofen$brood,labels=1:3)
head(lnitrofen)
lnitrofen$jconc <- lnitrofen$conc + rep(c(-10,0,10),50)
#+nitrodat
ggplot(lnitrofen, aes(x=jconc,y=live, shape=brood)) + 
       geom_point(position = position_jitter(w = 0, h = 0.5)) + 
       xlab("Concentration")
library(lme4)
glmod <- glmer(live ~ I(conc/300)*brood + (1|id), nAGQ=25, 
         family=poisson, data=lnitrofen)
summary(glmod)
formula <- live ~ I(conc/300)*brood + f(id, model="iid")
imod <- inla(formula, family="poisson", data=lnitrofen)
imod$summary.fixed
bri.hyperpar.summary(imod)
library(reshape2)
mf <- melt(imod$marginals.fixed)
cf <- spread(mf,Var2,value)
names(cf)[2] <- 'parameter'
#+nitropost
ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
       scales="free") + geom_vline(xintercept=0) + ylab("density")
#+nitrofix
bri.fixed.plot(imod)
bri.fixed.plot(imod,together=TRUE)
multeff <- exp(imod$summary.fixed$mean)
names(multeff) <- imod$names.fixed
multeff[-1]
sden <- data.frame(bri.hyper.sd(imod$marginals.hyperpar[[1]]))
#+nitropost2
ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
       xlab("linear predictor")
bri.hyperpar.plot(imod)
mden <- data.frame(inla.tmarginal(function(x) exp(1/sqrt(x)), 
        imod$marginals.hyperpar[[1]]))
#+nitromult
ggplot(mden,aes(x,y)) + geom_line() + ylab("density") + 
       xlab("multiplicative")
formula <- live ~ I(conc/300)*brood + f(id, model="iid")
imod <- inla(formula, family="poisson", data=lnitrofen, 
        control.compute=list(dic=TRUE))
formula <- live ~ log(conc+1)*brood + f(id, model="iid")
imod2 <- inla(formula, family="poisson", data=lnitrofen, 
         control.compute=list(dic=TRUE))
c(imod$dic$dic,  imod2$dic$dic)
mreff <- imod$summary.random$id$mean
qqnorm(mreff)
qqline(mreff)
formula <- live ~ I(conc/300)*brood + f(id, model="iid")
imod <- inla(formula, family="poisson", data=lnitrofen, 
        control.compute=list(cpo=TRUE))
#+nitrocpp
plot(log(imod$cpo$cpo),ylab="log(CPO)")
lnitrofen[which.min(imod$cpo$cpo),]
lnitrofen %>% filter(brood == 2, conc==310)
pit <- imod$cpo$pit
n <- length(pit)
uniquant <- (1:n)/(n+1)
logit <- function(p) log(p/(1-p))
#+nitrologit
plot(logit(uniquant), logit(sort(pit)), xlab="uniform quantiles", 
     ylab="Sorted PIT values")
abline(0,1)
which.max(pit)
sdu <- 0.3
pcprior <- list(prec = list(prior="pc.prec", param = c(3*sdu,0.01)))
formula <- live ~ I(conc/300)*brood + f(id, model="iid", 
           hyper = pcprior)
imod2 <- inla(formula, family="poisson", data=lnitrofen)
bri.hyperpar.summary(imod2)
lnitrofen$obsid <- 1:nrow(nitrofen)
formula <- live ~ I(conc/300)*brood + f(id, model="iid") + 
           f(obsid, model="iid")
imodo <- inla(formula, family="poisson", data=lnitrofen)
bri.hyperpar.summary(imodo)

#' # Binary GLMM

data(ohio, package="brinla")
table(ohio$smoke)/4
xtabs(resp ~ smoke + age, ohio)/c(350,187)
library(lme4)
modagh <- glmer(resp ~ age + smoke + (1|id), nAGQ=25, 
          family=binomial, data=ohio)
summary(modagh)
formula <- resp ~ age + smoke + f(id, model="iid")
imod <- inla(formula, family="binomial", data=ohio)
imod$summary.fixed
ilogit <- function(x) exp(x)/(1 + exp(x))
ilogit(imod$summary.fixed[1,c(3,4,5)])
exp(imod$summary.fixed[-1, c(3,4,5)])
bri.hyperpar.summary(imod)
exp(bri.hyperpar.summary(imod)[3:5])
table(xtabs(resp ~ id, ohio))
library(gridExtra)
#+ohiomarg
p1 <- ggplot(data.frame(imod$marginals.fixed[[1]]),aes(x,y)) +
      geom_line()+xlab("logit")+ylab("density")+ggtitle("Intercept")
p2 <- ggplot(data.frame(imod$marginals.fixed[[2]]),aes(x,y)) +
      geom_line()+xlab("logit")+ylab("density")+ggtitle("age")
p3 <- ggplot(data.frame(imod$marginals.fixed[[3]]),aes(x,y)) +
      geom_line()+xlab("logit")+ylab("density")+ggtitle("smoke")
sden <- data.frame(bri.hyper.sd(imod$marginals.hyperpar[[1]]))
p4 <- ggplot(sden,aes(x,y)) + geom_line() + xlab("logit") + 
      ylab("density")+ggtitle("SD(u)")
grid.arrange(p1,p2,p3,p4,ncol=2)
inla.pmarginal(0,imod$marginals.fixed$smoke)
formula <- resp ~ age + smoke + f(id, model="iid")
imod <- inla(formula, family="binomial", data=ohio, 
        control.compute=list(dic=TRUE,waic=TRUE))
formula <- resp ~ factor(age) + smoke + f(id, model="iid")
imod1 <- inla(formula, family="binomial", data=ohio, 
         control.compute=list(dic=TRUE,waic=TRUE))
formula <- resp ~ factor(age)*smoke + f(id, model="iid")
imod2 <- inla(formula, family="binomial", data=ohio, 
         control.compute=list(dic=TRUE,waic=TRUE))
c(imod$dic$dic, imod1$dic$dic, imod2$dic$dic)
c(imod$waic$waic, imod1$waic$waic, imod2$waic$waic)
exp(imod1$summary.fixed[-1, c(3,4,5)])
ohio$obst <- rep(1:4,537)
ohio$repl <- ohio$id + 1
formula <- resp ~ factor(age) + smoke + f(obst, model="ar1", 
           replicate = repl)
imod <- inla(formula, family="binomial", data=ohio)
exp(imod$summary.fixed[-1, c(3,4,5)])
bri.hyperpar.summary(imod)
#' Section 5.7.1 *Improving the approximation* was removed from the INLA
#' package due to stability concerns:
#' <https://groups.google.com/g/r-inla-discussion-group/c/GocEfBc5q00>
sessionInfo()
