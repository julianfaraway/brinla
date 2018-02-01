library(brinla)

# Splines as a Linear Mixed Model

library(INLA)
data(age.income, package = 'SemiPar')
str(age.income)
library(brinla)
X <- spline.mixed(age.income$age, degree = 2, type = 'TPB')$X
Z <- spline.mixed(age.income$age, degree = 2, type = 'TPB')$Z
age.income$ID <- 1:nrow(age.income)
formula <- log.income ~ -1 + X + f(ID, model = 'z', Z = Z, hyper = list(prec = list(param = c(1e-6, 1e-6))))
result.tpb <- inla(formula, data = age.income, control.predictor = list(compute = TRUE))
result.tpb$summary.fitted.values
bri.band.plot(result.tpb, x = age.income$age, alpha = 0.05, type = 'fitted', xlab = 'Age (years)', ylab = 'Log(annual income)', ylim = range(age.income$log.income))
points(age.income$age, age.income$log.income, cex = 0.5)
result.tpb$summary.random$ID[-(1:nrow(age.income)),]
X <- spline.mixed(age.income$age, type = 'OSS')$X
Z <- spline.mixed(age.income$age, type = 'OSS')$Z
age.income$ID <- 1:nrow(age.income)
formula <- log.income ~ -1 + X + f(ID, model = 'z', Z = Z, hyper = list(prec = list(param = c(1e-6, 1e-6))))
result.oss <- inla(formula, data = age.income, control.predictor = list(compute = TRUE))
B.tmp <- spline.mixed(age.income$age, type = 'OSS')$B
n.row <- nrow(B.tmp)
n.col <- ncol(B.tmp)
attributes(B.tmp) <- NULL
Bmat <- as(matrix(B.tmp, n.row, n.col), 'sparseMatrix')
Q <- as(spline.mixed(age.income$age, type = 'OSS')$Q, 'sparseMatrix')
formula <- y ~ -1 + f(x, model = 'generic0', Cmatrix = Q, hyper = list(prec = list(param = c(1e-6, 1e-6))))
data.inla <- list(y = age.income$log.income, x = 1:ncol(Bmat))
result.oss2 <- inla(formula, data = data.inla, control.predictor = list(A = Bmat, compute = TRUE))
result.oss2$summary.fitted.values[1:n.row, ]
library(brinla)
bri.band.plot(result.oss2, ind = 1:n.row, x = age.income$age, alpha = 0.05, type = 'fitted', xlab = 'Age (years)', ylab = 'Log(annual income)', ylim = range(age.income$log.income))
points(age.income$age, age.income$log.income, cex = 0.5)

# Analysis of Variance for Functional Data

data(DTI, package = 'refund')
str(DTI)
DTI.sub <- DTI[DTI$Nscans==4 & DTI$case == 1, ]
dat.v1 <- DTI.sub[DTI.sub$visit==1,]$rcst
plot(dat.v1[1,], type = "n", ylim = c(0.1, 1))
for(i in 1:dim(dat.v1)[1]) lines(dat.v1[i,])
y <- as.vector(t(DTI.sub$rcst))
ns <- dim(DTI.sub$rcst)[2]	  ## number of locations
ng <- length(unique(DTI.sub$visit))  ## number of visits (groups)
n <- length(unique(DTI.sub$ID))  ## number of subjects

D1 <- Matrix(rep(1,ng*n),ng*n,1)
A.mu <- kronecker(D1, Diagonal(n=ns, x=1))

D1 <- Diagonal(n = ns, x = 1)
D2 <- Diagonal(n = (ng-1), x = 1)
D3 <- Matrix(rep(0, ng-1), 1, ng-1)
D4 <- kronecker(rBind(D2, D3), D1)
A.a <- kronecker(Matrix(rep(1, n), n, 1), D4)

A <- cBind(A.mu, A.a)
mu <- 1:ns
alpha <- rep(1:ns, ng-1)
alpha.rep <- rep(1:(ng-1), each = ns)
mu2 <- c(mu, rep(NA, length(alpha)))
alpha2 <- c(rep(NA, length(mu)), alpha)
alpha2.rep <- c(rep(NA, length(mu)), alpha.rep)
data.inla <- list(y=y, mu=mu2, alpha=alpha2, alpha.rep=alpha2.rep)
formula <- y ~ -1 + f(mu, model = 'rw2', constr = FALSE, scale.model = T) + f(alpha, model = 'rw2', constr = FALSE, scale.model = TRUE, replicate = alpha.rep)
A1.lc <- kronecker(Matrix(rep(1,ng-1),ng-1,1), Diagonal(n=ns, x=1))
A2.lc <- Diagonal(n = (ng - 1)*ns, x = 1)
lc <- inla.make.lincombs(mu = A1.lc, alpha = A2.lc)
result <- inla(formula, data = data.inla, family = 'beta', control.predictor = list(A = A, compute = TRUE), control.compute = list(config = TRUE), lincomb = lc)
bri.band.ggplot(result, ind = 1:ns, type = 'lincomb') #1st visit
bri.band.ggplot(result, ind = 1:ns + ns, type = 'lincomb') #2nd visit
bri.band.ggplot(result, ind = 1:ns + 2*ns, type='lincomb') #3rd visit
bri.band.ggplot(result, name = 'mu', type = 'random') #4th visit
bri.band.ggplot(result, name='alpha', ind=1:ns, type='random')
bri.band.ggplot(result, name='alpha', ind=1:ns+ns, type='random')
bri.band.ggplot(result, name='alpha', ind=1:ns+2*ns, type='random')
res.exc1 <- excursions.brinla(result, name = 'alpha', ind = 1:ns, u = 0, type = '!=', alpha = 0.05, method = 'NIQC')
res.exc2 <- excursions.brinla(result, name = 'alpha', ind = 1:ns + ns, u = 0, type = '!=', alpha = 0.05, method = 'NIQC')
res.exc3 <- excursions.brinla(result, name = 'alpha', ind = 1:ns + 2*ns, u = 0, type = '!=', alpha = 0.05, method = 'NIQC')
bri.excursions.ggplot(res.exc1)
bri.excursions.ggplot(res.exc2)
bri.excursions.ggplot(res.exc3)

# Extreme Values

data(calder, package="brinla")
plot(Flow ~ WaterYear, calder)
library(INLA)
calder$year <- calder$WaterYear-1973
imod <- inla(Flow ~ 1+year, data=calder, family="gev",scale=0.1)
imod$summary.fixed
inla.pmarginal(0, imod$marginals.fixed$year)
library(brinla)
plot(bri.hyper.sd(imod$marginals.hyperpar$`precision for GEV observations`), type="l", xlab="SD", ylab="density")
plot(imod$marginals.hyperpar$`shape-parameter for gev observations`,type="l", xlim=c(-0.2,0.5), xlab="xi", ylab="density")
pgev <- function(y,xi,tau,eta,sigma=1){
  exp(-(1+xi*sqrt(tau*sigma)*(y-eta))^(-1/xi))
}
yr <- 1997-1973
maxflow <- 173.17
eta <- sum(c(1,yr)*imod$summary.fixed$mean)
tau <- imod$summary.hyperpar$mean[1]
xi <- imod$summary.hyperpar$mean[2]
sigma <- 0.1
pless <- pgev(maxflow, xi, tau, eta,sigma)
1-pless
1/(1-pless)
nsamp <- 999
imod <- inla(Flow ~ 1+year, data=calder, family="gev",scale=0.1, control.compute = list(config=TRUE))
postsamp <- inla.posterior.sample(nsamp, imod)
pps <- t(sapply(postsamp, function(x) c(x$hyperpar, x$latent[42:43])))
colnames(pps) <- c("precision","shape","beta0","beta1")
plot(shape*0.01 ~ precision, pps, ylab="shape")
imod <- inla(Flow ~ 1+year, data=calder, family="gev",scale=0.1, control.compute = list(config=TRUE),control.inla=list(int.strategy='grid', dz=0.2))
postsamp <- inla.posterior.sample(nsamp, imod)
pps <- t(sapply(postsamp, function(x) c(x$hyperpar, x$latent[42:43])))
colnames(pps) <- c("precision","shape","beta0","beta1")
plot(shape*0.01 ~ precision, pps, ylab="shape")
sigma <- 0.1
maxflow <- 173.17
retp = numeric(nsamp)
for(i in 1:nsamp){
  eta <- sum(c(1,yr)*pps[i,3:4])
  tau <- pps[i,1]
  xi <- 0.01*pps[i,2]
  pless <-  pgev(maxflow, xi, tau, eta,sigma)
  retp[i] <- 1/(1-pless)
}
quantile(retp, c(0.025, 0.5, 0.975))

# Density Estimation using INLA

library(brinla)
## Example 1: Density estimation Case 1
set.seed(123)
n <- 500
x <- rnorm(n)
x.den1 <- bri.density(x, cut = 0.3)
x.den2 <- density(x, bw = "SJ")

curve(dnorm(x, mean = 0, sd = 1), from = -4, to = 4, lwd = 3, lty = 3, col = "red", xlab = "x", ylab = "f(x)", cex.lab = 1.5, cex.axis = 1.5, ylim=c(0, 0.45))
lines(x.den1, col = "black", lty = 1, lwd = 3, ylim = c(0, 0.25))
lines(x.den1$x, x.den1$y.upper, lty = 2, lwd = 3)
lines(x.den1$x, x.den1$y.lower, lty = 2, lwd = 3)
lines(x.den2, col = "blue", lty = 4, lwd = 3)

## Example 2: Density estimation Case 2
set.seed(123)
n <- 1000
x <- c(rnorm(n/2, mean = -1.5, sd = 1), rnorm(n/2, mean = 2.5, sd = 0.75))
x.den1 <- bri.density(x, cut = 0.3)
x.den2 <- density(x, bw = "SJ")

curve(dnorm(x, mean = -1.5, sd = 1)/2 + dnorm(x, mean = 2.5, sd = 0.75)/2, from = -6, to = 6, lwd = 3, lty = 3, col = "red", xlab = "x", ylab = "f(x)", ylim=c(0, 0.3), cex.lab = 1.5, cex.axis = 1.5)
lines(x.den1, col = "black", lty = 1, lwd = 3)
lines(x.den1$x, x.den1$y.upper, lty = 2, lwd = 3)
lines(x.den1$x, x.den1$y.lower, lty = 2, lwd = 3)
lines(x.den2, col = "blue", lty = 4, lwd = 3)

# Quantile Regression Models

sessionInfo()
