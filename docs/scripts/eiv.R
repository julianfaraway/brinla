library(brinla)

# Introduction


# Classical Errors-in-Variables Models

set.seed(5)
n = 100
beta0 = 1
beta1 = 5
prec.y = 1
prec.u = 1
prec.x = 1
## generate the true unobserved predictor x
x <- rnorm(n, sd = 1/sqrt(prec.x))
## generate the observed predictor w with heteroscedastic error
d <- runif(n, min = 0.5, max = 1.5)
w <- x + rnorm(n, sd = 1/sqrt(prec.u * d))
## generate the response from the true model with the unobserved x
y <- beta0 + beta1*x + rnorm(n, sd = 1/sqrt(prec.y))
sim.data <- data.frame(y, w, d)
sim.inla <- inla(y ~ w, data = sim.data, family = "gaussian")
summary(sim.inla)
round(sim.inla$summary.fixed, 4)
## set initial values of hyperparameters
init.prec.u <- prec.u
init.prec.x <- var(w) - 1/prec.u
init.prec.y <- sim.inla$summary.hyperpar$mean

## set prior parameters
prior.beta = c(0, 0.0001)
prior.prec.u = c(10, 10)
prior.prec.x = c(10, 10)
prior.prec.y = c(10, 10)

# define the mec model formula
formula = y ~ f(w, model = "mec", scale = d, values = w, 
	hyper = list(
		beta = list(prior = "gaussian", param = prior.beta, fixed = FALSE),
		prec.u = list(prior = "loggamma", param = prior.prec.u, initial = log(init.prec.u), fixed = FALSE),
		mean.x = list(prior = "gaussian", initial = 0, fixed = TRUE),
		prec.x = list(prior = "loggamma", param = prior.prec.x, initial = log(init.prec.x), fixed = FALSE)))
sim.mec.inla <- inla(formula, data = sim.data, family = "gaussian",
	control.family = list(hyper = list(prec = list(param = prior.prec.y, initial = log(init.prec.y), fixed=FALSE))),
	control.fixed = list(mean.intercept = prior.beta[1], prec.intercept = prior.beta[2]))
summary(sim.mec.inla)
plot(w, y, xlim=c(-4, 4), col= "grey")
curve(1 + 5*x, -4,4, add=TRUE, lwd=2)
curve(sim.mec.inla$summary.fixed$mean + sim.mec.inla$summary.hyperpar$mean[2]*x, -4,4, add=TRUE, lwd=2, lty = 5)
curve(sim.inla$summary.fixed$mean[1] + sim.inla$summary.fixed$mean[2]*x, -4,4, add=TRUE, lwd=2, lty = 4)
data(framingham, package = "brinla")
fram <- subset(framingham, select = c("AGE", "SMOKE", "FIRSTCHD"))
fram$W1 <- log((framingham$SBP21 + framingham$SBP22)/2 - 50)
fram$W2 <- log((framingham$SBP31 + framingham$SBP32)/2 - 50)
fram$W <- (fram$W1 + fram$W2)/2
# Set prior paramters
prior.beta <- c(0, 0.01)
prior.lambda <- c(0, 0.01)
# estimate initial values of precision parameters 
prec.u <- 1/(var(fram$W1 - fram$W2)/2)
prec.x <- 1/(var((fram$W1 + fram$W2)/2) - (1/4)*(var((fram$W1 - fram$W2)/2)))
prior.prec.x <- c(prec.x, 1)
prior.prec.u <- c(prec.u, 1)
Y <- matrix(NA, 4*n, 3)
Y[1:n, 1] <- fram$FIRSTCHD
Y[n+(1:n), 2] <- rep(0, n)
Y[2*n+(1:n), 3] <- fram$W1
Y[3*n+(1:n), 3] <- fram$W2

beta.0 <- c(rep(1, n), rep(NA, n), rep(NA, n), rep(NA, n))
beta.x <- c(1:n, rep(NA, n), rep(NA, n), rep(NA, n))
idx.x <- c(rep(NA, n), 1:n, 1:n, 1:n)
weight.x <- c(rep(1, n), rep(-1, n), rep(1, n), rep(1,n))
beta.smoke <- c(fram$SMOKE, rep(NA, n), rep(NA, n), rep(NA,n))
beta.age <- c(fram$AGE, rep(NA, n), rep(NA, n), rep(NA,n))
lambda.0 <- c(rep(NA, n), rep(1, n), rep(NA, n), rep(NA, n))
lambda.smoke <- c(rep(NA, n), fram$SMOKE, rep(NA, n), rep(NA, n))
lambda.age <- c(rep(NA, n), fram$AGE, rep(NA, n), rep(NA, n))
Ntrials <- c(rep(1, n), rep(NA, n), rep(NA, n), rep(NA, n))

fram.jointdata <- data.frame(Y=Y, 
                         beta.0=beta.0, beta.x=beta.x, beta.smoke=beta.smoke, beta.age=beta.age, 
                         idx.x=idx.x, weight.x=weight.x,
                         lambda.0=lambda.0, lambda.smoke=lambda.smoke, lambda.age=lambda.age,
                         Ntrials=Ntrials)
fram.formula <- Y ~  f(beta.x, copy = "idx.x", hyper = list(beta = list(param = prior.beta, fixed = FALSE))) 
	+ f(idx.x, weight.x, model = "iid", values = 1:n, hyper = list(prec = list(initial = -15, fixed = TRUE))) 
	+ beta.0 - 1 + beta.smoke + beta.age + lambda.0 + lambda.smoke + lambda.age
			 
fram.mec.inla <- inla(fram.formula, Ntrials = Ntrials, data = fram.jointdata, family = c("binomial", "gaussian", "gaussian"),
       control.family = list(
         list(hyper = list()),
         list(hyper = list(
           prec = list(initial = log(prec.x),
                       param = prior.prec.x,
                       fixed = FALSE))),
         list(hyper = list(
           prec = list(initial=log(prec.u),
                       param = prior.prec.u,
                       fixed = FALSE)))),
       control.fixed = list(
           mean = list(beta.0=prior.beta[1], beta.smoke=prior.beta[1], beta.age=prior.beta[1], 
                     lambda.0=prior.lambda[1], lambda.smoke=prior.lambda[1], lambda.age=prior.lambda[1]),
           prec = list(beta.0=prior.beta[2], beta.smoke=prior.beta[2], beta.age=prior.beta[2], 
                     lambda.0=prior.lambda[2], lambda.smoke=prior.lambda[2], lambda.age=prior.lambda[2]))
)
fram.mec.inla <-inla.hyperpar(fram.mec.inla)
round(fram.mec.inla$summary.fixed, 4)
round(fram.mec.inla$summary.hyperpar, 4)
fram.naive.inla <- inla(FIRSTCHD~ SMOKE + AGE + W, family = "binomial", data = fram)
round(fram.naive.inla$summary.fixed, 4)

# Berkson Errors-in-Variables Models

data(bronch, package = "brinla")
bronch1 <- subset(bronch, dust >=1.28)
round(prop.table(table(bronch1$cbr)),4)
round(prop.table(table(bronch1$smoking)),4)
round(c(mean(bronch1$dust), sd(bronch1$dust)), 4)
round(cor(bronch1$dust, bronch1$expo), 4)
round(by(bronch1$dust, bronch1$smoking, mean), 4)
prior.beta <- c(0, 0.01) 
prec.u <- 1/1.3
prior.prec.u <- c(1/1.3, 0.01)
bronch.formula <- cbr ~  smoking + expo + f(dust, model="meb", hyper = list(beta = list(param = prior.beta, fixed = FALSE), prec.u = list(param = prior.prec.u, initial = log(prec.u), fixed = FALSE))) 
bronch.meb.inla <- inla(bronch.formula, data = bronch1, family = "binomial", control.fixed = list(mean.intercept = prior.beta[1], prec.intercept = prior.beta[2], mean = prior.beta[1], prec = prior.beta[2]))
bronch.meb.inla <- inla.hyperpar(bronch.meb.inla)
round(bronch.meb.inla$summary.fixed, 4)
round(bronch.meb.inla$summary.hyperpar, 4)
bronch.naive.inla <- inla(cbr ~ dust + smoking + expo, data = bronch1, family = "binomial",
control.fixed = list(mean.intercept = prior.beta[1], prec.intercept = prior.beta[2], mean = prior.beta[1], prec = prior.beta[2]))
round(bronch.naive.inla$summary.fixed, 4)
sessionInfo()
