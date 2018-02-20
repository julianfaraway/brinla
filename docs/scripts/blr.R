## Bayesian Linear Regression
library(INLA); library(brinla)
## ------------------------------------------------------------------

# Example 1: Multiple linear regression -- US air pollution data
# Read the data
data(usair, package = "brinla")

# Explortory data analysis
library(ggplot2); library(GGally)
pairs.chart <- ggpairs(usair[,-1], lower = list(continuous = "cor"), upper = list(continuous = "points", combo = "dot")) + ggplot2::theme(axis.text = element_text(size = 6)) 
pairs.chart
# Manuf and Pop are highly correlated: r=0.955; we keep only one into the model.

# multiple linear regression: frequentist approach
usair.formula1 <-  SO2 ~ negtemp + manuf + wind + precip + days
usair.lm1 <- lm(usair.formula1, data = usair)
round(coef(summary(usair.lm1)), 4)
# output the residual standard error
round(summary(usair.lm1)$sigma, 4)

# Bayesian Regression with INLA 
usair.inla1 <- inla(usair.formula1, data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
# summary(usair.inla1)
round(usair.inla1$summary.fixed, 4)
round(usair.inla1$summary.hyperpar, 4)

# compute the posterior estimate of sigma
inla.emarginal(fun = function(x) 1/sqrt(x), marg = usair.inla1$marginals.hyperpar$`Precision for the Gaussian observations`)
round(bri.hyperpar.summary(usair.inla1), 4)
# plot the posterior density of sigma
bri.hyperpar.plot(usair.inla1)

# Change priors
usair.inla2 <- inla(usair.formula1, data = usair, control.compute = list(dic = TRUE, cpo =TRUE), control.fixed = list(mean.intercept = 100, prec.intercept = 10^(-2), mean = list(negtemp = 2, wind = -3, default =0), prec = 1), control.family = list(hyper = list(prec = list(prior="gaussian", param =c(0,1)))))
summary(usair.inla2)

## Prediction in lm 
new.data <- data.frame(negtemp = c(-50, -60, -40), manuf = c(150, 100, 400), pop = c(200, 100, 300), wind = c(6, 7, 8), precip = c(10, 30, 20), days = c(20, 100, 40))
predict(usair.lm1, new.data, se.fit = TRUE)
## Prediction in inla
usair.combined <- rbind(usair, data.frame(SO2 = c(NA, NA, NA), new.data))
usair.link <- c(rep(NA, nrow(usair)), rep(1, nrow(new.data)))
usair.inla1.pred <- inla(usair.formula1, data = usair.combined, control.predictor = list(link = usair.link))
usair.inla1.pred$summary.fitted.values[(nrow(usair)+1):nrow(usair.combined),]

# variable selection with AIC
library(MASS)
usair.step <- stepAIC(usair.lm1, trace = FALSE)
usair.step$anova
# Final multiple regression model
usair.formula2 <- SO2 ~ negtemp + manuf + wind + precip
usair.lm2 <- lm(usair.formula2, data = usair)
round(coef(summary(usair.lm2)), 4)

# Final best model: the smallest DIC (=341.26) in all subsets regression
usair.inla3 <- inla(usair.formula2, data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
round(usair.inla3$summary.fixed, 4)
c(usair.inla1$dic$dic, usair.inla3$dic$dic)

## Posterior Predictive Checking
 
usair.inla3.pred <- inla(usair.formula2, data = usair, control.predictor = list(link = 1, compute = TRUE))
post.predicted.pval <- vector(mode = "numeric", length = nrow(usair))
for(i in (1:nrow(usair))) {
  post.predicted.pval[i] <- inla.pmarginal(q=usair$SO2[i], marginal = usair.inla3.pred$marginals.fitted.values[[i]])
}
hist(post.predicted.pval, main="", breaks = 10, xlab="Posterior predictive p-value")

# Cross-Validation Model Checking
sum(usair.inla3$cpo$failure)

par(mfrow = c(1, 2))
hist(usair.inla3$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(usair.inla3$cpo$pit))), usair.inla3$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(usair.inla3$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))


# Compare LPML for the full model and reduce model
LPML1 <- sum(log(usair.inla1$cpo$cpo)); LPML3 <- sum(log(usair.inla3$cpo$cpo))
c(LPML1, LPML3)

par(mfrow=c(1,1))
plot(usair.inla1$cpo$cpo, usair.inla3$cpo$cpo, xlab="CPO for the full model", ylab="CPO for the reduced model")
abline(0,1)

## Bayesian residual plots
par(mfrow=c(1,2))
bri.lmresid.plot(usair.inla3)
bri.lmresid.plot(usair.inla3, usair$negtemp, xlab = "negtemp", smooth =TRUE)

## --------------------------------------------------
## Robust Linear regression with t-distribution
## Usair data example
usair.inla4 <- inla(usair.formula2, family = "T", data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
# summary(usair.inla4)
round(usair.inla4$summary.fixed, 4)
round(usair.inla4$summary.hyperpar, 4)

## --------------------------------------------------
## Analysis of Variance
## Randomized Complete Block With Two Factors
## PainRelief data

data(painrelief, package = "brinla")
painrelief$PainLevel <- as.factor(painrelief$PainLevel)
painrelief$Codeine <- as.factor(painrelief$Codeine)
painrelief$Acupuncture <- as.factor(painrelief$Acupuncture)

painrelief.lm <- lm(Relief ~ PainLevel + Codeine*Acupuncture, data = painrelief)
summary(painrelief.lm)
anova(painrelief.lm)

# use INLA
painrelief.inla <- inla(Relief ~ PainLevel + Codeine*Acupuncture, data = painrelief)
#summary(painrelief.inla)
round(painrelief.inla$summary.fixed, 4)

est1 <- data.frame(x =c("2", "3", "4", "5", "6", "7", "8"),
                   Estimate = painrelief.inla$summary.fixed[c(2:8),1],
                   L = painrelief.inla$summary.fixed[c(2:8),3],
                   U = painrelief.inla$summary.fixed[c(2:8),5])
p1 <- ggplot(est1, aes(x = x, y = Estimate)) + geom_point(size = 5) + geom_errorbar(aes(ymax = U, ymin = L),size=1) + geom_hline(yintercept=0, size=1.5, col="black") + xlab("Pain Level")
p1


## --------------------------------------------------
## Ridge regression
data(frencheconomy, package = "brinla")
round(cor(frencheconomy[,-1]),4)

# Standardize the predictors
fe.scaled <- cbind(frencheconomy[, 1:2], scale(frencheconomy[, c(-1,-2)]))

# set priors
n <- nrow(frencheconomy)
# for the ridge regression, the priors of the the parameters have a common unknown variance 
# to implement this we have to use copies, and change the data set
fe.scaled$beta1 <- rep(1,n)
fe.scaled$beta2 <- rep(2,n)
fe.scaled$beta3 <- rep(3,n)

# this is the prior for the precision of beta
param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3)))

formula.ridge  = IMPORT ~ f(beta1, DOPROD, model="iid", values = c(1,2,3), hyper = param.beta) + f(beta2, STOCK, copy="beta1", fixed=T) + f(beta3, CONSUM, copy="beta1", fixed=T)
frencheconomy.ridge <- inla(formula.ridge, data=fe.scaled)
ridge.est <- rbind(frencheconomy.ridge$summary.fixed, frencheconomy.ridge$summary.random$beta1[,-1]) 
round(ridge.est,4)

# Comparing with Standard Bayesian Linear regression
formula <- IMPORT ~  DOPROD + STOCK + CONSUM
frencheconomy.inla <- inla(formula, data = fe.scaled, control.fixed = list(prec = 1.0e-3), control.family = list(hyper = param.beta))
round(frencheconomy.inla$summary.fixed, 4)

## Comparing with Ridge regression frequentist approach
library(MASS)
reg2 <- lm.ridge(IMPORT ~  DOPROD + STOCK + CONSUM, data = fe.scaled, lambda = seq(0, 1, length=100))

reg2.final <- lm.ridge(IMPORT ~  DOPROD + STOCK + CONSUM, data = fe.scaled, lambda = reg2$kHKB)
reg2.final

# --------------------------------------------------
# Linear regression with autoregressive errors: New Zealand unemployment data

# Read the data
data(nzunemploy, package = "brinla")
nzunemploy$time <- 1:nrow(nzunemploy)

# Plot the time series data by adult and youth
library(tidyr)
qplot(time, value, data = gather(nzunemploy[,c(2,3,5)], variable, value, -time), geom = "line") + geom_vline(xintercept = 90) + facet_grid(variable ~ ., scale = "free") + ylab("Unemployment rate") + theme_bw()

# Centering predictor
nzunemploy$centeredadult = with(nzunemploy, adult - mean(adult))

# Fit a standard linear regression with independent error
formula1 <- youth ~ centeredadult*policy
nzunemploy.inla1 <- inla(formula1, data= nzunemploy)
# summary(nzunemploy.inla1)
round(nzunemploy.inla1$summary.fixed, 4)
# Plot the Bayesian residuals
par(mfrow=c(1,1))
nzunemploy.res1 <- bri.lmresid.plot(nzunemploy.inla1, type="o")
# Plot the autocorrelation and partial autocorrelation
par(mfrow = c(2,1))
acf(nzunemploy.res1$resid, main = "")
acf(nzunemploy.res1$resid, type = "partial", main="")

# Fit a linear regression with AR(1) error
formula2 <- youth ~ centeredadult*policy + f(time, model = "ar1")
nzunemploy.inla2 <- inla(formula2, data = nzunemploy, control.family = list(hyper = list(prec = list(initial = 15, fixed = TRUE))))
# summary(nzunemploy.inla2)
round(nzunemploy.inla2$summary.fixed, 4)
round(nzunemploy.inla2$summary.hyperpar, 4)

# compare with the frequestist approach
library(nlme)
nzunemploy.gls <- gls(youth ~ centeredadult*policy, correlation = corAR1(form=~1), data = nzunemploy)
summary(nzunemploy.gls)

# plot the fitted lines
ggplot(nzunemploy, aes(centeredadult, youth)) + geom_point(aes(shape = factor(policy)), size = 3) + geom_abline(intercept = nzunemploy.inla2$summary.fixed$mean[1], slope = nzunemploy.inla2$summary.fixed$mean[2]) + geom_abline(intercept = nzunemploy.inla2$summary.fixed$mean[1] + nzunemploy.inla2$summary.fixed$mean[3], slope = nzunemploy.inla2$summary.fixed$mean[2] + nzunemploy.inla2$summary.fixed$mean[4]) + theme_bw()
  
sessionInfo()