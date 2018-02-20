### Packages we'll need ###
library(INLA)
library(MASS)
library(ggplot2)
library(grid)
library(gridExtra)
theme_set(theme_bw())


## Generalized Linear Models
library(INLA); library(brinla)
# Example 1: Logistic regression -- Low Birth Weight data
# Read the data
data(lowbwt, package = "brinla")
# Frequentist GLM method
lowbwt.glm1 <- glm(LOW ~ AGE + LWT + RACE + SMOKE + HT + UI + FTV, data=lowbwt, family=binomial())
# summary(lowbwt.glm1)
round(coef(summary(lowbwt.glm1)), 4)
# INLA Method
lowbwt.inla1 <- inla(LOW ~ AGE + LWT + RACE + SMOKE + HT + UI + FTV, data=lowbwt, family = "binomial", Ntrials = 1, control.compute = list(dic = TRUE, cpo = TRUE))
# summary(lowbwt.inla1)
# Output Posterior Estimates
round(lowbwt.inla1$summary.fixed, 4)

# Fit the reduced model
# Frequentist GLM method
lowbwt.glm2 <- glm(LOW ~ LWT + RACE + SMOKE + HT + UI, data=lowbwt, family=binomial())
summary(lowbwt.glm2)
# INLA Method
lowbwt.inla2 <- inla(LOW ~ LWT + RACE + SMOKE + HT + UI, data=lowbwt, family = "binomial", Ntrials = 1, control.compute = list(dic = TRUE, cpo = TRUE))
# summary(lowbwt.inla2)
round(lowbwt.inla2$summary.fixed, 4)
# Compare the DICs for two models: the smaller the better
c(lowbwt.inla1$dic$dic, lowbwt.inla2$dic$dic)

## Example 2: Analysis of Australia AIDS Data
# Read the data
data(AIDS, package = "brinla")

# Plot the data
par(mfrow=c(1,2))
plot(DEATHS ~ TIME, data=AIDS, ylim=c(0,60))
hist(AIDS$DEATHS, xlab="DEATHS", main="")


# frequentist method 
AIDS.glm1 <- glm(DEATHS ~ TIME, family=poisson(), data=AIDS) 
# summary(AIDS.glm1)
round(coef(summary(AIDS.glm1)), 4)

# INLA method
AIDS.inla1 <- inla(DEATHS ~ TIME, data = AIDS, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE)) 
# Output Posterior Estimates
round(AIDS.inla1$summary.fixed, 4)

# Plot the results
par(mfrow = c(1,1))
plot(DEATHS ~ TIME, data=AIDS, ylim=c(0,60))
lines(AIDS$TIME, AIDS.inla1$summary.fitted.values$mean, lwd=2)
lines(AIDS$TIME, AIDS.inla1$summary.fitted.values$"0.025quant", lwd=1, lty=2)
lines(AIDS$TIME, AIDS.inla1$summary.fitted.values$"0.975quant", lwd=1, lty=2)

# Poisson Regression with Transformed TIME: INLA method
AIDS.inla2 <- inla(DEATHS ~ log(TIME), data=AIDS, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE)) 
# Output Posterior Estimates
round(AIDS.inla2$summary.fixed, 4)
# Compare DICs
c(AIDS.inla1$dic$dic, AIDS.inla2$dic$dic)

# Poisson Regression with Transformed TIME: frequentist method
AIDS.glm <- glm(DEATHS ~ log(TIME), family=poisson(), data=AIDS) 
round(coef(summary(AIDS.glm)), 4)

# Plot the results
par(mfrow = c(1,1))
plot(DEATHS ~ log(TIME), data = AIDS, ylim=c(0,60))
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$mean, lwd=2)
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$"0.025quant", lwd=1, lty=2)
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$"0.975quant", lwd=1, lty=2)

# Example 3: Analysis for Count Data -- horseshoe crab data
# Read the data
data(crab, package = "brinla")
str(crab)

# Compute mean and variance for the response variable
round(c(mean(crab$SATELLITES), var(crab$SATELLITES)), 4)
with(crab, tapply(SATELLITES, COLOR, function(x){round(mean(x), 4)}))
with(crab, tapply(SATELLITES, COLOR, function(x){round(var(x), 4)}))

# Histograms the response variable "SATELLITES"
library(ggplot2); library(gridExtra)
p1 <- ggplot(crab, aes(x=SATELLITES)) + geom_histogram(binwidth=1, color="black")
p2 <- p1 + facet_wrap(~ COLOR, ncol=2)
grid.arrange(p1, p2, ncol = 2)

# Checking the correlation among predictors
round(cor(crab$WEIGHT, crab$WIDTH),4)
## Since the correlation is high (0.89), we use on WIDTH in the model

# Negative Binomial Regression for crab data
# Frequentist GLM method
data(crab, package = "brinla")
library(MASS)
crab.glm <- glm.nb(SATELLITES ~ COLOR + SPINE + WIDTH, data=crab)
# summary(crab.glm)
round(coef(summary(crab.glm)), 4)
# INLA method
crab.inla1 <- inla(SATELLITES ~ COLOR + SPINE + WIDTH, data = crab, family = "nbinomial", control.compute = list(dic = TRUE, cpo = TRUE))
# Output Posterior Estimates
round(crab.inla1$summary.fixed, 4)
# Output Posterior estimates for the hyperparameter n (size) plays the role of an dispersion parameter, which is (1/overdispersion)
round(crab.inla1$summary.hyperpar, 4)
# OutputPosterior estimates for overdispersion parameter
# mean
round(inla.emarginal(fun=function(x) 1/x, marg=crab.inla1$marginals.hyperpar[[1]]), 4)
# 0.025quant and 0.975quant
overdisp_post <- inla.tmarginal(fun=function(x) 1/x, marg=crab.inla1$marginals.hyperpar[[1]])
round(inla.emarginal(fun=function(x) x, marg=overdisp_post), 4)
round(inla.qmarginal(c(0.025, 0.975), overdisp_post), 4)

# Poisson Regression
# Frequentist GLM method	
crab.glm2 <- glm(SATELLITES ~ COLOR + SPINE + WIDTH, family = poisson(link = log), data = crab)
# summary(crab.glm2)
round(coef(summary(crab.glm2)), 4)
# INLA method
crab.inla2 <- inla(SATELLITES ~ COLOR + SPINE + WIDTH, data = crab, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE))
# Output Posterior Estimates
round(crab.inla2$summary.fixed, 4)
# Compare DICs
c(crab.inla1$dic$dic, crab.inla2$dic$dic)

## # Example 4: Analysis for Rate Data -- Car Insurance claims
data(Insurance, package = "MASS")
# Frequentist GLM method	
insur.glm <- glm(Claims ~ District + Group + Age + offset(log(Holders)), data = Insurance, family = poisson)
round(summary(insur.glm)$coefficients, 4)

# INLA method
insur.inla1 <- inla(Claims ~ District + Group + Age, data = Insurance, family = "poisson", E = Holders)
# Output Posterior Estimates
round(insur.inla1$summary.fixed, 4)

## Bayesian Pearson residual plot
par(mfrow=c(1,2))
insur.bresid <- bri.Pois.resid(insur.inla1, plot = TRUE)
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)
qqnorm(insur.bresid$resid); qqline(insur.bresid$resid)

## Gamma regression
library(faraway)
data(wafer)
# Plot the histogram of the data
par(mfrow = c(1,1))
hist(wafer$resist, prob=T, col="grey", xlab= "resist", main = "")
lines(density(wafer$resist), lwd=2)

# fit gamma regression using the frequentist approach
formula <- resist ~ x1 + x2 + x3 + x4
wafer.glm <- glm(formula, family = Gamma(link = log), data = wafer)
# summary(wafer.glm)
round(coef(summary(wafer.glm)), 4)
round(summary(wafer.glm)$dispersion, 4)

# INLA method
wafer.inla1 <- inla(formula, family = "gamma", data = wafer)
# summary(wafer.inla1)
round(wafer.inla1$summary.fixed, 4)
round(wafer.inla1$summary.hyperpar, 4)

# extract the posterior mean for the dispersion parameter
disp_post <- inla.tmarginal(fun=function(x) 1/x, wafer.inla1$marginals.hyperpar[[1]])
round(inla.emarginal(function(x) x, marg = disp_post),4)
# 0.025quant and 0.975quant
round(inla.qmarginal(c(0.025, 0.975), disp_post), 4)

# Comparing the linear regression with log tranformation.
wafer.inla2 <- inla(log(resist) ~ x1 + x2 + x3 + x4, family = "gaussian", data = wafer)
# summary(wafer.inla2)
round(wafer.inla2$summary.fixed, 4)


## Beta regression
# comparing different beta densities
x <- seq(0.000, 1, length.out =1000)
y1 <- dbeta(x, shape1 = 0.1*50, shape2 = 50*(1-0.1))     
y2 <- dbeta(x, shape1 = 0.25*5, shape2 = 5*(1-0.25)) 
y3 <- dbeta(x, shape1 = 0.5*25, shape2 = 25*(1-0.5)) 
y4 <- dbeta(x, shape1 = 0.75*5, shape2 = 5*(1-0.75))   
y5 <- dbeta(x, shape1 = 0.9*50, shape2 = 50*(1-0.9))  

# plot the densities
par(mfrow = c(1,1))
plot(x, y1, type ="n", xlim =c(0,1), ylim=c(0,10), ylab ="density")
lines(x, y1)
legend(0.07, 9, "(0.1,50)", bty = "n")
lines(x, y2)
legend(0.16, 3.2, "(0.25,5)", bty = "n")
lines(x, y3)
legend(0.39, 5, "(0.5,25)", bty = "n")
lines(x, y4)
legend(0.63, 3.2, "(0.75,5)", bty = "n")
lines(x, y5)
legend(0.715, 9, "(0.9,50)", bty = "n")

# load the GasolineYield data
library(betareg)	
data(GasolineYield)

# plot the histogram of the response variable
par(mfrow= c(1,1))
hist(GasolineYield$yield, prob=T, col="grey", xlab= "yield", main = "", xlim =c(0,1))
lines(density(GasolineYield$yield), lwd=2)

# fit beta regression with MLE
gas.glm <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)
# summary(gas.glm)
coef(summary(gas.glm))
# model checking plots
par(mfrow=c(2,3))
plot(gas.glm, which = 1:6)

# fit the beta model using INLA
gas.inla1 <- inla(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield, family = "beta", control.compute = list(dic = TRUE, cpo = TRUE))
# summary(gas.inla1)
round(gas.inla1$summary.fixed, 4)
round(gas.inla1$summary.hyperpar, 4)

# Check the Bayesian residuals for the INLA model
par(mfrow = c(1,1))
gas.inla1.resid <- bri.beta.resid(gas.inla1, plot = TRUE, ylim = c(-2,2))
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)


# Final model with smalllest DIC
gas.inla2 <- inla(yield ~ temp10 + temp, data = GasolineYield, family = "beta", control.compute = list(dic = TRUE, cpo = TRUE))
# summary(gas.inla2)
round(gas.inla2$summary.fixed, 4)
round(gas.inla2$summary.hyperpar, 4)

# Compare DICs
c(gas.inla1$dic$dic, gas.inla2$dic$dic)

# Check the Bayesian residuals for the final model
par(mfrow = c(1,1))
gas.inla2.resid <- bri.beta.resid(gas.inla2, plot = TRUE, ylim = c(-2,2))
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)

# -------------------
# Zero-inflated Models
# Read the articles data
data(articles, package = "brinla")
table(articles$art)
round(prop.table(table(articles$art)),3)

# Plot the histogram for the response variable
ggplot(articles, aes(art)) + geom_histogram()

# fit the zero-inflated Poisson Model using INLA
articles.inla1 <- inla(art ~ fem + mar + kid5 + phd + ment, data = articles, family = "zeroinflatedpoisson1", control.compute = list(dic = TRUE, cpo = TRUE))
# summary(articles.inla1)
round(articles.inla1$summary.fixed, 4)
round(articles.inla1$summary.hyperpar, 4)

# Poisson distribution assumes that its variance equals its mean, thus the estimated sˆe’s are too low and the inference from the PMM could be mistaken. 

# fit the zero-inflated Negative-Binomial Model using INLA
articles.inla2 <- inla(art ~ fem + mar + kid5 + phd + ment, data = articles, family = "zeroinflatednbinomial1", control.compute = list(dic = TRUE, cpo = TRUE))
# summary(articles.inla2)
round(articles.inla2$summary.fixed, 4)
round(articles.inla2$summary.hyperpar, 4)

# Compare DICs 
c(articles.inla1$dic$dic, articles.inla2$dic$dic)

# Compute the overdispersion parameter and its 95% confidence interval
round(articles.inla2$summary.hyperpar[1,],4)
overdisp_post <- inla.tmarginal(fun = function(x) 1/x, marg = articles.inla2$marginals.hyperpar[[1]])
round(inla.emarginal(fun=function(x) x, marg=overdisp_post), 4)
round(inla.qmarginal(c(0.025, 0.975), overdisp_post), 4)

sessionInfo()
