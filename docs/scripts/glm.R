library(brinla)

# GLMs


# %Bayesian Analysis


# Binary Responses

data(lowbwt, package = "brinla")
lowbwt.glm1 <- glm(LOW ~ AGE + LWT + RACE + SMOKE + HT + UI + FTV, data=lowbwt, family=binomial())
round(coef(summary(lowbwt.glm1)), 4)
lowbwt.inla1 <- inla(LOW ~ AGE + LWT + RACE + SMOKE + HT + UI + FTV, data=lowbwt, family = "binomial", Ntrials = 1, control.compute = list(dic = TRUE, cpo = TRUE))
round(lowbwt.inla1$summary.fixed, 4)
lowbwt.inla2 <- inla(LOW ~ LWT + RACE + SMOKE + HT + UI, data=lowbwt, family = "binomial", Ntrials = 1, control.compute = list(dic = TRUE, cpo = TRUE))
round(lowbwt.inla2$summary.fixed, 4)
paste("DIC for the model 1: ", round(lowbwt.inla1$dic$dic, 4), "; DIC for model 2: ", round(lowbwt.inla2$dic$dic, 4), sep = "")

# Count Responses

AIDS.inla1 <- inla(DEATHS ~ TIME, data = AIDS, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE)) 
round(AIDS.inla1$summary.fixed, 4)
plot(DEATHS ~ TIME, data=AIDS, ylim=c(0,60))
lines(AIDS$TIME, AIDS.inla$summary.fitted.values$mean, lwd=2)
lines(AIDS$TIME, AIDS.inla$summary.fitted.values$"0.025quant", lwd=1, lty=2)
lines(AIDS$TIME, AIDS.inla$summary.fitted.values$"0.975quant", lwd=1, lty=2)
AIDS.inla2 <- inla(DEATHS ~ log(TIME), data=AIDS, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE)) 
round(AIDS.inla2$summary.fixed, 4)
plot(DEATHS ~ log(TIME), data = AIDS, ylim=c(0,60))
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$mean, lwd=2)
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$"0.025quant", lwd=1, lty=2)
lines(log(AIDS$TIME), AIDS.inla2$summary.fitted.values$"0.975quant", lwd=1, lty=2)
paste("DIC for the model 1: ", round(AIDS.inla1$dic$dic, 4), "; DIC for model 2: ", round(AIDS.inla2$dic$dic, 4), sep = "")
AIDS.glm <- glm(DEATHS ~ log(TIME), family=poisson(), data=AIDS) 
round(coef(summary(AIDS.glm)), 4)
round(c(mean(crab$SATELLITES), var(crab$SATELLITES)), 4)
with(crab, tapply(SATELLITES, COLOR, function(x){
	sprintf("%1.4f(%1.4f)", mean(x), var(x))
}))
(p1 <- ggplot(crab, aes(x=SATELLITES)) + geom_histogram(binwidth=1, color="black"))
(p2 <- p1 + facet_wrap(~ COLOR, ncol=2))
round(cor(crab$WEIGHT, crab$WIDTH),4)
library(MASS)
crab.glm <- glm.nb(SATELLITES ~ COLOR + SPINE + WIDTH, data=crab)
round(coef(summary(crab.glm)), 4)
crab.inla1 <- inla(SATELLITES ~ COLOR + SPINE + WIDTH, data = crab, family = "nbinomial", control.compute = list(dic = TRUE, cpo = TRUE))
round(crab.inla1$summary.fixed, 4)
round(crab.inla1$summary.fixed, 4)
overdisp_post <- inla.tmarginal(fun=function(x) 1/x, marg=crab.inla1$marginals.hyperpar[[1]])
round(inla.emarginal(fun=function(x) x, marg=overdisp_post), 4)
round(inla.qmarginal(c(0.025, 0.975), overdisp_post), 4)
crab.inla2 <- inla(SATELLITES ~ COLOR + SPINE + WIDTH, data = crab, family = "poisson", control.compute = list(dic = TRUE, cpo = TRUE))
round(crab.inla2$summary.fixed, 4)

# Modeling Rates

insur.inla1 <- inla(Claims ~ District + Group + Age, data = Insurance, family = "poisson", E = Holders)
round(insur.inla1$summary.fixed, 4)
par(mfrow=c(1,2))
insur.bresid <- bri.Pois.resid.plot(insur.inla1)
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)
qqnorm(insur.bresid$resid); qqline(insur.bresid$resid)

# Gamma Regression for Skewed Data

library(faraway)
data(wafer)
hist(wafer$resist, prob=T, col="grey", xlab= "resist", main = "")
lines(density(wafer$resist), lwd=2)
formula <- resist ~ x1 + x2 + x3 + x4
wafer.glm <- glm(formula, family = Gamma(link = log), data = wafer)
round(coef(summary(wafer.glm)), 4)
round(summary(wafer.glm)$dispersion, 4)
wafer.inla1 <- inla(formula, family = "gamma", data = wafer)
round(wafer.inla1$summary.fixed, 4)
round(wafer.inla1$summary.hyperpar, 4)
round(inla.emarginal(function(x) 1/x, marg = wafer.inla1$marginals.hyperpar[[1]]),4)
wafer.inla2 <- inla(log(resist) ~ x1 + x2 + x3 + x4, family = "gaussian", data = wafer)
round(wafer.inla2$summary.fixed, 4)

# Proportional Responses

library(betareg)	
data(GasolineYield)
hist(GasolineYield$yield, prob=T, col="grey", xlab= "yield", main = "", xlim =c(0,1))
lines(density(GasolineYield$yield), lwd=2)
gas.glm <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)
coef(summary(gas.glm))
gas.inla <- inla(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield, family = "beta")
round(gas.inla$summary.fixed, 4)
round(gas.inla$summary.hyperpar, 4)
bri.beta.resid.plot(gas.inla, ylim = c(-2,2))
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)
gas.inla2 <- inla(yield ~ temp10 + temp, data = GasolineYield, family = "beta", control.compute = list(dic = TRUE, cpo = TRUE))
round(gas.inla2$summary.fixed, 4)
round(gas.inla2$summary.hyperpar, 4)
paste("DIC for the model 1: ", round(gas.inla1$dic$dic, 4), "; DIC for model 2: ", round(gas.inla2$dic$dic, 4), sep = "")
bri.beta.resid.plot(gas.inla2, ylim = c(-2,2))
abline(1.96, 0, lty=2); abline(-1.96, 0, lty=2)

# Modeling Zero-inflated Data

data(articles, package = "brinla")
table(articles$art)
round(prop.table(table(articles$art)),3)
articles.inla1 <- inla(art ~ fem + mar + kid5 + phd + ment, data = articles, family = "zeroinflatedpoisson1", control.compute = list(dic = TRUE, cpo = TRUE))
round(articles.inla1$summary.fixed, 4)
articles.inla2 <- inla(art ~ fem + mar + kid5 + phd + ment, data = articles, family = "zeroinflatednbinomial1", control.compute = list(dic = TRUE, cpo = TRUE))
round(articles.inla2$summary.fixed, 4)
c(articles.inla1$dic$dic, articles.inla2$dic$dic)
overdisp_post <- inla.tmarginal(fun = function(x) 1/x, marg = articles.inla2$marginals.hyperpar[[1]])
round(inla.emarginal(fun=function(x) x, marg=overdisp_post), 4)
round(inla.qmarginal(c(0.025, 0.975), overdisp_post), 4)

# %Gamma regression

sessionInfo()
