# Examples using INLA for Bayesian Regression Models

by [Julian Faraway](https://julianfaraway.github.io/)


You will need to [install the brinla R package](https://github.com/julianfaraway/brinla) to run the examples described here.

These examples do not correspond exactly with found in the [book](../index.md). Some come from earlier drafts,
some are not found in the book and some are experimental. 

# Worked Examples

- [Introduction](intro.md) - a simple example concerning Hubble's law.
- [Linear regression](chicago.md) - Multiple regression example using some Chicago Insurance data
- [One-way ANOVA](reeds.md) - just one random effect
- [Ridge regression](ridge.md) - Ridge regression with meat spectroscopy data
- [Gaussian Process regression](gpreg.md) - GP regression on fossil data
- [Confidence bands for smoothness](smoothband.md) - GP regression to determine uncertainty regarding smoothness
- [Non-stationary smoothing](nonstat.md) - GP regression with variable smoothing
- [Generalized Extreme Values](gev.md) - fitting maximum annual river flows
- [Define your own prior](prior.md) - using a half Cauchy prior for the SD of a random effect.
- [Linear regression with bounded parameters](frenchpres.md) - French presidential election example with slope parameters bounded in [0,1]

# Linear Mixed Model Examples

*See <https://github.com/julianfaraway/rexamples> for updated versions of these
examples. The estimation procedure in INLA has changed since the output
in the examples below was produced resulting in some noticeable changes.
I keep these here for comparison purposes*

These come from my [Extending Linear Models with R](https://julianfaraway.github.io/faraway/ELM/) book.
I demonstrate these methods for each of the examples in the
text. You'll need to read the text for more background on datasets and the interpretations
or you can just look at the help pages for the datasets. I've focussed attention on the
process for fitting the model and summaries. There's lots more you can do
so these analyses are far from complete.

- [Single Random Effect](oneway.md) - the `pulp` data
- [One Fixed and One Random Effect](rbd.md) - the `penicillin` data
- [Split Plots](split.md) - the `irrigation` data
- [Nested Effects](nested.md) - the `eggs` data
- [Crossed Effects](crossed.md) - the `abrasion` data
- [Multilevel Models](multilevel.md) - the `jsp` data
- [Longitudinal Models](longitudinal.md) - the `psid` data
- [Repeated Measures](repeated.md) - the `vision` data
- [Multiple Response Models](multiple.md) - the `jsp` data
