# Robust Bayesian Inference with Weak Instruments in Stata
## 1. Introduction
Instrumental variables remain a popular method for causal inference in econometrics.  They deliver identification with two key assumptions: 1) exogeneity of the instrument and 2) relevance for the endogenous explanatory variable.  It is well known that inference fails when the instruments are only slightly relevant, in which they are called weak.  Frequentist solutions involve testing for instrument strength or developing inference that is robust to instrument weakness.  Giannone, Lenza, and Primiceri (2025) propose a Bayesian alternative, quantifying uncertainty around the target parameter and the degree of endogeneity simultaneously.  This article introduces `bayesweakiv`, a Stata postestimation command that implements their method and demonstrates its use in estimating returns to schooling.

The article is organized as follows.  Section 2 describes the instrumental variables framework and the Bayesian method.  Section 3 documents the Stata command and describes how to use it.  Section 4 gives an example, applying the method to estimating the returns to schooling as in Angrist and Krueger (1991).

## 2. Setup and theory



## 3. The `bayesweakiv` command
The command `bayesweakiv` is implemented as a postestimation command for `ivreg2` and `ivregress`.  Calling `bayesweakiv` directly after an IV regression will use the relevant information returned from the IV regression command to draw from the posterior distribution of the parameters.
The command returns the posterior mean and a credible interval of prespecified coverage for two parameters: the coefficient $\beta$ for the endogenous regressor and the endogeneity parameter $\delta$.  Under the hood, the command reads the relevant information from the IV regression, passes it to a posterior sampling command, and then returns the means and credible intervals.  The intermediate step, the posterior sampling command, returns arrays of posterior draws of $\beta$ and $\delta$.  The .ado file for `bayesweakiv` can easily be altered to return other features of the posterior distributions, such as modes or density plots.

### Syntax
The syntax of the main command is:

`bayesweakiv [, reps(integer 20000) discard(integer 1000) level(real 95)]`

### Options

`reps(integer)` specifies the number of draws to take from the posterior distribution.  Its default is 20,000.

`discard(integer)` specifies the number of initial draws to discard.  Larger numbers will decrease the dependence on the initial conditions of the sampling algorithm.  The default is 1,000.  The resulting sample size from the posterior distribution will be `reps` minus `discard`.  For example, the default parameters yield a posterior sample of 19,000 draws.

`level(real)` chooses the nominal coverage of the returned credible interval.  It takes values in percent, so its default of 95 specifies a 95% coverage.


## 4. Returns to schooling
