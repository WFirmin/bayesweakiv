# Robust Bayesian Inference with Weak Instruments in Stata
## 1. Introduction
Instrumental variables remain a popular method for causal inference in econometrics.  They deliver identification with two key assumptions: 1) exogeneity of the instrument and 2) relevance for the endogenous explanatory variable.  It is well known that inference fails when the instruments are only slightly relevant, in which they are called weak.  Frequentist solutions involve testing for instrument strength or developing inference that is robust to instrument weakness.  Giannone, Lenza, and Primiceri (2025) propose a Bayesian alternative, quantifying uncertainty around the target parameter and the degree of endogeneity simultaneously.  This article introduces `bayesweakiv`, a Stata postestimation command that implements their method and demonstrates its use in estimating returns to schooling.

The article is organized as follows.  Section 2 describes the instrumental variables framework and the Bayesian method.  Section 3 documents the Stata command and describes how to use it.  Section 4 gives an example, applying the method to estimating the returns to schooling as in Angrist and Krueger (1991).

## 2. Setup and theory
This method considers the model

$y=x\beta+c\alpha+\delta\nu+\varepsilon$

$x=z\pi+c\rho+\nu$

with dependent variable $y\in\mathbb{R}^T$, endogenous regressor $x\in R^T$, excluded instruments $z\in R^{T\times k}$, included instruments $c\in R^{T\times l}$, and unobserved errors $\nu,\varepsilon\in R^T$.  The errors are uncorrelated with each other.  The parameters $\beta\in R$, $\delta\in R$, $\alpha\in R^l$, $\rho\in R^l$, and $\pi\in R^k$ are unknown.  The parameter on the endogenous regressor, $\beta$, is the main parameter of interest.  The $\delta$ parameter captures the level of endogeneity of $x$ and is thus of secondary interest.  For $\delta=0$, $x$ is exogenous and we can use ordinary least squares.

Economists commonly use two stage least squares (TSLS) to estimate the system, and the parameter estimates are typically asymptotically normally distributed.  However, when $\pi$ is very close to zero, inference breaks down.  Similarly, a simple Bayesian regression with flat priors leads to undesirable properties for inference. This is the problem of weak instruments.  

Giannone, Lenza, and Primiceri (2025) point out that with a correctly specified likelihood function, the issue must come from the prior.  They propose a nonflat prior for the first stage that induces a flat prior on the concentration parameter $\mu^2=\frac{\pi'z'z\pi}{\sigma_\nu^2}$, which captures the strength of the instruments.  The prior is agnostic about instrument strength, yielding weak instrument robust inference.  

The proposed prior places weight on different values of $\pi$ in different ways depending on the number of instruments.  For $k=2$, no adjustment is needed: a flat prior on $\pi$ induces a flat prior on $\mu^2$.  For $k=1$, mass is shifted away from zero.  When $k>2$, the prior shrinks the first stage coefficients $\pi$ towards zero.

This article presents code to implement their prior following a call to `ivregress` or `ivreg2`.  Pulling the data from the previous command, it knows the number of instruments and applies the appropriate prior.  It draws from the posterior using a Gibbs sampling algorithm provided by Giannone, Lenza, and Primiceri (2025).

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
