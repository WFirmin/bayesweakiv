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

### Stored Results
The `bayesweakiv` command stores the following in `e()`:

Scalars

`e(beta_mean)`: the sample mean of the posterior draws for $\beta$.

`e(beta_lo)`: the lower endpoint of the credible interval for $\beta$.

`e(beta_hi)`: the upper endpoint of the credible interval for $\beta$.

`e(delta_mean)`: the sample mean of the posterior draws for $\delta$.

`e(delta_lo)`: the lower endpoint of the credible interval for $\delta$.

`e(delta_hi)`: the upper endpoint of the credible interval for $\delta$.

Matrices

`e(beta_draws)`: the draws from the posterior distribution of $\beta$.

`e(delta_draws)`: the draws from the posterior distribution of $\delta$.



## 4. Returns to schooling

This section replicates selected specifications from:
 
> Angrist, J. D., and A. B. Krueger. 1991. Does compulsory school attendance affect schooling and earnings? *Quarterly Journal of Economics* 106(4): 979–1014.
 
Angrist and Krueger use quarter of birth as an instrument for years of schooling in a wage equation, exploiting the fact that students born earlier in the year are forced to stay in school longer due to compulsory attendance laws. The paper is a landmark application of IV methods and a canonical weak-instruments benchmark.

We download the data and replication files from https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive.  We run the do files with the simple addition of the weakbayesiv command after each 2SLS regression.  We use 200,000 draws of the posterior and discard the first 20,000 for each specification.

Angrist and Krueger perform the same analysis for several cohorts.  We focus on the first cohort, that analyzed in table 4, for simplicity.  The `bayesweakiv` command fits in with minimal changes.  The final lines of the replication code run four 2SLS specifications:

```
** Col 2 4 6 8 ***
ivregress 2sls LWKLYWGE YR20-YR28 (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
ivregress 2sls LWKLYWGE YR20-YR28 AGEQ AGEQSQ (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
ivregress 2sls LWKLYWGE YR20-YR28 RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT  (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
ivregress 2sls LWKLYWGE YR20-YR28 RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT AGEQ AGEQSQ (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
```

With the simple addition of `bayesweakiv` between each call to `ivregress`, we can produce robust credible intervals.  We add the code as follows, also specifying the number of repetitions and discarded draws to be 200,000 and 20,000.

```
** Col 2 4 6 8 ***
ivregress 2sls LWKLYWGE YR20-YR28 (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
bayesweakiv, reps(200000) discard(20000)
ivregress 2sls LWKLYWGE YR20-YR28 AGEQ AGEQSQ (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
bayesweakiv, reps(200000) discard(20000)
ivregress 2sls LWKLYWGE YR20-YR28 RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT  (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
bayesweakiv, reps(200000) discard(20000)
ivregress 2sls LWKLYWGE YR20-YR28 RACE MARRIED SMSA NEWENG MIDATL ENOCENT WNOCENT SOATL ESOCENT WSOCENT MT AGEQ AGEQSQ (EDUC = QTR120-QTR129 QTR220-QTR229 QTR320-QTR329 YR20-YR28)
bayesweakiv, reps(200000) discard(20000)
```

The output for the second 2SLS regression (column 4) and its posterior is
```
. ivregress 2sls LWKLYWGE YR20-YR28 AGEQ AGEQSQ (EDUC = QTR120-QTR129 QTR220-QTR
> 229 QTR320-QTR329 YR20-YR28)
note: QTR327 omitted because of collinearity.
note: QTR329 omitted because of collinearity.

Instrumental-variables 2SLS regression            Number of obs   =    247,199
                                                  Wald chi2(12)   =     104.08
                                                  Prob > chi2     =     0.0000
                                                  R-squared       =     0.1023
                                                  Root MSE        =     .61707

------------------------------------------------------------------------------
    LWKLYWGE | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
        EDUC |   .1310424    .033356     3.93   0.000     .0656659     .196419
        YR20 |  -.1134065   .0709263    -1.60   0.110    -.2524194    .0256064
        YR21 |  -.1081955   .0649955    -1.66   0.096    -.2355843    .0191934
        YR22 |   -.103873    .055718    -1.86   0.062    -.2130782    .0053322
        YR23 |  -.0938903   .0495893    -1.89   0.058    -.1910835    .0033029
        YR24 |   -.080653   .0425868    -1.89   0.058    -.1641217    .0028156
        YR25 |  -.0573875   .0337871    -1.70   0.089     -.123609    .0088341
        YR26 |  -.0427072   .0271859    -1.57   0.116    -.0959906    .0105761
        YR27 |  -.0188154   .0175163    -1.07   0.283    -.0531468     .015516
        YR28 |   .0003555   .0102907     0.03   0.972     -.019814     .020525
        AGEQ |    .140915   .0703863     2.00   0.045     .0029605    .2788696
      AGEQSQ |  -.0013605   .0007873    -1.73   0.084    -.0029035    .0001826
       _cons |   .1337523   1.652725     0.08   0.935    -3.105528    3.373033
------------------------------------------------------------------------------
Endogenous: EDUC
Exogenous:  YR20 YR21 YR22 YR23 YR24 YR25 YR26 YR27 YR28 AGEQ AGEQSQ QTR120
            QTR121 QTR122 QTR123 QTR124 QTR125 QTR126 QTR127 QTR128 QTR129
            QTR220 QTR221 QTR222 QTR223 QTR224 QTR225 QTR226 QTR227 QTR228
            QTR229 QTR320 QTR321 QTR322 QTR323 QTR324 QTR325 QTR326 QTR328

. bayesweakiv, reps(200000) discard(20000)
Finished drawing from posterior
Bayesian weak-IV posterior summaries
Credible intervals based on 180000 draws

beta (coefficient of interest)
  Posterior mean:    0.2421
  95% credible interval: [  -0.9971,    1.5647]

delta (endogeneity parameter)
  Posterior mean:   -0.1619
  95% credible interval: [  -1.4845,    1.0771]
```
The table shows the usual output of `ivregress` and is used in Angrist and Krueger (1991).  Calling `bayesweakiv` returns estimates for the returns to schooling (beta) and the level of endogeneity (delta).  For this specification, the posterior mean of 0.2421 nearly doubles the 2SLS estimate of 0.1310.  However, the robust credible interval is much wider, reflecting the issue of weak instruments.  Whereas the typical 95% confidence interval would be [0.0655, 0.1965], the robust credible interval contains this interval, zero, and more.  The credible interval for the endogeneity parameter is wide, offering little insight to the nature of the endogeneity.  The Applications.md file in this repository presents the results for all specifications and cohorts in Angrist and Krueger (1991).

