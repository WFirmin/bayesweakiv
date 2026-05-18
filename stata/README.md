# bayesweakiv — Stata
 
Stata implementation of the Bayesian IV estimator from:
 
> Giannone, D., Lenza, M., and Primiceri, G. E. (2024). *Bayesian Inference in IV Regressions with Possible Weak Instruments*. Working paper.
 
---
 
## Requirements
 
- Stata 17 or later
- No additional packages required
---
 
## Installation
 
```stata
net install bayesweakiv, from("https://raw.githubusercontent.com/WFirmin/bayesweakiv/main/stata")
```
 
To update an existing installation:
 
```stata
ado uninstall bayesweakiv
net install bayesweakiv, from("https://raw.githubusercontent.com/WFirmin/bayesweakiv/main/stata")
```
 
---
 
## Usage
 
`bayesweakiv` is a post-estimation command that must immediately follow `ivregress`:
 
```stata
ivregress 2sls depvar [controls] (endog = instruments) [if] [in]
bayesweakiv [, reps(#) discard(#) level(#)]
```
 
### Options
 
| Option | Default | Description |
|--------|---------|-------------|
| `reps(#)` | 20000 | Total number of Gibbs draws |
| `discard(#)` | 1000 | Number of burn-in draws to discard |
| `level(#)` | 95 | Credible interval level (%) |
 
### Saved results (`e()`)
 
| Scalar | Description |
|--------|-------------|
| `e(beta_mean)` | Posterior mean of β |
| `e(beta_lo)` | Lower bound of credible interval for β |
| `e(beta_hi)` | Upper bound of credible interval for β |
| `e(delta_mean)` | Posterior mean of δ (endogeneity parameter) |
| `e(delta_lo)` | Lower bound of credible interval for δ |
| `e(delta_hi)` | Upper bound of credible interval for δ |
 
The full posterior draw vectors `beta_draws` and `delta_draws` are also left
in memory as Stata matrices after the command runs, for custom post-processing
such as plotting posterior densities.
 
---
 
## Example
 
```stata
sysuse auto
ivregress 2sls price (mpg = weight length)
bayesweakiv, reps(20000) discard(1000) level(95)
 
* Access saved results
di "Posterior mean of beta: " %8.4f e(beta_mean)
di "95% credible interval:  [" %8.4f e(beta_lo) ", " %8.4f e(beta_hi) "]"
 
* Plot the posterior density
preserve
    drop _all
    svmat beta_draws, names(beta)
    kdensity beta1, title("Posterior density of {&beta}")
restore
```
 
---
 
## Limitations
 
- Supports a **single endogenous regressor** only.
- A constant is included automatically unless `noconstant` was specified in
  the preceding `ivregress` call.
- Requires Stata 17 or later due to use of Mata features introduced in that version.
---
 
## Citation
 
```
Giannone, D., Lenza, M., and Primiceri, G. E. (2024). Bayesian Inference in IV
Regressions with Possible Weak Instruments. Working paper.
 
Firmin, W. (2025). bayesweakiv: Bayesian IV Estimation Robust to Weak
Instruments [Stata and Julia]. https://github.com/WFirmin/bayesweakiv
```
