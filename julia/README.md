# bayesweakiv — Julia
 
Julia implementation of the Bayesian IV estimator from:
 
> Giannone, D., Lenza, M., and Primiceri, G. E. (2024). *Bayesian Inference in IV Regressions with Possible Weak Instruments*. Working paper.
 
---
 
## Requirements
 
- Julia 1.9 or later
- The following packages (install once with `Pkg.add`):
```julia
using Pkg
Pkg.add(["Distributions", "LinearAlgebra", "Statistics", "StatsBase", "PDMats", "StatsFuns"])
```
 
---
 
## Installation
 
This package is not yet registered in the Julia General registry. To use it,
download `bayesweakiv.jl` and load it directly:
 
```julia
include("bayesweakiv.jl")
```
 
---
 
## Usage
 
### Simulate data
 
```julia
Y, X, Z, C = drawData(T=250, beta=0.0, delta=0.75, k=5)
```
 
`C` should include a column of ones if you want a constant in the model.
`drawData` does not add one automatically.
 
### Run the sampler
 
```julia
betas, deltas = drawPost(Y, X, Z, C, 20_000)
```
 
Discard the first portion of draws as burn-in before computing summaries:
 
```julia
using Statistics
burn = 2000
println("Posterior mean β: ", mean(betas[burn+1:end]))
println("Posterior mean δ: ", mean(deltas[burn+1:end]))
println("95% CI for β:     ", quantile(betas[burn+1:end], [0.025, 0.975]))
```
 
---
 
## Function reference
 
### `drawData`
 
```julia
drawData(; T, beta, delta, k, sigma2_nu, sigma2_eps, l, alpha, rho) -> (Y, X, Z, C)
```
 
| Argument | Default | Description |
|----------|---------|-------------|
| `T` | 250 | Number of observations |
| `beta` | 0 | True structural coefficient |
| `delta` | 0.75 | True endogeneity parameter |
| `k` | 10 | Number of excluded instruments |
| `sigma2_nu` | 1 | Variance of first-stage error |
| `sigma2_eps` | 1 | Variance of structural error |
| `l` | 5 | Number of exogenous controls |
| `alpha` | `zeros(l)` | True coefficients on controls in structural equation |
| `rho` | `zeros(l)` | True coefficients on controls in first stage |
 
First-stage coefficients `π` are drawn as `π ~ N(0, s²I)` with
`s ~ Uniform(0, 0.25)`, giving weak-to-moderate instrument strength by default.
 
### `drawPost`
 
```julia
drawPost(Y, X, Z, C, M) -> (betas, deltas)
```
 
| Argument | Type | Description |
|----------|------|-------------|
| `Y` | T × 1 Matrix | Outcome variable |
| `X` | T × 1 Matrix | Endogenous regressor |
| `Z` | T × k Matrix | Excluded instruments |
| `C` | T × l Matrix | Exogenous controls (include a ones column for a constant) |
| `M` | Int | Number of Gibbs draws |
 
Returns length-`M` vectors `betas` and `deltas` of posterior draws for β and δ.
 
---
 
## A note on instrument count (`k`)
 
The sampler behaves differently depending on the number of excluded instruments:
 
- **`k > 2`**: `π` and `ρ` are drawn jointly in a single multivariate normal
  step, and the hyperparameter `γ²` is also sampled.
- **`k ≤ 2`**: `π` and `ρ` are drawn sequentially. When `k > 1`, a
  Metropolis–Hastings accept/reject step is applied to the draw of `π` to
  account for the geometry of the weak-instrument region.
In both cases the posterior for β is proper regardless of instrument strength.
 
---
 
## Citation
