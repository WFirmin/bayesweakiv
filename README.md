# bayesweakiv

**Bayesian inference for instrumental variables (IV) regression, robust to weak instruments.**

This repository provides Stata and Julia implementations of the Bayesian IV estimator developed in:

> Giannone, D., Lenza, M., and Primiceri, G., 2026. "Bayesian Inference in IV Regressions," Working Paper Series 3189, European Central Bank.

The estimator places a prior on the first-stage coefficients that shrinks toward zero and is calibrated to the data, making it robust to weak instruments without requiring the practitioner to pre-test or choose a threshold.

---

## Model

The estimator covers the triangular IV model:

```
X = Z * π + C * ρ + ν          (first stage)
Y = X * β + C * α + δ * ν + ε  (structural equation)
```

where:
- `Y` is the outcome variable
- `X` is the endogenous regressor
- `Z` is the matrix of excluded instruments
- `C` contains exogenous controls (including a constant)
- `β` is the structural coefficient of interest
- `δ` captures endogeneity (correlation between the first-stage error `ν` and the structural error `ε`)

The key parameter `δ` is jointly estimated with `β`, allowing the posterior for `β` to automatically account for the degree of instrument strength.

---

## Repository structure

```
bayesweakiv/
├── README.md           ← this file
├── LICENSE             ← MIT
├── CITATION.cff        ← structured citation
├── stata/
│   ├── README.md       ← Stata-specific instructions
│   ├── bayesweakiv.ado ← main command and sampler implementation
│   ├── bayesweakiv.pkg ← package manifest
│   ├── stata.toc       ← stata table of contents
│   └── example.do      ← worked example (planned)
└── julia/
    ├── README.md       ← Julia-specific instructions
    ├── bayesweakiv.jl  ← sampler implementation
    └── example.jl      ← worked example (planned)
```

---

## Quick start

### Stata

```stata
* Install (once)
net install bayesweakiv, from("https://raw.githubusercontent.com/WFirmin/bayesweakiv/main/stata")

* Run ivregress, then call bayesweakiv
sysuse auto
ivregress 2sls price (mpg = weight length)
bayesweakiv, reps(20000) discard(1000) level(95)
```

### Julia

This package is not yet registered in the Julia General registry. To use it,
download `bayesweakiv.jl` and load it directly:

```julia
include("bayesweakiv.jl")
betas, deltas = drawPost(Y, X, Z, C, 20000)
```

See the `stata/` and `julia/` subdirectories for full worked examples.

---

## Citation

If you use this software, please cite both the paper and the software:

**Paper:**
```
Giannone, D., Lenza, M., and Primiceri, G., 2026. "Bayesian Inference in IV Regressions," Working Paper Series 3189, European Central Bank.
```

**Software:**
```
[Names] (2026). bayesweakiv: Bayesian IV Estimation Robust to Weak Instruments
[Stata and Julia]. https://github.com/WFirmin/bayesweakiv
```

---

## License

MIT — see `LICENSE`.
