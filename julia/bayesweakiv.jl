using Random, Distributions, Statistics, StatsBase;
using PDMats, LinearAlgebra, StatsFuns;

##############################################################################
#  Bayesian IV estimation robust to weak instruments
#
#  Model:
#    X = Z * π  + C * ρ  + ν          (first stage)
#    Y = X * β  + C * α  + δ * ν + ε  (structural equation)
#
#  Reference:
#    Giannone, D., Lenza, M., and Primiceri, G. (2024).
#    "Bayesian Inference in IV Regressions"
##############################################################################


"""
    drawData(; T, beta, delta, k, sigma2_nu, sigma2_eps, l, alpha, rho) -> (Y, X, Z, C)

Simulate data from the triangular IV model used by Giannone, Lenza, and
Primiceri (2024).

# Keyword arguments
- `T::Int = 250`: number of observations.
- `beta::Real = 0`: true structural coefficient on the endogenous regressor.
- `delta::Real = 0.75`: true endogeneity parameter (loading of the first-stage
  error `ν` on the structural equation).
- `k::Int = 10`: number of excluded instruments. The first-stage coefficients
  `π` are drawn as `π ~ N(0, s² I)` where `s ~ Uniform(0, 0.25)`, producing
  weak-to-moderate instrument strength by default.
- `sigma2_nu::Real = 1`: variance of the first-stage error `ν`.
- `sigma2_eps::Real = 1`: variance of the structural error `ε`.
- `l::Int = 5`: number of exogenous controls in `C`.
- `alpha::Vector = zeros(l)`: true coefficients on `C` in the structural equation.
- `rho::Vector = zeros(l)`: true coefficients on `C` in the first stage.

# Returns
- `Y`: T × 1 outcome variable.
- `X`: T × 1 endogenous regressor.
- `Z`: T × k instrument matrix (columns drawn i.i.d. from N(0, I)).
- `C`: T × l control matrix (columns drawn i.i.d. from N(0, I)).

# Example
```julia
Y, X, Z, C = drawData(T=500, beta=1.0, delta=0.5, k=3)
```
"""
function drawData(; T=250, beta=0, delta=0.75, k=10, sigma2_nu=1, sigma2_eps=1, l=5, alpha=zeros(l), rho=zeros(l))
    Z = rand(MvNormal(zeros(k),I),T)'
    C = rand(MvNormal(zeros(l),I),T)'

    # First Stage:
    s = rand(Uniform(0,0.25))
    pie = rand(MvNormal(zeros(k),s^2*I))
    nu = rand(Normal(0,sigma2_nu),T)
    X = Z*pie + C*rho + nu

    # Second Stage:
    eps = rand(Normal(0,sigma2_eps),T)
    Y = X*beta + C*alpha + delta*nu + eps

    return Y,X,Z,C
end

"""
    forceSym(A) -> Matrix

Return `(A + A') / 2`, enforcing exact symmetry after floating-point
accumulation. Use before passing a nearly-symmetric matrix to `cholesky`.
"""
forceSym(A) = (A + A')/2

"""
    cholInv(cholA, v) -> Matrix

Solve the linear system `A \\ v` given the lower-triangular Cholesky factor
`cholA` of `A` (i.e. `A = cholA * cholA'`).

Avoids re-factorising `A` when the same factorisation is reused across
multiple right-hand sides.
"""
cholInv(cholA, v) = cholA' \ (cholA \ v)

"""
    efficientMvNormal(mu, cholSigma) -> Vector

Draw from N(`mu`, `Sigma`) given the lower-triangular Cholesky factor `cholSigma`
of the covariance matrix `Sigma` (i.e. `Sigma = cholSigma * cholSigma'`).

Avoids forming `Sigma` explicitly: draws `z ~ N(0, I)` and returns
`cholSigma * z + mu`.
"""
efficientMvNormal(mu, cholSigma) = cholSigma * rand(MvNormal(zeros(size(cholSigma,1)), I)) + mu

"""
    efficientMvNormalInv(mu, cholSigmaInv) -> Vector

Draw from N(`mu`, `Sigma`) given the lower-triangular Cholesky factor
`cholSigmaInv` of the **precision** matrix `Sigma⁻¹`
(i.e. `Sigma⁻¹ = cholSigmaInv * cholSigmaInv'`).

Avoids forming `Sigma` explicitly: draws `z ~ N(0, I)` and returns
`cholSigmaInv' \\ z + mu`. This is the form used in the Gibbs sampler
because the precision matrix arises naturally from the posterior conditionals.
"""
efficientMvNormalInv(mu, cholSigmaInv) = cholSigmaInv' \ rand(MvNormal(zeros(size(cholSigmaInv,1)), I)) + mu

"""
    drawPost(Y, X, Z, C, M) -> (betas, deltas)

Gibbs sampler for the Bayesian IV model of Giannone, Lenza, and Primiceri (2024).

Samples from the joint posterior of the structural parameters `(β, α, δ)` and
the first-stage parameters `(π, ρ)` via the following four-block scheme:

1. **β, α, δ, σ²_ε** — multivariate normal / inverse-gamma conditional,
   given the current residual `ν̂ = X − Zπ − Cρ` as a regressor.
2. **σ²_ν** — inverse-gamma conditional on the current `ν̂` and `π`.
3. **γ²** — inverse-gamma hyperprior on the scale of the instrument prior
   (sampled only when `k > 2`).
4. **π, ρ** — multivariate normal conditional. When `k > 2`, `π` and `ρ`
   are drawn jointly. When `k ≤ 2`, they are drawn sequentially; a
   Metropolis–Hastings accept/reject step is applied to `π` when `k > 1`
   to account for the non-standard geometry of the weak-instrument region.

The prior on the first-stage coefficients is
`π | γ², σ²_ν ~ N(0, γ² σ²_ν Ω⁻¹)` where `Ω = Z'Z / T`.
This shrinks `π` toward zero adaptively, making the posterior for `β` proper
even when instruments are arbitrarily weak.

# Arguments
- `Y`: T × 1 outcome variable.
- `X`: T × 1 endogenous regressor.
- `Z`: T × k excluded instruments.
- `C`: T × l exogenous controls (include a column of ones for a constant).
- `M`: total number of Gibbs draws to collect.

# Returns
- `betas`:  length-`M` vector of posterior draws for β.
- `deltas`: length-`M` vector of posterior draws for δ.

Discard an initial burn-in segment before computing posterior summaries,
e.g. `betas[1001:end]`.

# Example
```julia
Y, X, Z, C = drawData(T=250, beta=0.0, delta=1.0, k=2)
betas, deltas = drawPost(Y, X, Z, C, 20_000)

burn = 2000
using Statistics
println("Posterior mean β: ", mean(betas[burn+1:end]))
println("Posterior mean δ: ", mean(deltas[burn+1:end]))
```
"""
function drawPost(Y,X,Z,C,M)
    # pie | gamma2, sigma2_nu ~ N(0,gamma2*sigma2_nu*Omega)
    
    T,l = size(C)
    k = size(Z,2)
    regularize=0.0000001

    # initial calculations:
    invOmega = Z'Z / T
    #if k == 1; invOmega = fill(invOmega, 1,1); end;
    cross = [X Z C Y]'*[X Z C Y]
    P = zeros(l+k,l+k)
    P[1:k,1:k] = invOmega

    # initial values: use TSLS or some random defaults (like gamma2=0.1)
    gamma2 = 0.1
    coef0First = cholInv(cholesky(forceSym(cross[2:end-1,2:end-1] + T*P/gamma2)).L, cross[2:end-1,1])
    pie, rho = coef0First[1:k], coef0First[k+1:end]

    betas = []
    deltas = []

    
    xTilde = [X C X - Z*pie - C*rho]
    xTilde2 = xTilde'xTilde + regularize*I

    # initialize
    for m in 1:M
        # block 1: beta, alpha, delta, sigma2_eps
        xTilde[:,end] = X - Z*pie - C*rho # only updates last column
        xTilde2[:,end] = xTilde'xTilde[:,end] # this and next line only update last col and row
        xTilde2[end,:] = xTilde[:,end]'xTilde
        xTilde2[end,end] += regularize

        xTilde2Chol = cholesky(xTilde2).L
        xiHat = cholInv(xTilde2Chol, xTilde'Y)
        sHat = sum((Y - xTilde*xiHat).^2)
        
        sigma2_eps = rand(InverseGamma((T-l-4)/2, sHat/2))
        coef1 = efficientMvNormalInv(xiHat, xTilde2Chol/sqrt(sigma2_eps))
        beta, alpha, delta = coef1[1], coef1[2:end-1], coef1[end]

        # block 2: sigma2_nu
        if k > 2
            sigma2_nu = rand(InverseGamma((T+k-2)/2, (xTilde2[end,end] - regularize)/2 + pie'invOmega*pie/(2*gamma2)))
        else 
            sigma2_nu = rand(InverseGamma(T/2, (xTilde2[end,end] - regularize)/2))
        end

        # block 4: gamma2
        # this is the fourth block in the original paper 
        # we can move it to the third position without consequence
        if k > 2
            gamma2 = rand(InverseGamma((k-2)/2, pie'invOmega*pie/(2*sigma2_nu)))
        end

        # block 3: pi, rho
        if k > 2
            cholA = cholesky(forceSym((1+sigma2_nu/sigma2_eps*delta^2)*cross[2:end-1,2:end-1] + P/gamma2)).L
            omegaHat = cholInv(cholA, cross[2:end-1, 1]*(1 + sigma2_nu/sigma2_eps*delta*(beta+delta)) - sigma2_nu/sigma2_eps*delta*(cross[2:end-1,end] - cross[2:end-1,k+2:end-1]*alpha))
            omega = efficientMvNormalInv(omegaHat, cholA / sqrt(sigma2_nu))
            pie, rho = omega[1:k], omega[k+1:end]
        else 
            # rho:
            cholA = cholesky((1+sigma2_nu/sigma2_eps*delta^2)*cross[k+2:end-1,k+2:end-1]).L
            rhoHat = cholInv(cholA, 
                            cross[k+2:end-1,1]*(1 + sigma2_nu/sigma2_eps*delta*(beta+delta)) 
                            - cross[k+2:end-1,2:k+1]*pie*(1 + sigma2_nu/sigma2_eps*delta^2)
                            - sigma2_nu/sigma2_eps*delta*(cross[k+2:end-1,end] - cross[k+2:end-1,k+2:end-1]*alpha))
            rho = efficientMvNormalInv(rhoHat, cholA / sqrt(sigma2_nu))

            # pi:
            cholA = cholesky((1+sigma2_nu/sigma2_eps*delta^2)*cross[2:k+1,2:k+1]).L
            piBar = cholInv(cholA, 
                            cross[2:k+1,1]*(1 + sigma2_nu/sigma2_eps*delta*(beta+delta))
                            - cross[2:k+1,k+2:end-1]*(rho + sigma2_nu/sigma2_eps*delta*(rho*delta-alpha))
                            -sigma2_nu/sigma2_eps*delta*cross[2:k+1,end])
            pieCandidate = efficientMvNormalInv(piBar, cholA / sqrt(sigma2_nu))
            if k == 2
                pie = pieCandidate 
            else 
                p = min(1, abs(pieCandidate[1])/abs(pie[1]))
                pie = rand(Bernoulli(p)) ? pieCandidate : pie
            end
        end

       
        # save
        push!(betas, beta)
        push!(deltas, delta)
    end
    return betas, deltas
end
