program define bayesweakiv, eclass
    version 17.0
    syntax [, reps(integer 20000) discard(integer 1000) level(real 95)]

	/* ---------------------------------------------------------
       1. Verify that last command was a supported IV estimator
    --------------------------------------------------------- */
local cmd "`e(cmd)'"
if !inlist("`cmd'", "ivregress", "ivreg2") {
	di as error "bayesweakiv must follow ivregress or ivreg2"
	exit 198
}

    /* ---------------------------------------------------------
       2. Extract variables — branches by command
    --------------------------------------------------------- */

/* --- Dependent Variable --- */
local depvar `e(depvar)'
/* --- Endogenous regressor --- */
if "`cmd'" == "ivregress" {
	local endogvars `e(endog)'
}
else {
	local endogvars `e(instd)'
}

/* --- Excluded instruments (Z) and controls (C) --- */
if "`cmd'" == "ivregress" {
	local exogr `e(exogr)'
	local exog  `e(exog)'
	local exogr_clean
	foreach v of local exogr {
		if substr("`v'",1,2) != "o." local exogr_clean `exogr_clean' `v'
	}
	local exog_clean
	foreach v of local exog {
		if substr("`v'",1,2) != "o." local exog_clean `exog_clean' `v'
	}
	local Zvars    : list exog_clean - exogr_clean
	local controls `exogr_clean'
	local addcons  = ("`e(constant)'" != "noconstant")
}
else {
	local exexog `e(exexog)'
	local inexog `e(inexog)'
	local exexog_clean
	foreach v of local exexog {
		if substr("`v'",1,2) != "o." local exexog_clean `exexog_clean' `v'
	}
	local inexog_clean
	foreach v of local inexog {
		if substr("`v'",1,2) != "o." local inexog_clean `inexog_clean' `v'
	}
	local Zvars    `exexog_clean'
	local controls `inexog_clean'
	local addcons  = ("`e(noconstant)'" == "")
}

tempvar to_use 
gen `to_use' = e(sample)



/* ---------------------------------------------------------
       3. Call Mata Gibbs sampler
    --------------------------------------------------------- */

	mata: bayesweakiv_draws("`depvar'","`endogvars'","`Zvars'","`controls'", `addcons', `reps', `discard', "`to_use'")
	mata: printf("Finished drawing from posterior\n")

    /* Results returned in Mata global matrices:
         beta_draws
         delta_draws
    */


    /* Posterior summaries */
    mata: st_numscalar("beta_mean", mean(st_matrix("beta_draws")))
    mata: st_numscalar("delta_mean", mean(st_matrix("delta_draws")))

    local alpha = (100 - `level')/2
    mata: st_numscalar("beta_lo", quantile(st_matrix("beta_draws"), `alpha'/100))
    mata: st_numscalar("beta_hi", quantile(st_matrix("beta_draws"), 1 - `alpha'/100))
    mata: st_numscalar("delta_lo", quantile(st_matrix("delta_draws"), `alpha'/100))
    mata: st_numscalar("delta_hi", quantile(st_matrix("delta_draws"), 1 - `alpha'/100))

    /* ---------------------------------------------------------
       4. Display results
    --------------------------------------------------------- */
    di as text "Bayesian weak-IV posterior summaries"
	di as text "Credible intervals based on `=`reps' - `discard'' draws"
    di ""

    di as result "beta (coefficient of interest)"
    di as text    "  Posterior mean: " %9.4f scalar(beta_mean)
    di as text    "  `level'% credible interval: [" ///
        %9.4f scalar(beta_lo) ", " %9.4f scalar(beta_hi) "]"
    di ""

    di as result "delta (endogeneity parameter)"
    di as text    "  Posterior mean: " %9.4f scalar(delta_mean)
    di as text    "  `level'% credible interval: [" ///
        %9.4f scalar(delta_lo) ", " %9.4f scalar(delta_hi) "]"

    /* ---------------------------------------------------------
       5. Return results in e()
    --------------------------------------------------------- */
    ereturn scalar beta_mean = scalar(beta_mean)
    ereturn scalar beta_lo   = scalar(beta_lo)
    ereturn scalar beta_hi   = beta_hi
	ereturn matrix beta_draws = beta_draws
    ereturn scalar delta_mean = scalar(delta_mean)
    ereturn scalar delta_lo   = scalar(delta_lo)
    ereturn scalar delta_hi   = scalar(delta_hi)
	ereturn matrix delta_draws = delta_draws
    ereturn local cmd "bayesweakiv"
end


mata:

real matrix solveFromChol(real matrix L, real matrix v)
{
    return(solveupper(L', solvelower(L,v)))
}
real matrix efficientMvNormalInv(real matrix mu, real matrix cholSigmaInv)
{
	return(solveupper(cholSigmaInv', rnormal(rows(mu),1,0,1)) + mu)
}

/* --------------------------------------------------------------------------
   pop_zscore()
   Population-normalised zscore: subtract mean, divide by sqrt(mean(dev^2)).
   Matches MATLAB's zscore(X, 1) — divides by T, not T-1.
   Works column-by-column for matrices.
-------------------------------------------------------------------------- */
real matrix pop_zscore(real matrix X)
{
    real matrix result, mu
    real scalar j, sd

    mu     = mean(X)
    result = X :- (J(rows(X), 1, 1) * mu)
    for (j = 1; j <= cols(X); j++) {
        sd = sqrt(mean(result[., j]:^2))
        if (sd > 0) {
            result[., j] = result[., j] :/ sd   // standardise normally
        }
        else {
            result[., j] = X[., j]              // constant column — restore original
        }
    }
    return(result)
}


void bayesweakiv_draws(string scalar Yname,
                       string scalar Xname,
                       string scalar Znames,
                       string scalar Cnames,
					   real scalar addcons,
                       real scalar M,
					   real scalar disc,
					   string scalar to_usename)
{
	to_use = st_data(., to_usename)

	y = select(st_data(.,Yname), to_use)
	sy = sqrt(variance(y))
	y = (y :- mean(y)) / sy
	
	x = select(st_data(.,Xname), to_use)
	sx = sqrt(variance(x))
    x = (x :- mean(x)) / sx
	
    z = pop_zscore(select(st_data(.,tokens(Znames)), to_use))
	
	
	//c = st_matrix(Cname)
	c = (Cnames != "" ? select(st_data(., tokens(Cnames)), to_use) : J(rows(y), 0, .))
    if (addcons) c = (c, J(rows(y), 1, 1))
	if (rows(c) > 0 & cols(c) > 0) {
		c = pop_zscore(c) 
	}
	// Dimensions:
	
	T = rows(z)
	k = cols(z)
	l = cols(c)
	
	// Precompute and name certain matrices:
	
	P = makesymmetric((x, c, z, y)'*(x, c, z ,y))
	ww = P[|2,2\l+k+1,l+k+1|]
	wx = P[|2,1\l+k+1,1|]
	wy = P[|2,l+k+2\l+k+1,l+k+2|]
	wc = P[|2,2\l+k+1,l+1|]
	xcz2 = P[|1,1\l+k+1,l+k+1|]
	xczy = P[|1,l+k+2\l+k+1,l+k+2|]
	yy = P[l+k+2,l+k+2]
	cx = P[|2,      1      \ l+1,   1|]        // l x 1
	cz = P[|2,      l+2    \ l+1,   l+k+1|]   // l x k
	cy = P[|2,      l+k+2  \ l+1,   l+k+2|]   // l x 1
	cc = P[|2,      2      \ l+1,   l+1|]      // l x l
	zx = P[|l+2,    1      \ l+k+1, 1|]        // k x 1
	zz = P[|l+2,    l+2    \ l+k+1, l+k+1|]   // k x k
	zy = P[|l+2,    l+k+2  \ l+k+1, l+k+2|]   // k x 1
	zc = P[|l+2,    2      \ l+k+1, l+1|]      // k x l
	// Storage:
	
	BETA = J(M, 1, .)
	DELTA = J(M, 1, .)
	
	// Starting values:
	
	gam2draw = 0.1 
	zzT = (z'*z)/T 
	invVprior = J(l+k,l+k,0)
	invVprior[|l+1,l+1\l+k,l+k|] = zzT 
	pidraw = cholsolve(ww+invVprior/gam2draw, wx)
	uu = (1 \ -pidraw)'*xcz2*(1 \ -pidraw)
	
	// Clean up:
	x = y = z = c = J(0,0,.)
	
	for (i=1; i<=M; i++) {
		
		// Block 1
		R = (I(l+1) \ J(k,l+1,0)), (1 \ -pidraw)
		xcu2 = R'*xcz2*R + diag((0.0001^2, J(1,l,0), 0.0001^2))
		Lxcu2 = cholesky(xcu2)
		bhatols = solveFromChol(Lxcu2, R'*xczy)
		ssr = yy-xczy'*R*bhatols 
		s2epsdraw = 1 / rgamma(1,1,(T-4-l)/2, 2/ssr)
		betadeltadraw = sqrt(s2epsdraw)*solveupper(Lxcu2', rnormal(l+2, 1, 0, 1)) + bhatols
		betadraw = betadeltadraw[1]
		BETA[i] = betadraw 
		deltadraw = betadeltadraw[rows(betadeltadraw)]
		DELTA[i] = deltadraw 
		
		// Block 2
		if (k > 2) {
			s2udraw = 1 / rgamma(1,1,(T+k-2)/2, 2/(uu+pidraw'*invVprior*pidraw/gam2draw))
		} else {
			s2udraw = 1 / rgamma(1,1,T/2, 2/uu)
		}
		
		// Block 3
		if (k > 2) {
			gam2draw = 1 / rgamma(1,1, (k-2)/2, 2/(pidraw[l+1..l+k]'*zzT*pidraw[l+1..l+k]/s2udraw))
		}
		
		// Block 4 
		if (k > 2) {
			alpha = (deltadraw^2*s2udraw+s2epsdraw)/s2epsdraw 
			B = wx-(wy-wx*(betadraw+deltadraw)-wc*betadeltadraw[2..rows(betadeltadraw)-1])*deltadraw*s2udraw/s2epsdraw 
			A = alpha*ww + invVprior/gam2draw 
			LA = cholesky(A) 
			pihat = solveFromChol(LA, B)
			pidraw = sqrt(s2udraw) * solveupper(LA', rnormal(l+k, 1, 0, 1)) + pihat
			uu = (1 \ -pidraw)'*xcz2*(1 \ -pidraw)
		} else {
			aux = (deltadraw^2 * s2udraw + s2epsdraw) / s2epsdraw
			Bxi = cx - cz*pidraw[l+1..l+k] - (cy - cx*(betadraw+deltadraw) - cc*betadeltadraw[2..length(betadeltadraw)-1] + cz*pidraw[l+1..l+k] * deltadraw) * deltadraw * s2udraw/s2epsdraw
			Axi = aux*cc
			LAxi = cholesky(Axi)
			xihat = solveFromChol(LAxi, Bxi)
			xidraw = sqrt(s2udraw) * solveupper(LAxi', rnormal(l,1,0,1)) + xihat
			pidraw[1..l] = xidraw 
			
			Bal = zx - zc*pidraw[1..l] - (zy - zx *(betadraw+deltadraw) - zc * betadeltadraw[2..length(betadeltadraw)-1] + zc * pidraw[1..l] * deltadraw) * deltadraw * s2udraw / s2epsdraw
			Aal = aux*zz
			LAal = cholesky(Aal)
			alhat = solveFromChol(LAal, Bal)
			aldraw = sqrt(s2udraw) * solveupper(LAal', rnormal(k,1,0,1)) + alhat 
			p = min((1, exp(-(k-2)/2*(log(aldraw'*zz*aldraw) - log(pidraw[l+1..l+k]'*zz*pidraw[l+1..l+k])))))
			if (runiform(1,1) < p) {
				pidraw[l+1..l+k] = aldraw
			}
			uu = (1 \ -pidraw)'*xcz2*(1 \ -pidraw)
		}
		
	}
	st_matrix("beta_draws", BETA[(disc+1)..length(BETA)] * sy/sx)
    st_matrix("delta_draws", DELTA[(disc+1)..length(DELTA)] * sy/sx)
}
