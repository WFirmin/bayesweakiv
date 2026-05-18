program define bayesweakiv, eclass
    version 17.0
    syntax [, reps(integer 20000) discard(integer 1000) level(real 95)]

    /* ---------------------------------------------------------
       1. Verify that last command was ivregress
    --------------------------------------------------------- */
    if "`e(cmd)'" != "ivregress" {
        di as error "bayesweakiv must follow ivregress"
        exit 198
    }

    /* ---------------------------------------------------------
       2. Extract variables from last ivregress call
    --------------------------------------------------------- */
	
    tempname Y X Z C
	tempvar exogr exog exogr_clean exog_clean Zvars

    /* Dependent variable */
    mata: st_matrix("Y", st_data(., "`e(depvar)'"))

    /* Endogenous regressors */
    mata: st_matrix("X", st_data(., "`e(endog)'"))
	
	/* Exogenous variables */
	local exogr `e(exogr)'
	local exog `e(exog)'
	/* Remove omitted variables (those starting with o.) */
	local exogr_clean
	foreach v of local exogr {
		if substr("`v'",1,2) != "o." {
			local exogr_clean `exogr_clean' `v'
		}
	}
	local exog_clean
	foreach v of local exog {
		if substr("`v'",1,2) != "o." {
			local exog_clean `exog_clean' `v'
		}
	}
	local Zvars : list exog_clean - exogr_clean
	mata: st_matrix("Z", st_data(., "`Zvars'"))
	if "`exogr_clean'" != "" {
		mata: st_matrix("C", st_data(., "`exogr_clean'"))
	}
	else {
		mata: st_matrix("C", J(rows(st_matrix("X")),0,.))
	}
	if "`e(constant)'" != "noconstant" {
		mata: st_matrix("C", (st_matrix("C"), J(rows(st_matrix("C")),1,1)))
	}

    /* ---------------------------------------------------------
       3. Call Mata Gibbs sampler
    --------------------------------------------------------- */
    mata: bayesweakiv_draws("Y","X","Z","C", `reps', `discard')

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
    di as text "Credible intervals based on `reps' draws"
    di ""

    di as result "beta (coefficient of interest)"
    di as text    "  Posterior mean: " %9.4f beta_mean
    di as text    "  `level'% credible interval: [" ///
        %9.4f beta_lo ", " %9.4f beta_hi "]"
    di ""

    di as result "delta (endogeneity parameter)"
    di as text    "  Posterior mean: " %9.4f delta_mean
    di as text    "  `level'% credible interval: [" ///
        %9.4f delta_lo ", " %9.4f delta_hi "]"
		
    /* ---------------------------------------------------------
       5. Return results in e()
    --------------------------------------------------------- */
    ereturn clear
    ereturn scalar beta_mean = beta_mean
    ereturn scalar beta_lo   = beta_lo
    ereturn scalar beta_hi   = beta_hi
    ereturn scalar delta_mean = delta_mean
    ereturn scalar delta_lo   = delta_lo
    ereturn scalar delta_hi   = delta_hi
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

void bayesweakiv_draws(string scalar Yname,
                       string scalar Xname,
                       string scalar Zname,
                       string scalar Cname,
                       real scalar M,
					   real scalar disc)
{
    Y = st_matrix(Yname)
    X = st_matrix(Xname)
    Z = st_matrix(Zname)
    C = st_matrix(Cname)

    T = rows(Y)
    k = cols(Z)
    l = cols(C)
	reg = 1e-7

    /* Initial Calculations */
	invOmega = (Z'Z) / T
	cross = makesymmetric((X, Z, C, Y)' * (X, Z, C, Y))
	P = J(k+l, k+l, 0)
    P[1..k,1..k] = invOmega
    
    /* Initial values */
    gamma2 = 0.1
	coef0 = cholsolve(cross[|2,2\k+l+1,k+l+1|] + T*P/gamma2, cross[|2,1\k+l+1,1|])
    
    pi = coef0[1..k]
    rho = coef0[(k+1)..(k+l)]

    /* Storage */
    beta_draws  = J(M,1,.)
    delta_draws = J(M,1,.)

    /* Construct xTilde */
    xTilde = (X, C, X - Z*pi - C*rho)
	xTilde2 = xTilde' * xTilde + reg*I(l+2)
	

    for (m=1; m<=M; m++) {

        /* -------------------------
           Block 1: beta, alpha, delta
        ------------------------- */
        xTilde[.,cols(xTilde)] = X - Z*pi - C*rho
		xTilde2[.,cols(xTilde2)] = xTilde' * xTilde[.,cols(xTilde)]
		xTilde2[cols(xTilde2),.] = xTilde2[.,cols(xTilde2)]'
		xTilde2[rows(xTilde2),cols(xTilde2)] = xTilde2[rows(xTilde2),cols(xTilde2)] + reg
		
		xTilde2Chol = cholesky(xTilde2)
		xiHat = solveFromChol(xTilde2Chol, xTilde'Y)
        resid = Y - xTilde*xiHat
        sHat = quadcross(resid,resid)

        sigma2_eps = 1 / rgamma(1,1,(T-l-4)/2, (sHat/2)^(-1))
		coef1 = efficientMvNormalInv(xiHat, xTilde2Chol/sqrt(sigma2_eps))
		

        beta  = coef1[1]
        alpha = coef1[|2\l+1|]
        delta = coef1[l+2]

        /* -------------------------
           Block 2: sigma2_nu
        ------------------------- */
        if (k > 2) {
			piInvOmegaPi = quadcross(pi, invOmega*pi)
			sigma2_nu = 1/rgamma(1,1,(T+k-2)/2, 
				2/((xTilde2[cols(xTilde),cols(xTilde)]-reg)
                + piInvOmegaPi/gamma2))
        }
        else {
			sigma2_nu = 1/rgamma(1,1,T/2,
				2/(xTilde2[cols(xTilde),cols(xTilde)]-reg))
        }

        /* -------------------------
           Block 4: gamma2
        ------------------------- */
        if (k > 2) {
			gamma2 = 1/rgamma(1,1,(k-2)/2,
				2*sigma2_nu/piInvOmegaPi)
        }

        /* -------------------------
           Block 3: pi, rho
        ------------------------- */
		if (k > 2) {
			cholA = cholesky((1 + sigma2_nu/sigma2_eps * delta^2) * cross[|2,2\k+l+1,k+l+1|] + P/gamma2)
			omegaHat = solveFromChol(cholA, cross[|2,1\k+l+1,1|] * (1 + sigma2_nu/sigma2_eps*delta*(beta+delta)) ///
				- sigma2_nu/sigma2_eps*delta * (cross[|2,cols(cross)\k+l+1,cols(cross)|] - cross[|2,k+2\k+l+1,k+l+1|]*alpha))
			omega = efficientMvNormalInv(omegaHat, cholA/sqrt(sigma2_nu))
			pi = omega[|1\k|]
			rho = omega[|k+1\k+l|]
		}
		else {
			cholA = cholesky((1 + sigma2_nu/sigma2_eps * delta^2) * cross[|k+2,k+2\k+l+1,k+l+1|])
			rhoHat = solveFromChol(cholA, cross[|k+2,1\k+l+1,1|]*(1 + sigma2_nu/sigma2_eps*delta*(beta+delta)) ///
				- cross[|k+2,2\k+l+1,k+1|]*pi*(1 + sigma2_nu/sigma2_eps*delta^2) ///
				- sigma2_nu/sigma2_eps*delta*(cross[|k+2,k+l+2\k+l+1,k+l+2|] - cross[|k+2,k+2\k+l+1,k+l+1|]*alpha))
			rho = efficientMvNormalInv(rhoHat, cholA/sqrt(sigma2_nu))
			
			cholA = cholesky((1+sigma2_nu/sigma2_eps*delta^2)*cross[|2,2\k+1,k+1|])
			piBar = solveFromChol(cholA, cross[|2,1\k+1,1|]*(1 + sigma2_nu/sigma2_eps*delta*(beta+delta)) ///
				- cross[|2,k+2\k+1,k+l+1|]*(rho + sigma2_nu/sigma2_eps*delta*(rho*delta-alpha)) ///
				- sigma2_nu/sigma2_eps*delta*cross[|2,k+l+2\k+1,k+l+2|])
			piCandidate = efficientMvNormalInv(piBar, cholA/sqrt(sigma2_nu))
			if (k == 2) {
				pi = piCandidate 
			}
			else {
				p = min((1, abs(piCandidate)/abs(pi)))
				if (runiform(1,1) <= p) {
					pi = piCandidate 
				}
			}
		}
		

        /* Save draws */
        beta_draws[m]  = beta
        delta_draws[m] = delta
    }

    st_matrix("beta_draws", beta_draws[(disc+1)..length(beta_draws)])
    st_matrix("delta_draws", delta_draws[(disc+1)..length(beta_draws)])
}


void draw_data(real scalar delta,
				real scalar k)
{
	Z = rnormal(250,k,0,1)
	s = runiform(1,1)/4
	nu = rnormal(250,1,0,1)
	st_matrix("X", Z*rnormal(k,1,0,s) + nu)
	st_matrix("Y", delta*nu + rnormal(250,1,0,1))
	st_matrix("Z", Z)
	st_matrix("C", J(250,1,1))
}

void addFstat(real matrix Z,
			real matrix C,
			real matrix X,
			real matrix Fstats)
{
	piRho = cholsolve((Z, C)'*(Z, C), (Z, C)'*X )
	piRhoR = cholsolve(C'*C, C'*X )
	SSRu = sum((X - (Z, C)*piRho):^2)
	SSRr = sum((X - C*piRhoR):^2)
	F = (SSRr-SSRu)/cols(Z) / (SSRu/(rows(X)-cols(Z)-cols(C)))
	st_matrix("Fstats", Fstats\F)
	
}

void update_CI(real matrix betas,
				real matrix uppers,
				real matrix lowers,
				real scalar level)
{
	st_matrix("uppers", uppers \ quantile(betas, 1 - (1 - level/100)/2))
	st_matrix("lowers", lowers \ quantile(betas, (1 - level/100)/2))
}

void power(real matrix domain,
			real matrix uppers,
			real matrix lowers,
			real matrix Fstats)
{
	fVals = (0,2,4,6,8,10,20)
	counts = J(length(domain),6,0)
	for(i=1;i<=6;i++){
		mask = (Fstats :>= fVals[i]) :& (Fstats :< fVals[i+1])
		U = select(uppers, mask)
		L = select(lowers, mask)
		for(j=1;j<=length(domain);j++){
			counts[j,i] = sum((L :<= domain[j]) :& (U :>= domain[j]))
		}
		counts[.,i] = counts[.,i] / length(U)
	}
	st_matrix("counts",counts)
	
	
}
end

