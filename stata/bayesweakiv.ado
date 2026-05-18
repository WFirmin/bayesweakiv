program define bayesweakiv, eclass
    version 17.0
    syntax [, reps(integer 20000) discard(integer 1000) level(real 95)]

	* Compile Mata functions if not already loaded
    if (0 == findfile("lbayesweakiv.mlib")) {
        quietly do `"`c(sysdir_plus)'b/bayesweakiv.mata"'
    }

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
