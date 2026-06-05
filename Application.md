# Replications
## Summary
The following table describes each application.  The Citation column gives author initials and publication year.  Author names are listed out fully in their section.  The variables Y, X, and Z describe the dependent variable, endogenous regressor, and instrumental variables respectively.  These descriptions leave many specifics in the original papers.  k reports the number of instruments.  Obs. reports the number of observations for cross-sectional data, and when panel data is used, the dimensions of the panel data.  The Data column describes what one observation is.  The Clust. column indicates whether standard errors are clustered.
| Citation | Y | X | Z | k | Obs. | Data | Clust. |
|----------|---|---|---|---|-------|------|------------|
|AK (1991) | Weekly wage | Years of education | (Quarter x Year) of birth dummies | 27 | 250,000-500,000 for each cohort dataset | Census respondents | No |
| ACS (2014) | Economic output | Public investment | Gov. dismissal due to mafia activity | 2 | 95 x 10 | Italian provinces, annual | Yes |
| Y (2014) | TFP growth | Labor-income-share-weighted sum of change in national employment shares | Various macro instruments | 59 | 60 x 24 | Industry, year | No |
| D (2015) | Avg. premium | Ease of manipulating subsidy | Regulation changes | 1 | 34 x 5 | Region, year | Yes |


## Angrist and Krueger (1991)
 
This section replicates selected specifications from:
 
> Angrist, J. D., and A. B. Krueger. 1991. Does compulsory school attendance affect schooling and earnings? *Quarterly Journal of Economics* 106(4): 979–1014.
 
Angrist and Krueger use quarter of birth as an instrument for years of schooling in a wage equation, exploiting the fact that students born earlier in the year are forced to stay in school longer due to compulsory attendance laws. The paper is a landmark application of IV methods and a canonical weak-instruments benchmark.

We download the data and replication files from https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive.  We run the do files with the simple addition of the weakbayesiv command after each 2SLS regression.  We use 200,000 draws of the posterior and discard the first 20,000 for each specification.

 
### Notation
 
Results are reported for specifications drawn from Tables 4, 5, and 6 of the original paper. The notation **T(C)** denotes:
 
- **T** — the table number in Angrist and Krueger (1991)
- **C** — the column number in table T

 
### Results
 
Estimates of the return to schooling (coefficient on years of education). TSLS confidence intervals and Bayesian credible intervals are both 95%. A ⚠️ flag marks specifications where the two methods give notably different intervals, indicating that instrument weakness is influencing the TSLS result.  The Bayes estimate is the mode of the posterior.

| Spec | TSLS | TSLS 95% CI | Bayes | Bayes 95% CI | Note |
|------|-----:|-------------|------:|--------------|------|
4(2) | 0.0769    | [0.0475, 0.1063] | 0.0744 | [0.0357,0.1127] ||
4(4) | 0.1310    | [0.0655, 0.1965] | 0.2074 | [-0.8361,1.2898]|⚠️|
4(6) | 0.0669    | [0.0373, 0.0965] | 0.0658 | [0.0270, 0.1042]||
4(8) | 0.1007    | [0.0352, 0.1662] | 0.4414 | [-1.5846, 2.1307]|⚠️|
5(2) | 0.0891    | [0.0575, 0.1207] | 0.0955 | [0.0569, 0.1356] ||
5(4) | 0.0760    | [0.0192, 0.1328] | 0.0723 | [-0.9355,1.0522] |⚠️|
5(6) | 0.0806    | [0.0485, 0.1127] | 0.0871 | [0.0478,0.1278] ||
5(8) | 0.0600    | [0.0014, 0.1186] | -0.0042| [-0.6460, 0.7589]|⚠️|
6(2) | 0.0553    | [0.0283, 0.0823] | 0.0423 | [0.0017, 0.0819]||
6(4) | 0.0948    | [0.0511, 0.1385] | 0.1385 | [0.0613, 0.2333]||
6(6) | 0.0393    | [0.0109, 0.0677] | 0.0155 | [-0.0293,0.0577]|⚠️|
6(8) | 0.0779    | [0.0311, 0.1247] | 0.1253 | [0.0226,0.2563]||

---
## Acconcia et al. (2014)

This paper estimates the output multiplier from spending cuts.  We report the original results from the first TSLS specification in Table 4.  Since the Bayesian method does not allow clustering of standard errors, we also report TSLS confidence intervals without clustering.

| Method | Estimate | 95% CI |
|--------|----------|--------|
|TSLS | 1.4565 | [0.4922,2.4208]|
|TSLS (No Clustering) |1.4565 | [0.5220,2.3910]|
|Bayesian | 1.5261 | [0.6089,2.7551] |

## Young (2014)
This paper considers many instruments, but selects four sets of instruments with an F test and runs TSLS using each set.  Of the four sets, only one provides a statistically significant estimate.  Three estimates are negative while one (statistically insignificant) estimate is positive.  The bayesian method has the advantage of providing valid inference without a selection step.

## Decarolis (2015)
