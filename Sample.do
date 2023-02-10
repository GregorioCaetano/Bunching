qui {
*************************************************************************************
* This is the only piece of code you should need to change
*************************************************************************************
local path = "/Users/gregoriocaetano/Dropbox/My Mac (TCB699021)/Downloads/Code/"	// add your path here
local outcome = "Y"						// outcome variable
local treatment = "X"					// treatment variable
local controls = "ChildMale ChildWhite ChildBlack ChildHispanic"			// list of controls
local boot_reps = 100					// number of bootstrap iterations
local K = 5 							// number of clusters of controls for the estimation of the expectation, and number of cluster indicators for the specification of controls
local Kd = 1							// number of clusters of Z for \delta(Z). 
local method = "symmetric" 				// correction method, see below for options:
* "naiveN" (uncorrected, no controls)
* "naive" (uncorrected, with controls including `K' cluster indicators)
* "tobit" (Tobit correction, with controls including `K' cluster indicators)
* "het_tobit" (semiparametric Tobit correction, with controls including `K' cluster indicators)
* "symmetric" (tail symmetry correction, with controls including `K' cluster indicators)
*************************************************************************************
* Housekeeping macros
*************************************************************************************
set seed 1
capture program drop ccn_dist
cd "`path'"
*************************************************************************************
* Create clusters to discretize Zs, denoted Z_`K'. Skip this step if dataset already includes them.
*************************************************************************************
local maxK = max(`K',`Kd')
use "SampleData.dta", clear
noisily display "Creating clusters indicators..."
cluster wardslinkage `controls', measure(Gower)  		// hierarquical clustering on original Zs
forvalues k=1(1)`maxK' {
	cluster generate Z_`k' = groups(`k'), ties(more)
}
noisily display "clusters done"
keep `outcome' `treatment' `controls' Z_* 
order `outcome' `treatment' `controls' Z_* 
save "SampleDataWithClusters.dta", replace
*************************************************************************************
*************************************************************************************
*************************************************************************************
*************************************************************************************
* Beginning estimation...
** Get point estimates
global boot = 0
use "SampleDataWithClusters.dta", clear		
ccn_dist `method' `K' `Kd' `outcome' `treatment' "`controls'"
local beta = round(`r(BETA)',0.001)
local delta = round(`r(DELTA)',0.001)
local F = round(`r(F)',0.001)
*** need the point estimates of each \delta in order to report correct p-value of the F-test under bootstrap, see below
if `Kd'>1 {
	foreach type in `method' {
		forvalues kk=1(1)`Kd' {
			global d`kk'_`type' = r(DELTA_`kk')
		}
	}
}
** Get bootstrapped standard errors 
global boot = 1
matrix actual_estimates_mat = J(`boot_reps', 2,. )
forvalues i=1(1)`boot_reps' {
	noisily display "Iteration `i'"
 	use "SampleDataWithClusters.dta", clear
	bsample 
	ccn_dist `method' `K' `Kd' `outcome' `treatment' "`controls'"
	matrix actual_estimates_mat[`i',1] = r(BETA)
	if `Kd'==1 {
		matrix actual_estimates_mat[`i',2] = r(DELTA)	
	}
	else {
		matrix actual_estimates_mat[`i',2] = r(F)
	}
}
*** If Kd=1, then report \delta(SE). Otherwise, report F-stat(p-value) of the test of whether \delta(Z)=0 for all clusters of Z.
if `Kd'==1 {
	mat colnames actual_estimates_mat = beta delta
	clear
	svmat actual_estimates_mat, names(col)
	sum beta
	local se_beta = round(`r(sd)',0.001)
	sum delta
	local se_delta = round(`r(sd)',0.001)
	** Report estimates
	noisily display "beta(SE) = `beta'(`se_beta')"
	if "`method'"=="tobit" | "`method'"=="het_tobit" | "`method'"=="symmetric" {
		noisily display "delta(SE) = `delta'(`se_delta')"
	}
}
else {
	mat colnames actual_estimates_mat = beta F
	clear
	svmat actual_estimates_mat, names(col)
	sum beta
	local se_beta = round(`r(sd)',0.001)
	gen diff_F=abs(F - `F')
	xtile xtile_F = F, nq(`boot_reps' )
	sort diff_F
	gen pvalue_F = 1 - (xtile_F[1]/(`boot_reps' ))
	sum pvalue_F
	local pvalue_delta = round(`r(mean)',0.001)
	** Report estimates
	noisily display "beta(SE) = `beta'(`se_beta')"
	if "`method'"=="tobit" | "`method'"=="het_tobit" | "`method'"=="symmetric" {
		noisily display "F(pvalue) = `F'(`pvalue_delta')"
	}
}


}
