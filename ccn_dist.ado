qui {
 capture program drop ccn_dist
 noisily display "`model_type'..."
 program define ccn_dist, rclass
	local model_type = "`1'"
	local K = "`2'"
	local Kd = "`3'"	
	local dep_var = "`4'"
	local treat_var = "`5'"
	local controls = "`6'"'	
	preserve
	capture rename `treat_var' X
	capture rename `dep_var' Y
	gen cens_ind = (X==0)
	sum X
	local samp = r(N)
	if ("`model_type'" == "naiveN") {
		sum X
		local samp = r(N)
		reg Y X
		return scalar ESAMP = e(N)
		return scalar SAMP = `samp'		
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA = _b[X]
		return scalar SE = _se[X]	
		return scalar F = .
		return scalar df = .
		return scalar df_r = .			
	}	
	if ("`model_type'" == "naive") {
		sum X
		local samp = r(N)
		reg Y X i.Z_`K' `controls', noconstant
		return scalar ESAMP = e(N)
		return scalar SAMP = `samp'		
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA = _b[X]
		return scalar SE = _se[X]		
		return scalar F = .
		return scalar df = .
		return scalar df_r = .			
	}	

	if ("`model_type'" == "tobit") {
		tobit X i.Z_`K', ll(0) noconstant
		matrix coeff_mat = e(b)
		matrix var_mat = e(V)
		predict imr_term, e(.,0)
		local sigma = coeff_mat[1,`K'+2]^0.5
		
		gen reg_term = imr_term*cens_ind + X
		levelsof Z_`K', local(levels)
		
		foreach i of local levels {
			sum imr_term if Z_`K' == `i'
			local imr_term_`i' = r(mean)
		}
	}
	if ("`model_type'" == "het_tobit") {
		*Run Tobit
		gen imr_term = .
		levelsof Z_`K', local(levels)
		foreach i of local levels {
			capture tobit X if Z_`K' == `i', ll(0)
			if _rc==0 {
				matrix coeff_mat_`i' = e(b)
				matrix var_mat_`i' = e(V)
				predict imr_temp, e(.,0)
				replace imr_term = imr_temp if Z_`K' == `i'
				drop imr_temp
				local sigma_`i' = coeff_mat_`i'[1,2]^0.5
			}
			else {
				drop if Z_`K'==`i'
			}
		}
		*get an imr_term (cens_exp) for each keep
		levelsof Z_`K', local(levels)
		foreach i of local levels {
			sum imr_term if Z_`K' == `i'
			local imr_term_`i' = r(mean)
		}
		gen reg_term = imr_term*cens_ind + X
	}

	if ("`model_type'" == "symmetric") {
		*generate a variable that tells whether cluster is censored more than 50%
		levelsof Z_`K', local(levels)
		foreach i of local levels {
		local switch_`i' = 0
			sum cens_ind if Z_`K' == `i', d
			if (r(p50) == 1) {
				local switch_`i' = 1
			}
		}
		*Get the Censored Expectations
		*Step 1 : get the censorted expectations assuming symmetry.
		*These will not be correct within a cluster if the censoring is more than 50%
		gen cens_exp = .
		foreach i of local levels {
			*Step 1 -- find F_{X|Z=z)(0) -- this is just the proportion of observations with Z=z who have X_i = 0
			sum cens_ind if Z_`K' == `i'
			local hatF_0_`i' = r(mean)
			local op_hatF_0_`i' = 1- `hatF_0_`i''

			*Step 2: find the 1-hatF_0_i percentile among X such that Z = i
			cumul X if Z_`K' == `i', gen(Xcumul_`i') equal
			sum Xcumul_`i' if Z_`K' == `i'
			local Xcumul_min = r(min)
			if `op_hatF_0_`i'' > `Xcumul_min' { 
				sum X if Z_`K' == `i' & Xcumul_`i' <= `op_hatF_0_`i''
				local step2_`i' = r(max)
			}
			if `op_hatF_0_`i'' <= `Xcumul_min' {
				local step2_`i' = 0
			}
			*Step 3: find the average value of X conditional on Z >= step2 and Z = i
			sum X if X >= `step2_`i'' & Z_`K' == `i'
			local step3_`i' = r(mean)

			*Step 4: Calculate the key censored expectation, conditional on Z = i
			local cens_expec_`i' = `step2_`i''- `step3_`i''
			replace cens_exp = `cens_expec_`i'' if Z_`K' == `i'

		}
		*also want switch codes for cases where the below qreg does not work.
		foreach i of local levels {
			capture qreg X if Z_`K' == `i', q(`op_hatF_0_`i'')
			if (_rc != 0) {
				local switch_`i' = 1
				local qreg_switch_`i' = 1
			}
		}
		*get total number of switches
		local tot_switch = 0
		foreach i of local levels {
			local tot_switch = `tot_switch' + `switch_`i''
		}
		return scalar CLUS_SWITCH = `tot_switch'
		if (`tot_switch'>0) {
			*Run Het Tobit
			gen imr_term = .
			foreach i of local levels {
				if (`switch_`i'' == 1) {
					capture tobit X if Z_`K' == `i', ll(0)
					if _rc==0 {
						matrix coeff_mat_`i' = e(b)
						matrix var_mat_`i' = e(V)
						predict imr_temp, e(.,0)
						replace imr_term = imr_temp if Z_`K' == `i'
						drop imr_temp
						local sigma_`i' = coeff_mat_`i'[1,2]^0.5
					}
					else {
						drop if Z_`K'==`i'
					}
					*get an imr_term (cens_exp) for each keep
					sum imr_term if Z_`K' == `i'
					local imr_term_`i' = r(mean)
				}
			}
			levelsof Z_`K', local(levels) // do it again because some clusters have few observations and het_tobit does not work.
		}
		gen reg_term = cens_exp*cens_ind + X
	}
	*Corrected Models
	if ("`model_type'" == "tobit") | ("`model_type'" == "het_tobit") | ("`model_type'" == "symmetric") {
		if `Kd'>1 {			
			noisily reg Y X `controls' i.Z_`K' i.Z_`Kd'#c.reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar SAMP = `samp'		
			return scalar BETA = _b[X]
			if `K'==1 {
				return scalar DELTA = _b[1.Z_`Kd'#c.reg_term]
			}
			else {
				return scalar DELTA = .	
			}
			local klist=""
			levelsof Z_`Kd', local(levelsZ)
			foreach kk of local levelsZ {
				capture return scalar DELTA_`kk'=_b[`kk'.Z_`Kd'#c.reg_term]
				if _rc~=0 {
					return scalar DELTA_`kk'=.
				}
			}
			* Implement bootstrapped F-test of whether \delta(X)=0 for all clusters of X
			if "$boot"=="1" { 
				local dlist=""	
				foreach kk of local levelsZ {
					capture display _b[`kk'.Z_`Kd'#c.reg_term]
					if _rc~=0 {
						local dlist = "`dlist'"
					}
					else {
						local dlist = "`dlist' (_b[`kk'.Z_`Kd'#c.reg_term]=${d`kk'_`model_type'})"
					}
				}
				test `dlist'			
			}
			else {
				testparm i.Z_`Kd'*#c.reg_term
			}
			return scalar F = `r(F)'
			return scalar df = `r(df)'
			return scalar df_r = `r(df_r)'	
		}
		else {
			reg Y X `controls' i.Z_`K' reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar SAMP = `samp'	
			return scalar BETA = _b[X]
			return scalar DELTA = _b[reg_term]	
			return scalar F = .
			return scalar df = .
			return scalar df_r = .
		}
	}
	noisily display "`model_type' done"
	restore
 end
}
