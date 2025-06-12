*! version 1.1 10February2022 
*! Stefan Tübbicke, Institute for Employment Research
*! email: stefan.tuebbicke@iab.de


program define ebct, eclass 

	version 14.0
	syntax varlist [if] [in] , ///
		treatvar(varname) ///
		[basew(varname)]  ///
		[samplew(varname)] ///
		[p(integer 1)] ///
		[out(varname)] ///
		[est(string)] ///
		[cv] ///
		[bandw(real 0.0)] ///
		[bootstrap] ///
		[reps(real 100)] ///
		[graph] ///
		[generate(string)] ///
		[seed(integer 2222)]
	 

	/*Organize sample & base weights*/
	tempvar r2_u F_u p_u sampleind baseweight sampleweight

	if "`samplew'"==""{
		qui gen `sampleweight'=1 `if' `in'
	}

	else qui gen `sampleweight'=`samplew'  `if' `in'

	if "`basew'"==""{
		qui gen `baseweight'=1 `if' `in'
	}

	else qui gen `baseweight'=`basew' `if' `in'

	if `bandw'!=0 local display_help "user-specified"
	if `bandw'==0 & "`cv'"=="" local display_help "Fan/Gijbels (1992) plug-in"
	if `bandw'==0 & "`cv'"=="cv" local display_help "cross-validated"


	/*Obtain initial measures of covariate balance*/
	qui reg `treatvar' `varlist' `if' `in' [aw=`sampleweight']
	scalar `r2_u'=e(r2)
	scalar `F_u'=e(F)
	scalar `p_u'=1-F(e(df_m),e(df_r),e(F))

	*Mark estimation sample
	qui gen `sampleind'=(e(sample)==1 & `baseweight'!=.)

	/*remove colinear variables from varlist*/
	_rmcoll `varlist' if `sampleind', forcedrop
	local varlist "`r(varlist)'"

	*set degree of polynomial in T to be rendered uncorrelated with X
	if "`p'"=="" local p=1
	if `p'<1 local p=1
	if `p'>3 local p=3



	di as text "Note: Balancing weights are estimated such that T,...,T^P is uncorrelated with X, P=`p'"


	***Estimate balancing weights using EBCT, and 
	***the DRF/its derivatives using local polynomial regression
	
	
	tempvar weight
	qui gen `weight'=.
	
	tempvar grid
	qui gen `grid'=.
	
	qui sum `treatvar' `if' `in'

	qui replace `grid' = r(min)+(_n-1)*(r(max)-r(min))/49 in 1/50
	
	tempname gridmat
	mkmat `grid', matrix(`gridmat')
	
	mat `gridmat'=`gridmat'[1..50,1]
		
	
	ebct_est `varlist', treatvar(`treatvar') sampleind(`sampleind') ///
		baseweight(`baseweight') sampleweight(`sampleweight') ///
		out(`out') est(`est') bandw(`bandw') `cv' p(`p') myseed(`seed') ///
		weightname(`weight') at(`gridmat')
	
	
	if "`out'"!=""{

		*Save DRF/derivative, gridoints and the kernel bandwidth
		tempname AA AAA
		mat `AA'=e(b)
		tempvar ehelp
		gen `ehelp'=e(sample)
		mat `AAA'=e(gridpoints)'

		local estimand=e(estimand)
		local bandwidth=e(bandwidth)
	}

	/*Obtain post-weighting balance statistics and print table*/
	tempvar r2_w F_w p_w obs
	qui reg `treatvar' `varlist' [aw=`weight']
	scalar `r2_w'=e(r2)
	scalar `F_w'=e(F)
	scalar `p_w'=1-F(e(df_m),e(df_r),e(F))
	scalar `obs'=e(N)

	mat balance=(`r2_u',`F_u', `p_u' \ `r2_w',`F_w', `p_w')'
	mat rownames balance = R-squared F-statistic p-value
	mat colnames balance = "before balancing" "after balancing"

	/*print table*/
	di as text ""
	di as text "###########################################"
	di as text "# Summary statistics on balancing quality #"
	di as text "###########################################"

	matlist balance, rowtitle("") border(rows)  twidth(16) format(%16.3f) ///
		nodotz title("Summary statistics from a (weighted) regression of T on X:")


	/*Generate flexible balance indicators*/
	tempname r2_bef r2_aft p_bef p_aft

	local firstvar : word 1 of `varlist'

	
	foreach var of varlist `varlist'{

		qui reg `var' c.`treatvar'##c.`treatvar'##c.`treatvar'##c.`treatvar'
		scalar `r2_bef'=e(r2)
		scalar `p_bef'=1-F(e(df_m), e(df_r), e(F))

		qui reg `var' c.`treatvar'##c.`treatvar'##c.`treatvar'##c.`treatvar' [aw=`weight']
		scalar `r2_aft'=e(r2)
		scalar `p_aft'=1-F(e(df_m), e(df_r), e(F))


		if "`var'"=="`firstvar'" mat r2_balance=[`r2_bef', `p_bef', `r2_aft', `p_aft']
		else mat r2_balance=[r2_balance\ `r2_bef', `p_bef', `r2_aft', `p_aft']
	}

	#del;
	mat colnames r2_balance = 
	"Before Balancing:R-square" 
	"Before Balancing:p-value" 
	"After Balancing:R-square" 
	"After Balancing:p-value"
	;
	#del cr

	mat rownames r2_balance = `varlist'

	*Print table
	di as text ""
	di as text "###########################################"
	di as text "#     Details on balance and weights      #"
	di as text "###########################################"

	matlist r2_balance, rowtitle("") border(rows)  twidth(16) format(%9.3f) ///
		nodotz showcoleq(combined) ///
		title("Balancing statistics from a regression of each element of X on a fourth-order polynomial in T:")

	di as text ""	
	
	qui sum `weight', d

	tempname dist
	mat `dist'=[r(min)\ r(p1)\ r(p5)\ r(p10)\ r(p25)\ r(p50)\ r(p75)\ r(p90)\ r(p95)\ r(p99)\ r(max)]
	
	#del;
	mat rownames `dist' = 
	"min" 
	"1st percentile" 
	"5th percentile"
	"10th percentile"
	"25th percentile"
	"Median"
	"75th percentile"
	"90th percentile"
	"95th percentile"
	"99th percentile"
	"max"
	;
	#del cr	
	
	mat colnames `dist' = "value"
	
	matlist `dist', rowtitle("") border(rows)  twidth(16) format(%9.3f) ///
		nodotz showcoleq(combined) ///
		title("Distribution of estimated weights:")
	
	
	di as text ""
	if "`estimand'"=="drf" di as text "Estimating the DRF using local linear regressions with `display_help' bandwidth."
	if "`estimand'"=="derivative" di as text "Estimating the derivatives of the DRF using local linear regressions with `display_help' bandwidth."


	if "`out'"!="" & "`bootstrap'"==""{
		di as text "For the estimation of standard errors and confidence intervals, specify the bootstrap option."
	}


	*Estimation of standard errors using the bootstrap routine
	if "`out'"!="" & "`bootstrap'"!=""{

		di as text "Obtaining Standard Errors using the Bootstrap."
		di as text "This may take a while..."
		di as text ""

		bootstrap _b, seed(`seed') reps(`reps') nowarn notable noheader ///
			nolegend: ebct_est `varlist', treatvar(`treatvar') ///
			sampleind(`sampleind') baseweight(`baseweight') ///
			sampleweight(`sampleweight') out(`out') p(`p') est(`est') ///
			myseed(`seed') weightname(`weight') at(`gridmat') bandw(`bandwidth')
			
		tempname DD
		mat `DD'=e(se)
	}

	/*Post results in e()*/
	if "`out'"==""{
		ereturn post 
		ereturn matrix balance=balance
		ereturn matrix balance_detail=r2_balance
		ereturn local cmdname="ebct"
		ereturn scalar N=`obs'
	}

	if "`out'"!=""{
		ereturn post `AA' , esample(`ehelp')
		ereturn matrix gridpoints=`AAA'
		ereturn local estimand="`estimand'"
		ereturn scalar bandwidth=`bandwidth'
		ereturn local estimator="Local linear regression"
		ereturn local kernel="Epanechnikov"
		ereturn matrix balance=balance
		ereturn matrix balance_detail=r2_balance
		ereturn local cmdname="ebct"
		ereturn scalar N=`obs'


	}

	if "`out'"!="" & "`bootstrap'"!=""{
		ereturn matrix se=`DD'

	}


	/*generate temporary variables if needed later on*/

	if "`graph'"!="" | "`generate'"!=""{

		tempname yhat xgrid 
		
		mat `yhat'=e(b)'
		mat `xgrid'=e(gridpoints)'

		tempvar estimates gridpoints
		qui gen `estimates'=.
		qui gen `gridpoints'=.
		
		forvalue i=1/50{
			qui replace `estimates'=`yhat'[`i',1] if _n==`i'
			qui replace `gridpoints'=`xgrid'[`i',1] if _n==`i'
		}

		
		if "`bootstrap'"!=""{


			tempname sehat 
			mat `sehat'=e(se)'

			tempvar stderr 
			qui gen `stderr'=.
		
			forvalue i=1/50{
				qui replace `stderr'=`sehat'[`i',1] if _n==`i'
			}
		
			tempvar ci_low ci_high
			qui gen `ci_low'=`estimates'-1.96*`stderr'
			qui gen `ci_high'=`estimates'+1.96*`stderr'


		}


	}


	/*Produce graph if requested*/

	if "`out'"!="" & "`graph'"!=""{


		local xttl : var label `treatvar'
	
		if `"`xttl'"' == "" {
			local xttl `treatvar'
		}


		local yttl: var label `out'

		if `"`yttl'"' == "" {
			local yttl `out'
		}

		
		
		*Do not draw confidence intervals without bootstrap
		if "`bootstrap'"==""{


			if "`estimand'"=="drf"{
				#del;
				twoway (line `estimates' `gridpoints', 
				lcolor(black) lpattern(solid)) , 
				ytitle("`yttl'") xtitle("t = `xttl'")
				legend(on) 
				legend(label(1 "E[Y(t)]"))
				scheme(s2mono) graphregion(color(white)) bgcolor(white)
				;
				#del cr
			}
			

		
			if "`estimand'"=="derivative"{
				#del;
				twoway (line `estimates' `gridpoints',
				lcolor(black) lpattern(solid)),
				ytitle("`yttl'") xtitle("t = `xttl'") 
				legend(on) 
				legend(label(1 "dE[Y(t)]/dt"))
				scheme(s2mono) graphregion(color(white)) bgcolor(white)
				;
				#del cr
			}
		}																			

		*Also draw confidence intervals with bootstrap
		if "`bootstrap'"!=""{
		
		
			if "`estimand'"=="drf"{
				#del;
				twoway (line `estimates' `gridpoints', lcolor(black) lpattern(solid))  
				(rarea `ci_low' `ci_high' `gridpoints', color(gs8%30)), 
				ytitle("`yttl'") xtitle("t = `xttl'") ///
				legend(order(1 "E[Y(t)]" 2 "95% Confidence Interval")) 
				scheme(s2mono) graphregion(color(white)) bgcolor(white)
				;
				#del cr
			}

			if "`estimand'"=="derivative"{
				#del;
				twoway (line `estimates' `gridpoints', lcolor(black) lpattern(solid))  
				(rarea `ci_low' `ci_high' `gridpoints', color(gs8%30)), 
				ytitle("`yttl'") xtitle("t = `xttl'") 
				legend(order(1 "dE[Y(t)]/dt" 2 "95% Confidence Interval")) 
				scheme(s2mono) graphregion(color(white)) bgcolor(white)
				;
				#del cr
			}
	
		}

	}
	
	
	
	*Produce variables if requested
	if "`generate'"!=""{
	
		di as text ""
		di as text "generating variables `generate'_*"
	
		cap drop `generate'_weight
		qui gen `generate'_weight=`weight'
		label var `generate'_weight "Estimated EBCT weights"
	
		if "`generate'"!="" & "`out'"!=""{

			cap drop `generate'_estimates
			qui gen `generate'_estimates=`estimates'

			if "`est'"=="drf" | "`est'"=="" label var `generate'_estimates "Estimated DRF"
			if "`est'"=="derivative" label var `generate'_estimates "Estimated derivative of the DRF"

			cap drop `generate'_grid
			qui gen `generate'_grid=`gridpoints'
			label var `generate'_grid "Gridpoints"
		}
		
		if "`generate'"!="" & "`out'"!="" & "`bootstrap'"!=""{

			cap drop `generate'_stderr
			cap drop `generate'_ci_low
			cap drop `generate'_ci_up
					
			qui gen `generate'_stderr=`stderr'
			qui gen `generate'_ci_low=`ci_low'
			qui gen `generate'_ci_up=`ci_high'
					
				if "`est'"=="drf" {
					label var `generate'_stderr "Std. err. of the estimated DRF"
					label var `generate'_ci_low "Lower bound of 95% CI for the DRF"
					label var `generate'_ci_up "Upper bound of 95% CI for the DRF"

				}
					
				if "`est'"=="derivative" {
					label var `generate'_stderr "Std. err. of the estimated derivative of the DRF"
					label var `generate'_ci_low "Lower bound of 95% CI for the derivative of the DRF"
					label var `generate'_ci_up "Upper bound of 95% CI for the derivative of the DRF"

				}				
					
					
		}
	
	

	}

end






