*! version 1.1 10February2022 
*! Stefan Tübbicke, Institute for Employment Research
*! email: stefan.tuebbicke@iab.de

program define ebct_est, eclass

	version 14.0
	syntax varlist, ///
		treatvar(varname) ///
		[sampleind(varname)] ///
		[baseweight(varname)]  ///
		[sampleweight(varname)] ///
		[out(varname)] ///
		[est(string)] ///
		[bandw(real 0.0)] ///
		[cv] ///
		[p(integer 1)] ///
		[myseed(integer 2222)] ///
		[weightname(string)] ///
		[at(string)]
	
	
	di as text ""
	di as text "Estimating balancing weights..."
	di as text ""


	mata: m_ebct`p'("`treatvar'", "`varlist'", "`sampleind'", "`baseweight'", "`sampleweight'", "`weightname'")
	
	
	tempvar gridvar
	
	qui gen `gridvar'=.
	
	tempname gridat
	mat `gridat'=`at'
	
	forvalues i=1/50{
	qui replace `gridvar'=`gridat'[`i',1] if _n==`i'
	}
		
	*Estimate DRF/derivative

	if "`out'"!=""{
	
		qui ebct_lpoly `out' `treatvar' [aw=`weightname'], degree(1) adoonly bandw(`bandw') `cv' my2seed(`myseed') at2(`gridvar')

		if "`est'"=="derivative"{
			local estimand derivative
		}

		if "`est'"=="" | "`est'"=="drf"{
			local estimand drf
		}


	tempvar sampleindhelp
	gen  `sampleindhelp'=`sampleind'

	tempname BB CC
	mat `BB'=r(`estimand')'
	ereturn post `BB', esample(`sampleindhelp')

	mat `CC'=r(gridpoints)
	ereturn matrix gridpoints=`CC'

	ereturn local estimand="`estimand'"
	ereturn local bandwidth=r(width)
	
	}

end




* slightly altered version of locpoly ado package [version 2.1.3  20nov2006 (SJ3-4: st0053; SJ6-4: st0053_3)]
program ebct_lpoly, rclass sortpreserve
	version 14.0

	syntax varlist(min=2 max=2 numeric)	///
		[if] [in] [aw fw] [,			///
		Degree(integer 0)		///
		ADOonly				///
		bandw(real 0.0) ///
		cv ///
		my2seed(integer 2222) ///
		at2(varname) ///
	]

	local kernel epanechnikov
	marksample touse
	
	* handle weights
	tempvar wgt
	gen double `wgt' = 1
	if `"`weight'"' != "" {
		qui replace `wgt' `exp' if `touse'
		if `"`weight'"' == "aweight" {
			summ `wgt' if `touse', meanonly
			qui replace `wgt' = `wgt'*r(N)/r(sum)
		}
	}
	
	
	tokenize `varlist'
	local y `1'
	local x `2'

	/*Employ user-specified bandwidth,
	obtain bandwidth via IMSE optimal bandwidth or cross-validation*/
	
	/*If bandwidth is specified, use it for subsequent regressions*/
	if `bandw'!=0{
		local bw=`bandw'
	}

	/*Obtain MISE optimal bandwidth akin to lpoly if nothing is specified*/
	if `bandw'==0 & "`cv'"==""{

		tempvar low high adj_range

		qui sum `x'
		local range=r(max)-r(min)

		qui gen `low'=r(min)+0.05*`range'
		qui gen `high'=r(max)-0.05*`range'

		qui gen `adj_range'=(inrange(`x', `low', `high'))

		qui sum `adj_range'
		local share=r(mean)
		local num=r(N)

		qui reg `y' c.`x'##c.`x'##c.`x'##c.`x' [aw=`wgt']

		tempvar secderiv
		qui gen `secderiv'=2*_b[c.`x'#c.`x']+6*_b[c.`x'#c.`x'#c.`x']*`x'+12*_b[c.`x'#c.`x'#c.`x'#c.`x']*`x'^2 

		local sigma=e(rmse)


		tempvar secderiv_sq
		gen `secderiv_sq'=`adj_range'*`secderiv'^2

		sum `secderiv_sq'
		local msdsq=r(mean)


		local bw= ((15*`share'*`sigma'^2)/(`num'*`msdsq'))^(1/5)
	
	}
	
	/*use two-fold cross-validation to obtain optimal bandwidth if option "cv" is chosen*/
	if `bandw'==0 & "`cv'"=="cv" {
	
		preserve
		set seed `my2seed'
		tempvar 2sample
		gen `2sample'=(runiform()>0.5)

		cap drop _id*
		teffects nnmatch (`y' `x') (`2sample'), gen(_id) nneighbor(1)

		qui sum `x' 
		local nobs=r(N)
		local halfspread=(r(max)-r(min))/2

		tempvar _mat_y _mat_t _mat_w
		gen `_mat_y'=.
		gen `_mat_t'=.
		gen `_mat_w'=.

	
		forvalue i=1/`nobs'{

			sum _id1 if _n==`i', meanonly
			local nnid=r(mean)

			sum `y' if _n==`nnid', meanonly
			replace `_mat_y'=r(mean) if _n==`i'

			sum `x' if _n==`nnid', meanonly
			replace `_mat_t'=r(mean) if _n==`i'

			sum `wgt' if _n==`nnid', meanonly
			replace `_mat_w'=r(mean) if _n==`i'
			
		}
		
		drop if `2sample'==0

		tempname cvres
		mat `cvres'=[.,.]

		forvalue i=20/80{

			local width=(`halfspread'/100)*`i'

			tempvar cvhyat
			lpoly `_mat_y' `_mat_t' [aw=`_mat_w'], degree(1) at(`_mat_w') ///
				gen(`cvyhat') bw(`width') nograph

			tempvar sqerr
			gen `sqerr'=(`cvyhat'-`y')^2
			sum `sqerr' [aw=`wgt']

			mat `cvres'=[`cvres'\ `width', r(mean)]

		}

		svmat `cvres'
	
		sum `cvres'2
		sum `cvres'1 if `cvres'2==r(min), meanonly

		local bw=r(mean)

		restore

	}
	


	*obtain estimates of the DRF or its derivative
	tempvar xgrid yhat dhat
	qui gen double `xgrid' = .
	qui gen double `yhat' = .	
	qui gen double `dhat' = .

	*define smoothing grid
	local n = 50
	
	qui replace `xgrid' = `at2' 
	
	ebct_Lpwork2 `y' `x' if `touse', xgrid(`xgrid') yhat(`yhat') dhat(`dhat') ///
		n(`n') h(`bw') p(`degree') k(`kernel') wgt(`wgt')	
		
	
	tempvar gridpoints _drf _deriv

	qui gen `_deriv'=`dhat'
	qui mkmat `_deriv'
	mat `_deriv'=`_deriv'[1..`n',1]
	
	qui gen `_drf'=`yhat'
	qui mkmat `_drf'
	mat `_drf'=`_drf'[1..`n',1]
	
	gen `gridpoints'=`xgrid'
	qui mkmat `gridpoints'
	mat `gridpoints'=`gridpoints'[1..`n',1]

	return matrix drf=`_drf'
	return matrix derivative=`_deriv'
	return matrix gridpoints=`gridpoints'

	qui count if `yhat' < . 
	local ngrid = r(N)
	
	return local kernel `"`kernel'"'
	return scalar width = `bw'
	return scalar ngrid = `ngrid'
	return scalar degree = `degree'

end


program ebct_Lpwork2, rclass
	syntax varlist(min=2 max=2 numeric)	///
		[if],  				///
		xgrid(varname) 			///
		yhat(varname)			///
		dhat(varname)        ///
		n(integer)			///
		wgt(varname)	///
		[ p(integer 0)			///
		h(real 0.0)			///
		k(string) ]

	tokenize `varlist'	
	local y `1'
	local x `2'

	marksample touse
	
	tempvar arg karg 
	forvalues j = 1/`p' {
		tempvar x`j'
		qui gen double `x`j'' = .
		local xs `xs' `x`j''
	}
	qui gen double `arg' = .
	qui gen double `karg' = .
	
	
	forvalues i = 1/`n' {
		qui replace `arg' = (`x' - `xgrid'[`i'])/`h' if `touse'
		ebct_GetK2 `arg' `karg' `k' `touse' `wgt'
		
		forvalues j = 1/`p' {
			qui replace `x`j'' = (`h'*`arg')^`j' if `touse'	
		}	
		qui{
		cap regress `y' `xs' [iw = `karg'] if `touse', noheader
		if !_rc {
			if _b[_cons] < . {
				qui replace `yhat' = _b[_cons] in `i'
				qui replace `dhat' = _b[`xs'] in `i'
			}
		}
		}
		
		if _rc==2000 { 
			di as err "There are gaps in the support of the treatment variable. Either increase the bandwidth or restrict the sample to a region of thicker support."
			exit 2000
		}
	
	}		
	
	
	end


program ebct_GetK2
	args arg karg kern touse wgt

	qui replace `karg' = .
	
	// epanechnikov
	local con1 = 3/(4*sqrt(5))
	local con2 = sqrt(5)
	qui replace `karg' = `wgt'* `con1'*(1-((`arg')^2/5)) /*
			*/ if abs(round(`arg',1e-8)) <= `con2' & `touse'
end


/*define mata functions*/

mata:

	/*define evaluator function*/
	void ebcteval(real scalar todo, real vector b, real scalar N, real matrix GT, real vector Q1, val, grad, hess)
		{   
			real vector  W, eGb
			real scalar sumeGb
	
			/*compute weights for given parameter vector b*/
			eGb=Q1:*exp(GT*b')
			sumeGb=sum(eGb)
	
			W=eGb/sumeGb

			/*compute objective function for given parameter vector b*/	
			val=-log(sumeGb)
	
			/*compute gradient for given parameter vector b*/
			grad=-cross(W,GT)
		}


/*Estimate EBCT weights with p=1*/		
	void m_ebct1(string scalar treatvar, string scalar varlist, string scalar sampleind, string scalar baseweight, string scalar sampleweight, string scalar ebctweight)
	{

		real vector  T, Q1, Q2, T0, b_opt, W 
		real scalar N, S
		real matrix X, X0, GT
	
		/*Obtain matrices from Stata*/
		T=.
		st_view(T, ., tokens(treatvar), sampleind)

		X=.
		st_view(X, ., tokens(varlist), sampleind)
	
		Q1=.
		st_view(Q1, ., tokens(baseweight), sampleind)
	
		Q1=Q1/sum(Q1)
	
		Q2=.
		st_view(Q2, ., tokens(sampleweight), sampleind)
	
		Q2=Q2/sum(Q2)
	
		/*Save number of obs*/
		N=rows(T)
	
		/*Normalize T and X to mean zero*/
		T0=T-J(N,1,N*mean(Q2:*T))
	
		X0=X-cross(J(cols(X),N,1),diag(N*mean(Q2:*X)))

		/*Save (transposed) G Matrix*/
		GT=(T0,X0,diag(T0)*X0)

		/*Initialize optimization*/
		S  = optimize_init()

		/*Feed in scalars/matrices*/
		optimize_init_argument(S, 1, N)
		optimize_init_argument(S, 2, GT)
		optimize_init_argument(S, 3, Q1)

		/*Set ebcteval as relevant evaluator function*/
		optimize_init_evaluator(S, &ebcteval())

		/*Choose type GF1*/
		optimize_init_evaluatortype(S, "gf1")
	
		/*Choose optimization technique*/
		optimize_init_technique(S,"nr")
	
		/*Choose optimization technique*/
		optimize_init_singularHmethod(S, "hybrid")

		/*Start at b=0*/
		optimize_init_params(S, J(1, cols(GT), 0))

		/*Run optimization and save parameter vector as b_opt*/
		b_opt = optimize(S)

		/*Generate weight vector with optimal parameters b_opt*/
		W=Q1:*exp(GT*b_opt')/sum(Q1:*exp(GT*b_opt'))

		/*replace temporary weight variable in Stata*/
		st_store(., tokens(ebctweight), sampleind, W)


	}


	
/*Estimate EBCT weights with p=2*/		
void m_ebct2(string scalar treatvar, string scalar varlist, string scalar sampleind, string scalar baseweight, string scalar sampleweight, string scalar ebctweight)
	{

		real vector  T, Q1, Q2, T0, T0_2, b_opt, W
		real scalar N, S
		real matrix X, X0, GT
	
		/*Obtain matrices from Stata*/
		T=.
		st_view(T, ., tokens(treatvar), sampleind)

		X=.
		st_view(X, ., tokens(varlist), sampleind)
	
		Q1=.
		st_view(Q1, ., tokens(baseweight), sampleind)
	
		Q1=Q1/sum(Q1)
	
		Q2=.
		st_view(Q2, ., tokens(sampleweight), sampleind)
	
		Q2=Q2/sum(Q2)
	
		/*Save number of obs*/
		N=rows(T)
	
		/*Normalize T and X to mean zero*/
		T0=T-J(N,1,N*mean(Q2:*T))
	
		X0=X-cross(J(cols(X),N,1),diag(N*mean(Q2:*X)))
	
		/*Vector of squares in T*/
		T0_2=diag(T)*T-J(N,1,N*mean(Q2:*diag(T)*T))

		/*Save (transposed) G Matrix*/
		GT=(T0,T0_2,X0,diag(T0)*X0,diag(T0_2)*X0)

		/*Initialize optimization*/
		S  = optimize_init()

		/*Feed in scalars/matrices*/
		optimize_init_argument(S, 1, N)
		optimize_init_argument(S, 2, GT)
		optimize_init_argument(S, 3, Q1)

		/*Set ebcteval as relevant evaluator function*/
		optimize_init_evaluator(S, &ebcteval())

		/*Choose type GF1*/
		optimize_init_evaluatortype(S, "gf1")
	
		/*Choose optimization technique*/
		optimize_init_technique(S,"nr")
	
		/*Choose optimization technique*/
		optimize_init_singularHmethod(S, "hybrid")

		/*Start at b=0*/
		optimize_init_params(S, J(1, cols(GT), 0))

		/*Run optimization and save parameter vector as b_opt*/
		b_opt = optimize(S)

		/*Generate weight vector with optimal parameters b_opt*/
		W=Q1:*exp(GT*b_opt')/sum(Q1:*exp(GT*b_opt'))

		/*replace temporary weight variable in Stata*/
		st_store(., tokens(ebctweight), sampleind, W)
		
	}


/*Estimate EBCT weights with p=3*/		
void m_ebct3(string scalar treatvar, string scalar varlist, string scalar sampleind, string scalar baseweight, string scalar sampleweight, string scalar ebctweight)
	{

		real vector  T, Q1, Q2, T0, T0_2, T0_3, b_opt, W
		real scalar N, S
		real matrix X, X0, GT
	
		/*Obtain matrices from Stata*/
		T=.
		st_view(T, ., tokens(treatvar), sampleind)

		X=.
		st_view(X, ., tokens(varlist), sampleind)
	
		Q1=.
		st_view(Q1, ., tokens(baseweight), sampleind)
	
		Q1=Q1/sum(Q1)
	
		Q2=.
		st_view(Q2, ., tokens(sampleweight), sampleind)
	
		Q2=Q2/sum(Q2)
	
		/*Save number of obs*/
		N=rows(T)
	
		/*Normalize T and X to mean zero*/
		T0=T-J(N,1,N*mean(Q2:*T))
	
		X0=X-cross(J(cols(X),N,1),diag(N*mean(Q2:*X)))
	
		/*Vector of demeaned squares in T*/
		T0_2=diag(T)*T-J(N,1,N*mean(Q2:*diag(T)*T))

		/*Vector of demeaned qubic terms in T*/
		T0_3=diag(T)*diag(T)*T-J(N,1,N*mean(Q2:*diag(T)*diag(T)*T))

		/*Save (transposed) G Matrix*/
		GT=(T0,T0_2,T0_3,X0,diag(T0)*X0,diag(T0_2)*X0,diag(T0_3)*X0)

		/*Initialize optimization*/
		S  = optimize_init()

		/*Feed in scalars/matrices*/
		optimize_init_argument(S, 1, N)
		optimize_init_argument(S, 2, GT)
		optimize_init_argument(S, 3, Q1)

		/*Set ebcteval as relevant evaluator function*/
		optimize_init_evaluator(S, &ebcteval())

		/*Choose type GF1*/
		optimize_init_evaluatortype(S, "gf1")
	
		/*Choose optimization technique*/
		optimize_init_technique(S,"nr")
	
		/*Choose optimization technique*/
		optimize_init_singularHmethod(S, "hybrid")

		/*Start at b=0*/
		optimize_init_params(S, J(1, cols(GT), 0))

		/*Run optimization and save parameter vector as b_opt*/
		b_opt = optimize(S)

		/*Generate weight vector with optimal parameters b_opt*/
		W=Q1:*exp(GT*b_opt')/sum(Q1:*exp(GT*b_opt'))
				
		/*replace temporary weight variable in Stata*/
		st_store(., tokens(ebctweight), sampleind, W)
		
	}


end
