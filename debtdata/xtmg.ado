***************************************************************************************************
*!version 1.1.0 2May2013
*!xtmg   - 	Estimating panel time series models with heterogeneous slopes 
* 		by Markus Eberhardt, School of Economics, University of Nottingham
* 		For feedback please email me at markus.eberhardt@nottingham.ac.uk
* 		Visit http://sites.google.com/site/medevecon/ for macro panel data 
*			and other Stata routines
***************************************************************************************************
* xtmg
***************************************************************************************************
* Known Bugs: 	
*		-Need to drop groups which do not enter regression for FM/DOLS
*		-Not convinced the bwidth selection method works for DOLS (vis-a-vis dic).
*
* Revisions made:
* v.1.0.1 (7Feb2011)
*		-corrected some typos which prevented xtmg from running under -varabbrev off-
*		-changed parts of aug-routine first stage to mata to avoid large matsize
*		-nocons now works in the standard MG estimator
*		-added error messages for cce and aug options if -nocons- selected
*
* v.1.0.2 (3Jan2012)
*		-cleaned up some labels for regression output (cross-section averages, trends, etc.)
*		-minor changes to the presentation of the results
*		-added option to compute predicted values based on heterogeneous coefficients
*
* v.1.0.3 (2Aug2012)
*       -corrected a bug for CCEMG whereby groups with insufficient observations were previously
*        not dropped from the sample: each individual/group regression requires (2k+2) time-series 
*        observations.
*
* v.1.1.0 (2May2013)
*       -added group mean FMOLS and DOLS to the suite of estimators. Necessitates move to Stata 11.
*
***************************************************************************************************
* Acknowledgements:
*
* This routine builds to a considerable extent on the existing code for the Swamy RCM estimator 
* (xtrc), the Pesaran, Shin and  Smith (1999) Pooled Mean Group estimator written by Edward F. 
* Blackburne III and Mark W. Frank (xtpmg), and the Westerlund (2007) error correction
* cointegration test (xtwest) written by Damiaan Persyn and Joakim Westerlund. 
* Thanks to Kit Baum and a Stata Journal reviewer for useful comments, help and support. 
* Any remaining errors are my own.
*
* The revised 2013 version wraps user-written routines lrcov and cointreg for the single time 
* series. These commands were written by Qunyong Wang and Na Wu (see Stata Journal Vol.12(3), 2012).
*
***************************************************************************************************

capture program drop xtmg
program define xtmg, eclass prop(xt) sort
version 11

      syntax varlist(ts fv) [if] [in] [, I(varname) T(varname) TREND noCONstant AUGment IMPose CCE ///
		FULL Level(cilevel) ROBUST RES(namelist) PRED(namelist) FM DOLS ///
		dlead(integer 1) dlag(integer 1) dic(string) KERNel(string) bwidth(integer 1) ///
		bmeth(string) ]

      _xt, i(`i') t(`t')
      local ivar "`r(ivar)'"
      local tvar "`r(tvar)'"
	marksample touse `dvar'
	markout `touse' `offset' `t' 
	markout `touse' `ivar', strok
	
	if ("`dols'" != "" | "`fm'" != "" ) {
	tsreport if `touse', report panel
	if r(N_gaps) {
		di as error _col(2) ""
		di as error _col(2) "For FMOLS or DOLS the sample may not contain gaps."
		exit
		}
	}
	
	if ("`dols'" != "" | "`fm'" != "" ) {
	qui xtsum `tvar' if `touse'
	local N `r(n)'
	local T `r(Tbar)'
	if int(`T')*`N' ~= r(N) {
		di as error _col(2) ""
		di as error _col(2) "For FMOLS or DOLS the panel must be balanced."
		exit
		}
	}
	
quietly{

/* Tokenize varlist and determine dimensions */
		tokenize `varlist'
            local dep "`1'"
		local depname "`1'"
            mac shift
            local ind "`*'"
            noi _rmcoll `ind' if `touse', `constant'
            local ind "`r(varlist)'"
            local p : word count `ind'
		local rhs = `p'
		if ("`augment'" != "") & ("`cce'" != "" | "`dols'" != "" | "`fm'" != "") {
			display as error _col(2) ""
			display as error _col(2) "You can only select one of options -cce-, -fm-, -dols- or -augment-."
			exit
		} 
		if ("`cce'" != "") & ("`dols'" != "" | "`fm'" != "") {
			display as error _col(2) ""
			display as error _col(2) "You can only select one of options -cce-, -fm- or  -dols-."
			exit
		} 
		if ("`impose'" != "") & ("`augment'" == "") {
			display as error _col(2) ""
			display as error _col(2) "You cannot select option -impose- without -augment-."
			exit
		} 
		if ("`cce'"!="" & "`constant'" != "") {
			display as error ""
			display as error _col(2) "The CCEMG routine requires a group-specific intercept."
			display as error _col(2) "Please drop the -nocons- option."
			exit
		}
		if ("`augment'"!="" & "`constant'" != "") {
			display as error ""
			display as error _col(2) "The Augmented MG estimator requires a group-specific intercept."
			display as error _col(2) "Please drop the -nocons- option."
			exit
		}


		if ("`constant'" == "") { 
			local rhs = `rhs'+1 
		}
		if ("`trend'" != "") { 
			local rhs = `rhs'+1 
		}
		if ("`augment'" != "" & "`impose'"== "") { 
			local rhs = `rhs'+1 
		}
		else local rhs = `rhs'

        tempvar t T Tall
        sort `touse' `ivar' 
        by `touse' `ivar': gen int `t' = _n if `touse'
        by `touse' `ivar': gen int `T' = _N if `touse'
		count if `touse' 
		local nobso = r(N)
		
		if ("`cce'"!="") {
			by `touse' `ivar' : replace `touse' = 0 if `T'[_N] <= (2*`rhs')+1
			replace `T' = . if `touse'==0
			count if `touse' 
			local nobs = r(N)
		}
		else {
			by `touse' `ivar' : replace `touse' = 0 if `T'[_N] <= `rhs'
			replace `T' = . if `touse'==0
			count if `touse' 
			local nobs = r(N)
		}

		if `nobs' < `nobso' {
			noi di as text _col(2) "Note: " as res `nobso'-`nobs' as text " obs. dropped (panels too small)" 
		}

		if ("`kernel'"=="") local kernel="bartlett"
		if ("`bmeth'"=="") local bmeth="neweywest"
		

/* Create group variable and time trend */
		tempvar g gall timevar  mm
		egen `g' = group(`ivar') if `touse'
		summ `g' if `touse'
		local ng = r(max)
		egen `gall' = group(`ivar')
		summ `gall'
		local ngall = r(max)
			if ("`dols'" != "" | "`fm'" != "" ) & (`ngall'>`ng') {
				display as error ""
				di as error _col(2) "Due to a bug in the code the Group Mean FMOLS/DOLS estimators cannot" _n ///			
					_col(2) "deal with panel groups which do not enter the regression. Please drop" _n ///
					_col(2) "the following groups (-drop if `ivar'==...-) from the dataset. Note that " _n ///
					_col(2) "these observations would *not* enter the regression anyway."
				noi tab `ivar' if missing(`g') & !missing(`gall')
				exit
			}	
		summ `tvar' 
        local Tall = r(max)
		gen `timevar'_t=`tvar'-r(min)
		summ `timevar'_t
		local gmax = r(max)+1

            summarize `T' if `touse' & `ivar'~=`ivar'[_n-1], meanonly
            local n = r(N)
            local g1 = r(min)
            local g2 = r(mean)
            local g3 = r(max)

		if `c(matsize)' <`ng'{
			noi display in smcl  _col(2) "{help matsize##|_new:matsize} " ///
			as error "must be at least as large as the number" _n ///
			as error "of panels in the current estimation sample (`ng')."
			exit 908
		}
		


/* Compute country-estimates */
		tempname bols vari bbar vce  depT indT tmp tmp2 tmp3 linear r r1 yhat1 sig2 beta names yhat yhat1 sample
		generate `r'=.
		generate `yhat'=.
		if ("`cce'" != ""){
			if ("`trend'" != ""){
			mat `bols' = J(`ng', 2*(`rhs'-1)+1, 0)
			mat `vari' = J(`ng', 2*(`rhs'-1)+1, 0)
			mat `bbar' = J(`ng', 2*(`rhs'-1)+1, 0)
			mat `vce' = J(2*(`rhs'-1)+1, 2*(`rhs'-1)+1, 0)
			}			
			else {
			mat `bols' = J(`ng', 2*(`rhs'-1)+2, 0)
			mat `vari' = J(`ng', 2*(`rhs'-1)+2, 0)
			mat `bbar' = J(`ng', 2*(`rhs'-1)+2, 0)
			mat `vce' = J(2*(`rhs'-1)+2, 2*(`rhs'-1)+2, 0)
			}
		}
		else {
			mat `bols' = J(`ng', `rhs', 0)
			mat `vari' = J(`ng', `rhs', 0)
			mat `bbar' = J(`ng', `rhs', 0)
			mat `vce' = J(`rhs', `rhs', 0)
		}

/* Country regressions start here */
		
		local i = 1

		
	/* Augmented MG estimator: first stage */
		if ("`augment'" != ""){
			tsset `ivar' `tvar'
			tempvar yr year ddep 
			tempname comdyn zeros comdyn2  aug `aug'_c
			tab `tvar', gen(`yr')
			forvalues l=1/`gmax'{
				tsset `ivar' `tvar'
				gen `year'`l'=d.`yr'`l'
			}
			gen `ddep'=d.`dep' if `touse'
			reg `ddep' d.(`ind') `year'2-`year'`gmax', nocons
			mat `comdyn'=e(b)'
			mat comdyn=`comdyn'
			mat `zeros'=J(1,1,0)
			if ("`impose'" == ""){
				if ("`trend'" != ""){
					mat `comdyn'=`zeros' \ `comdyn'[(`rhs'-2)...,1]
				} 
				else mat `comdyn'=`zeros' \ `comdyn'[(`rhs'-1)...,1]
			}
			if ("`impose'" != ""){
				if ("`trend'" != ""){
					mat `comdyn'=`zeros' \ `comdyn'[(`rhs'-1)...,1]
				} 
				else mat `comdyn'=`zeros' \ `comdyn'[(`rhs')...,1]
			}
			mata: `comdyn' = st_matrix("`comdyn'")
			mata: `comdyn2' = J(`ngall',1,1)#`comdyn'
			gen `aug'_c=0
			mata: st_store(.,"`aug'_c",`comdyn2')
			sort `ivar' `tvar'
		}
		if ("`augment'" != "" & "`impose'" != "") {
			tempvar depaug
			gen `depaug'=`dep'-`aug'_c if `touse'
		}
		if ("`augment'" != "" & "`impose'" != ""){
				if ("`trend'" != ""){		
					reg `depaug' `ind' `timevar'_t if `touse' & `g'==`i', `constant'
				}
				if ("`trend'" == ""){	
					reg `depaug' `ind' if `touse' & `g'==`i', `constant' 
				} 
		}
		if ("`augment'" != "" & "`impose'" == ""){
				if ("`trend'" != ""){
					reg `dep' `ind' `aug'_c `timevar'_t if  `touse' & `g'==`i', `constant'
				}
				if ("`trend'" == ""){
					reg `dep' `ind' `aug'_c if `g'==`i', `constant'
				}
		}
		

	/* CCEMG estimator: group 1 */

		if ("`cce'" != ""){
			tempvar indT depT
			tokenize `varlist'
			sort `tvar' `ivar' 
			local m = 1
			tempvar `depT'_`m'
			by `tvar': egen `depT'_``m''=mean(``m'') if `touse' 
		      local coefs : word count `varlist'
			forvalues m = 2/`coefs'{
				tempvar `indT'_`m'
				by `tvar': egen `indT'_``m''=mean(``m'') if `touse'  	
			}
			sort `ivar' `tvar'
			if ("`trend'" != ""){
				reg `dep' `ind' `timevar'_t `depT'_* `indT'_* if `touse' & `g'==1, `constant'
			}
			else {
				reg `dep' `ind' `depT'_* `indT'_* if `touse' & `g'==1, `constant'
			}
	
		}

		
	/* Group Mean FM estimator: group 1 */
		if ("`fm'" != ""){
			if ("`trend'" != ""){
				cointreg `dep' `ind' if `touse' & `g'==1, est(fmols) eqtrend(1) kernel(`kernel') bmeth(`bmeth') 
				local waldtrend=0
				test linear
				if r(p)<(100-`level')/100{
					 local waldtrend=`waldtrend'+1
				}
			}
			else cointreg `dep' `ind' if `touse' & `g'==1, est(fmols) eqtrend(0) kernel(`kernel') bmeth(`bmeth') 
		}

		
	/* Group Mean DOLS estimator: group 1 */
		if ("`dols'" != "") {
			if ("`trend'" != ""){		
				cointreg `dep' `ind' if `touse' & `g'==1, est(dols) eqtrend(1) dic(`dic') kernel(`kernel')  bmeth(`bmeth') dlead(`dlead') dlag(`dlag')
				local waldtrend=0
				test linear
				if r(p)<(100-`level')/100{
					 local waldtrend=`waldtrend'+1
				}
			}
			else cointreg `dep' `ind' if `touse' & `g'==1, est(dols) eqtrend(0) dic(`dic') kernel(`kernel')  bmeth(`bmeth')  dlead(`dlead') dlag(`dlag')
		}


	/* Standard MG estimator: group 1 */

		if  ("`cce'" == "" & "`augment'" == "" & "`fm'"=="" & "`dols'"==""){
			if ("`trend'" != ""){
				reg `dep' `ind' `timevar'_t if `touse' & `g'==1, `constant'
			}
			else {
				reg `dep' `ind' if `touse' & `g'==1, `constant'
			}
		}
		
		
	/* Predict residuals for group 1 */
		if ("`fm'" != "" | "`dols'"!="") {
			mat `tmp' = get(_b)	
			local colstr : colnames `tmp'
			gen linear=`timevar'_t if `g'==1
			mat colnames `bols' = `colstr'
			mat `tmp2' = get(VCE)	
			predict double `yhat1' if `touse' & `g'==1
			drop linear
			gen `r1' = `dep'-`yhat1' if `touse' & `g'==1
			replace `r'=`r1' if `touse' & `g'==1
			replace `yhat'=`yhat1' if `touse' & `g'==1
			mat `bols'[1, 1] = `tmp'
			mat `vari'[1, 1] =  vecdiag(`tmp2')
			mat `bbar' = `bols'[1, 1...]
			}
		else {
			mat `tmp' = get(_b)	
			local colstr : colnames `tmp'
			mat colnames `bols' = `colstr'
			mat `tmp2' = get(VCE)	
			if ("`trend'" != ""){
				test `timevar'_t
				local waldtrend=0
					if r(p)<(100-`level')/100{
						 local waldtrend=`waldtrend'+1
					}
			}
			predict double `r1' if `touse' & `g'==1, res
			replace `r'=`r1' if `touse' & `g'==1
			predict double `yhat1' if `touse' & `g'==1
			replace `yhat'=`yhat1' if `touse' & `g'==1
			mat `bols'[1, 1] = `tmp'
			mat `vari'[1, 1] =  vecdiag(`tmp2')
			mat `bbar' = `bols'[1, 1...]
		}

		
	/* Groups 2 to N */

		local i = 2
		while `i' <= `ng'{
			tempvar `r'`i' `yhat'`i'

			
	/* AMG estimator: groups 2 to N */
			if ("`augment'" != "" & "`impose'" == "") {
				if ("`trend'" != ""){
					reg `dep' `ind' `aug'_c `timevar'_t if `touse' & `g'==`i', `constant'
				}
				else reg `dep' `ind' `aug'_c if  `touse' & `g'==`i'
			}
			if ("`augment'" != "" & "`impose'" != "") {
				if ("`trend'" != ""){		
					reg `depaug' `ind' `timevar'_t if `touse' & `g'==`i', `constant'
				}
				else reg `depaug' `ind' if `touse' & `g'==`i', `constant'
			}

			
	/* CCEMG estimator: groups 2 to N */
			if ("`cce'" != ""){
				if ("`trend'" != ""){		
					reg `dep' `ind' `timevar'_t `depT'_*  `indT'_* if `touse' & `g'==`i', `constant'
				}
				else reg `dep' `ind' `depT'_*  `indT'_* if `touse' & `g'==`i', `constant'
			}

			
	/* Standard MG estimator: groups 2 to N */
		if  ("`cce'" == "" & "`augment'" == "" & "`fm'"=="" & "`dols'"==""){
				if ("`trend'" != ""){		
					reg `dep' `ind' `timevar'_t if `touse' & `g'==`i', `constant'
				}
				else reg `dep' `ind' if `touse' & `g'==`i', `constant'
			}

			
	/* Group Mean FM estimator: groups 2 to N */
		if ("`fm'" != "") {
			if ("`trend'" != ""){		
					cointreg `dep' `ind' if `touse' & `g'==`i', est(fmols) eqtrend(1) kernel(`kernel') bmeth(`bmeth') 
					test linear
					if r(p)<(100-`level')/100{
						 local waldtrend=`waldtrend'+1
					}
				}
				else cointreg `dep' `ind' if `touse' & `g'==`i', est(fmols) eqtrend(0) kernel(`kernel')  bmeth(`bmeth') 
			}
		
	
	/* Group Mean DOLS estimator: groups 2 to N */
		if ("`dols'" != "") {
				if ("`trend'" != ""){		
				cointreg `dep' `ind' if `touse' & `g'==`i', est(dols) eqtrend(1) dic(`dic')  kernel(`kernel')  bmeth(`bmeth')  dvar(`ind') dlead(`dlead') dlag(`dlag')
				test linear
					if r(p)<(100-`level')/100{
						 local waldtrend=`waldtrend'+1
					}
				}
				else cointreg `dep' `ind' if `touse' & `g'==`i', est(dols) eqtrend(0) dic(`dic')  kernel(`kernel')  bmeth(`bmeth')  dvar(`ind') dlead(`dlead') dlag(`dlag')
			}
			
			
	/* Compute residuals and predicted values: groups 2 to N */
	if ("`fm'" != "" | "`dols'"!="") {
			mat `tmp' = get(_b)	
			mat `tmp2' = get(VCE)	
			gen linear=`timevar'_t if `g'==`i'
			if ("`trend'" != ""){
				qui test linear
				if r(p)<(100-`level')/100{
					 local waldtrend=`waldtrend'+1
				}
			}
			predict double `yhat'`i' if `touse' & `g'==`i'
			replace `yhat'=`yhat'`i' if `touse' & `g'==`i'
			gen `r'`i' = `dep'-`yhat'`i' if `touse' & `g'==`i'
			replace `r'=`r'`i' if `touse' & `g'==`i'
			drop linear
			mat `bols'[`i', 1] = `tmp'
			mat `vari'[`i', 1] =  vecdiag(`tmp2')
			mat `bbar' = `bbar' + `bols'[`i', 1...]
		}
		else{ 
			predict double `r'`i' if `touse' & `g'==`i', res
			mat `tmp' = get(_b)	
			mat `tmp2' = get(VCE)	
			if ("`trend'" != ""){
				qui test `timevar'_t
				if r(p)<(100-`level')/100{
					 local waldtrend=`waldtrend'+1
				}
			}
			replace `r'=`r'`i' if `touse' & `g'==`i'
			predict double `yhat'`i' if `touse' & `g'==`i'
			replace `yhat'=`yhat'`i' if `touse' & `g'==`i'
			mat `bols'[`i', 1] = `tmp'
			mat `vari'[`i', 1] =  vecdiag(`tmp2')
			mat `bbar' = `bbar' + `bols'[`i', 1...]
			}
			local i = `i' + 1
		}
		
	/* End of group 2 to N estimation and thus all country regressions */

		local ngg = 1/`ng'
		mat `bbar' = `bbar' * `ngg'
		mat colnames `bbar' = `colstr'
		mat colnames `vce' = `colstr'
		mat rownames `vce' = `colstr'
		if ("`res'" != ""){		
			capture confirm new variable `res'
			if _rc!=0{
				display as error  _col(2) "Variable `res' to hold residuals already exists."
				display as error  _col(2) "Either drop this variable or specify another name for residuals."
				exit
			}
			gen `res'=`r' if `touse'
		}
		if ("`pred'" != ""){		
			capture confirm new variable `pred'
			if _rc!=0{
				display as error  _col(2) "Variable `pred' to hold predicted values already exists."
				display as error  _col(2) "Either drop this variable or specify another name for residuals."
				exit
			}
			gen `pred'=`yhat' if `touse'
		}
		replace `r'=`r'^2 if `touse'
		sum `r' if `touse'
		scalar `sig2'=r(mean)

/* Standard errors */
		tempname names nn rnn tmp vce Xm stei tbetas Xn Xq
		local nn = rowsof(`bols')
		local names: colfullnames `bbar'
		scalar `rnn'=1/`nn'
		matrix `tmp'=`bols'-J(`nn',1,1)#(J(1,`nn',`rnn')*`bols')
		matrix coleq `tmp'=:
		matrix roweq `tmp'=:
		matrix `vce'=`tmp''*`tmp'/(`nn'*(`nn'-1))
		matrix rownames `vce'=`names'
		matrix colnames `vce'=`names'
		matrix colnames `vari'=`names'
		mata: `Xm' = st_matrix("`vari'")
		mata: `Xn' = st_matrix("`bols'")
		mata: `Xm' = sqrt(`Xm')
		mata: st_matrix("`stei'", `Xm')
		matrix colnames `stei'=`names'
		mata: `Xq' = `Xn' :/ `Xm'
		mata: st_matrix("`tbetas'", `Xq')
		matrix colnames `tbetas'=`names'

/* Robust means and DOLS/FMOLS inference */
	tempname rb rV temp est vcv names vcvm colms tb
	if ("`robust'" != ""){		
		local colms = colsof(`bols')
		svmat double `bols', names(`temp')
		matrix `rb'=J(1,`colms',0)
		matrix `rV'=J(`colms',`colms',0)
		local names: colfullnames `bols'
		forvalues i=1/`colms'{
			tempname `est'`i' `vcv'`i'
			rreg `temp'`i'
			matrix `est'`i'=e(b)
			matrix `vcv'`i'=e(V)
			matrix `rb'[1,`i']=`est'`i'[1,1]
			matrix `rV'[`i',`i']=`vcv'`i'[1,1]
		}
		matrix colnames `rb'=`names'
		matrix colnames `rV'=`names'
		matrix rownames `rV'=`names'
		mat `beta'=`rb'
		mat `vcvm'=`rV'
	}
	else {
		tempname beta2 sngg Xq tba tb tse tse2 rV2 colms vcvm beta
		if ("`fm'"!="" | "`dols'"!=""){
			local sngg = 1/sqrt(`ng')
			mata: `beta2' = st_matrix("`bbar'")
			mata: `Xq'   = st_matrix("`tbetas'")
			mata: `tba'  = colsum(`Xq')
			mata: `tb' = `tba' :* `sngg'
			mata: `tse' = `beta2' :/ `tb'
			mata: `tse2' = `tse' :^2
			mata: st_matrix("`tse2'", `tse2'')
			local names: colfullnames `bols'
			local colms = colsof(`bols')
			matrix `rV2'=J(`colms',`colms',0)
				forvalues i=1/`colms'{
				matrix `rV2'[`i',`i']=`tse2'[`i',1]
			}
			matrix colnames `bbar'=`names'
			matrix colnames `rV2'=`names'
			matrix rownames `rV2'=`names'
			matrix `vcvm'=`rV2'
		}
		else {
			mat `vcvm'=`vce'
		}
			mat `beta'=`bbar'
	}


/* Post the results */
		tempvar waldt 
		ereturn post `beta' `vcvm', obs(`nobs') depname(`depname') esample(`touse')
		capture test `ind', min `constant'
		ereturn scalar sigma=`sig2'
		if ("`trend'" != ""){
			local waldtrend = `waldtrend'/`ng'
			ereturn scalar trend_sig = `waldtrend'
			capture test `ind' `timevar'_t, min `constant'
			if _rc == 0 {
				ereturn scalar chi2 = r(chi2)
				ereturn scalar df_m = r(df)
				}
			else ereturn scalar df_m = 0
		}
		}
		capture test `ind', min `constant'
		if _rc == 0 {
			ereturn scalar chi2 = r(chi2)
			ereturn scalar df_m = r(df)
		}
		else  ereturn scalar df_m = 0
		ereturn scalar g_min  = `g1'
		ereturn scalar g_avg  = `g2'
		ereturn scalar g_max  = `g3'
		ereturn scalar N_g = `ng'
		ereturn matrix tbetas = `tbetas'
		ereturn matrix stebetas = `stei'
		ereturn matrix betas = `bols'
		ereturn matrix varbetas = `vari'
		ereturn local chi2type "Wald"
		if ("`augment'" != "" & "`impose'" != ""){
			ereturn local depvar "adjusted `depname'"
		}
		else	ereturn local depvar "`depname'"		
		if ("`augment'" != "") ereturn local title2 "AMG"
		if ("`cce'" != "") ereturn local title2 "CCEMG"
		if ("`augment'" == "" & "`cce'" == "") ereturn local title2 "MG"
		ereturn local title "Mean Group type estimation"
		ereturn local tvar "`tvar'"
		ereturn local ivar "`ivar'"
		ereturn local cmd "xtmg"



/* Make the result presentation more pleasing to the eye */

	mat bb=e(b)
	mat bbetas=e(betas)


	/* MG */
	if ("`augment'" == "" & "`cce'"=="" & "`fm'"=="" & "`dols'"=="") {
		if ("`trend'" != ""){
			mat colnames bb = `ind' trend _cons
			mat colnames bbetas = `ind' trend _cons
		}
		else{
			mat colnames bb = `ind' _cons 
			mat colnames bbetas = `ind' _cons
		}
	}

	/* GM FMOLS */
	if ("`fm'"!="") {
		if ("`trend'" != ""){
			mat colnames bb = `ind' trend _cons
			mat colnames bbetas = `ind' trend _cons
		}
		else{
			mat colnames bb = `ind' _cons 
			mat colnames bbetas = `ind' _cons
		}
	}

		/* GM DOLS */
	if ("`dols'"!="") {
		if ("`trend'" != ""){
			mat colnames bb = `ind' trend _cons
			mat colnames bbetas = `ind' trend _cons
		}
		else{
			mat colnames bb = `ind' _cons 
			mat colnames bbetas = `ind' _cons
		}
	}
	
	/* AMG */
	if ("`augment'" != "") {
		if ("`impose'" =="") {
			if ("`trend'" != ""){
				mat colnames bb = `ind' c_d_p trend _cons
				mat colnames bbetas = `ind' c_d_p trend _cons
				}
			else{
				mat colnames bb = `ind' c_d_p _cons 
				mat colnames bbetas = `ind' c_d_p _cons 
				}
			}
		else{
			if ("`trend'" != ""){
				mat colnames bb = `ind' trend _cons
				mat colnames bbetas = `ind' trend _cons
			}
			else{
				mat colnames bb = `ind' _cons 
				mat colnames bbetas = `ind' _cons 
			}
		}
	}


	/* CCEMG */
	if ("`cce'" != ""){
		local x
		foreach var of local ind{
			local x `x' `var'_avg
		}

		if ("`trend'" != ""){
			mat colnames bb = `ind' trend  `dep'_avg `x' _cons
			mat colnames bbetas = `ind' trend  `dep'_avg `x' _cons
		}
		else{
			mat colnames bb = `ind'  `dep'_avg `x' _cons
			mat colnames bbetas = `ind'  `dep'_avg `x' _cons
		}
	}

	ereturn repost b=bb, rename



display ""
display ""
if ("`dols'" != ""){
	display in gr "Pedroni (2001) Group Mean Dynamic OLS (DOLS) estimator"
display ""
}
if ("`fm'" != ""){
	display in gr "Pedroni (2000, 2001) Group Mean Fully-Modified OLS (FMOLS) estimator"
display ""
}
if ("`cce'" != ""){
	display in gr "Pesaran (2006) Common Correlated Effects Mean Group estimator"
display ""
}
if ("`augment'" != "") {
	display in gr "Augmented Mean Group estimator (Bond & Eberhardt, 2009; Eberhardt & Teal, 2010)"
	display ""
	if ("`impose'" !="") {
		display in gr as text "Common dynamic process " in ye "imposed" in gr " with unit coefficient " 
		display in gr "    Dependent variable" in ye " adjusted `depname'"
	}
	else {
		display in gr as text "Common dynamic process " in ye "included" in gr " as additional regressor" 
	}
}
if ("`augment'" == "" & "`cce'"=="" & "`dols'"=="" & "`fm'"=="") {
	display in gr "Pesaran & Smith (1995) Mean Group estimator"
	display ""
}
display in gr "All coefficients represent averages across groups (group variable: " in ye "`ivar'" in gr ")"
if ("`robust'" != ""){		
	display in gr "Coefficient averages computed as" in ye " outlier-robust" in gr " means (using rreg)" 
}
else {
	display in gr "Coefficient averages computed as " in ye "unweighted" in gr " means" 
}

_crcphdr
if ("`dols'" !="" | "`fm'" !="" ){
	di as text "{hline 13}{c TT}{hline 64}"
      	di as text %12s "`dep'" " {c |} " _c 
	di as text _col(21) ///
		`"Coef.   Std.Err.    panel-t   P>|z|   [`=strsubdp("`level'")'% Conf. Interval]"'	
        di as text "{hline 13}{c +}{hline 64}"
	tempname b se tstat pval cv coefs names name col sei 
      mat `b' = e(b)
      mat `se' = e(V)
      scalar `cv' = invnorm(1 - ((100-`level')/100)/2)
      local names : colnames bbetas
      local coefs : word count `names'
     local coefs = `coefs'-1
	  forvalues m = 1/`coefs' {
            	local col = 17
            	local name : word `m' of `names'
            	di as text %12s abbrev("`name'",12) " {c |}" as result _col(`col') %9.7g `b'[1,`m'] _c    
                  local col = `col' + 11
                  if (sqrt(`se'[`m',`m']) > 0 & sqrt(`se'[`m',`m']) < .) {
						di as res _col(`col') %9.7g sqrt(`se'[`m',`m']) "   " _c
						di as result %8.2f (`b'[1,`m']/sqrt(`se'[`m',`m'])) "   " _continue 
							scalar `pval'= 2*(1 - normal(abs(`b'[1,`m']/sqrt(`se'[`m',`m']))))
						di as result %5.3f `pval' "    " _c
						di as result %8.7g ( `b'[1,`m'] - `cv'*sqrt(`se'[`m',`m'])) "  " _continue
						di as result %8.7g ( `b'[1,`m'] + `cv'*sqrt(`se'[`m',`m'])) _continue
						di   
                  }
                  else {
                        di as text _col(36) ///
							".        .       .            .           ."
                  }
	}
	di as text "{hline 13}{c BT}{hline 64}"
}
else _coef_table, level(`level')
display in gr "Root Mean Squared Error (sigma): " in ye %5.4f sqrt(`sig2') 
if ("`robust'" != ""){
	display in gr "(RMSE uses residuals from group-specific regressions: unaffected by 'robust')."
	}
if ("`augment'" != "" & "`impose'"==""){
	display in gr "Variable " in ye "c_d_p" in gr " refers to the common dynamic process."
}
if ("`cce'" !=""){
	display in ye "Cross-section averaged regressors" in gr " are marked by the suffix " in ye "avg" in gr "." 
}
if ("`dols'" !=""){
	display in gr "DOLS lag-length selection: determined by " in ye "`dic'" in gr ". " 
	display in gr "DOLS panel t-statistics constructed as average t-statistics divided by N^(1/2). " 
}
if ("`fm'" !=""){
	display in gr "FMOLS panel t-statistic is the sum of all N t-statistics divided by N^(1/2). " 
}
if ("`fm'" !="" | "`dols'" !=""){
	display in gr "LR covariance estimated with " in ye "`kernel'" in gr " kernel, bandwidth selection: " in ye "`bmeth'" in gr "."
}
if ("`res'" != "") {
	display in gr "Residual series based on country regressions stored in variable: " in ye "`res'"
}
if ("`pred'" != "") {
	display in gr "Predicted values (yhat) based on country regressions stored in variable: " in ye "`pred'"
}
if ("`trend'" != ""){
	display in gr "Variable " in ye "trend" in gr " refers to the group-specific linear trend terms."
	display in gr "These deterministics are " in ye "not" in gr " part of the heterogeneous cointegrating equation."
	display in gr "Share of group-specific trends significant at " (100-`level') "% level: " ///
		in ye %5.2f `waldtrend' 
	}
*if ("`trend'" != "") & ("`cce'" != "") {
*	display in gr "Note that augmentation of the CCEMG estimator should account for the impact of"
*	display in gr "group-specific linear trends. The latter are unidentified."
*}
if ("`dols'" !="" | "`fm'" !="" ){
	display in gr " "
	display in gr "Note: The FMOLS/DOLS routines are based on" in ye " lrcov " in gr "and "  in ye "cointreg" in gr " written by "
	display in gr "      Wang Qunyong & Wu Na (Stata Journal, 2012). The usual disclaimers apply."
	display in gr " "
}
if ("`full'" != ""){
	di
	di
	di _col(25) "Group-specific coefficients"
	di as text "{hline 78}"
if ("`dols'" !="" | "`fm'" !="" ){
	di as text _col(21) ///
		`"Coef.   Std.Err.  panel-t   P>|z|     [`=strsubdp("`level'")'% Conf. Interval]"'	
	}
	else di as text _col(21) ///
		`"Coef.   Std. Err.      z    P>|z|    [`=strsubdp("`level'")'% Conf. Interval]"'	
	di as text "{hline 13}{c TT}{hline 64}"

	tempname b se tstat pval cv coefs names name col sei 
      mat `b' = e(betas)
      mat `se' = e(stebetas)
	  mat `tstat' = e(tbetas)	
      scalar `cv' = invnorm(1 - ((100-`level')/100)/2)
      local names : colnames bbetas
      local coefs : word count `names'
      forvalues i = 1/`ng'{
      	di as text %12s "Group `i'" " {c |} "  
        di as text "{hline 13}{c +}{hline 64}"
		forvalues m = 1/`coefs' {
            	local col = 17
            	local name : word `m' of `names'
            	di as text %12s abbrev("`name'",12) " {c |}" as result _col(`col') %9.7g `b'[`i',`m'] _c    
                  local col = `col' + 11
                  if (`se'[`i',`m'] > 0 & `se'[`i',`m'] < .) {
                  	di as res _col(`col') %9.7g `se'[`i',`m'] "   " _c
                  	di as result %6.2f `tstat'[`i',`m'] "   " _continue 
                 		scalar `pval'= 2*(1 - normal(abs(`tstat'[`i',`m'])))
				di as result %5.3f `pval' "    " _c
				di as result %9.7g ( `b'[`i',`m'] - `cv'*`se'[`i',`m']) "   " _continue
				di as result %9.7g ( `b'[`i',`m'] + `cv'*`se'[`i',`m']) _continue
                  	di   
                  }
                  else {
                        di as text _col(36) ///
				".        .       .            .           ."
                  }
           }
	if (`i' < `ng') {
      	di as text "{hline 13}{c +}{hline 64}"
      }
	else {
      	di as text "{hline 13}{c BT}{hline 64}"
		}
}

}


end

exit

