******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                Appendix Figure D-2 - Sample reduction 
******************************************************************************
*                       Created: 25th February 2021
*                       Revised: 21st April 2022
******************************************************************************
*                          Markus Eberhardt
******************************************************************************
*            School of Economics, University of Nottingham,
*            University Park, Nottingham NG7 2RD, England
*              email: markus.eberhardt@nottingham.ac.uk
*            web: https://sites.google.com/site/medevecon/
******************************************************************************


clear all 
set matsize 11000

global path0 "D:/Dropbox"
global path0 "/Users/lezme/Dropbox"

global pspath "$path0/Curious GMM/Papaioannou replication/Data"
global path "$path0/Curious GMM/Acemoglu Replication"

global dofolder "$path/do_files"
global otherdata "$path/Data"
global outputfolder "$path/Output"
global rubbishfolder "$path/Rubbish"
global texfolder "$path/Tex"

global xtmgpath = "D:/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"
global xtmgpath = "/Users/lezme/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"


***
*** BMR - Panel (a)
***

use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 
xtreg y demBMR yy* tradewb ginv, fe
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum
sort wbcode
gen demevent2 = 0
by wbcode: replace demevent2 = 1 if demBMR[_n-1]==0 & demBMR==1
gen revevent2 = 0
by wbcode: replace revevent2 = 1 if demBMR[_n-1]==1 & demBMR==0
by wbcode: egen democ_sample = sum(demevent2) if e(sample)
by wbcode: egen autoc_sample = sum(revevent2) if e(sample)
by wbcode: egen always_dem = sum(demBMR) if e(sample)
replace always_dem = always_dem/csum
gen never_dem = 1 if always_dem==0 & e(sample)
replace never_dem = 0 if missing(never_dem) & e(sample)
replace always_dem = 0 if always_dem!=1 & e(sample)
* always_dem indicates the countries which are democracies throughout the sample
gen dem_sample = 0 if e(sample)
replace dem_sample = 1 if (democ_sample!=0 | autoc_sample!=0) & e(sample)
replace dem_sample = 0 if (wbcode=="COG" | wbcode=="CYP")
* if a country has no democratisation or reversal it is either always a demo or never
xtreg y demBMR yy* tradewb ginv if dem_sample==1, fe 
* dem_sample indicates countries which did *switch* into democracy or autocracy during the sample
drop c csum

label var never_dem "Country never a democracy during sample period"
label var always_dem "Country always a democracy during sample period"
label var dem_sample "Country transitioned into or from democracy during sample period"


sort year
by year: egen yT = mean(y) if never_dem==1
by year: egen ginvT = mean(ginv) if never_dem==1
by year: egen tradewbT = mean(tradewb) if never_dem==1

by year: egen yNT = mean(yT)
by year: egen ginvNT = mean(ginvT)
by year: egen tradewbNT = mean(tradewbT)

sort wbcode2 year
tsset wbcode2 year
gen dginv = d.ginv
gen dtradewb = d.tradewb
gen dy = d.y
gen ddemBMR = d.demBMR

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddemBMR dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	


xtmg y demBMR ginv tradewb ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1
drop if !e(sample)
tab year if e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
mat list e(betas)
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum

run "$xtmgpath/xtmg.ado"

matrix CKdyn=J(16,14,.)

local i = 1
cls
tsset wbcode2 year
forvalues k=2007(-1)1992{
	**** Plain vanilla model
	qui: xtmg y demBMR  ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
		if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',1]=`k'
	* Total number of observations
	mat CKdyn[`i',2]=e(N)
	* Number of Countries
	mat CKdyn[`i',3]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (plain vanilla)
	mat CKdyn[`i',4]=b_1
	mat CKdyn[`i',5]=se_1	

	**** Plain vanilla model (MG)
	qui: xtmg y demBMR  ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (MG)
	mat CKdyn[`i',6]=b_1
	mat CKdyn[`i',7]=se_1	
	
	**** Model with covariates
	qui: xtmg y demBMR ginv tradewb ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
	if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',8]=`k'
	* Total number of observations
	mat CKdyn[`i',9]=e(N)
	* Number of Countries
	mat CKdyn[`i',10]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (add covariates)
	mat CKdyn[`i',11]=b_1
	mat CKdyn[`i',12]=se_1	

	**** Model with covariates (MG)
	qui: xtmg y demBMR ginv tradewb ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (MG w/ covariates)
	mat CKdyn[`i',13]=b_1
	mat CKdyn[`i',14]=se_1	

	di in gr "Result #" in ye `i' in gr ": Completed results for end year " in ye `k' in gr ". N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list CKdyn
svmat CKdyn, names(CKdyn_)
rename CKdyn_1 CKdyn_endyear
rename CKdyn_2 CKdyn_obs
rename CKdyn_3 CKdyn_N
rename CKdyn_4 CKdyn_b_lr_dem
rename CKdyn_5 CKdyn_se_lr_dem
rename CKdyn_6 CKdyn_b_dem_mg
rename CKdyn_7 CKdyn_se_dem_mg
rename CKdyn_8 CKdyn_endyear_cov
rename CKdyn_9 CKdyn_obs_cov
rename CKdyn_10 CKdyn_N_cov
rename CKdyn_11 CKdyn_b_lr_dem_cov
rename CKdyn_12 CKdyn_se_lr_dem_cov
rename CKdyn_13 CKdyn_b_lr_dem_covmg
rename CKdyn_14 CKdyn_se_lr_dem_covmg

preserve
keep CKdyn_*
gen model = _n if !missing(CKdyn_endyear)
drop if _n>30
order model
sort CKdyn_endyear
save  "$otherdata/202010308_CK_lr_covar_reduction_time_BMR.dta", replace
restore


*** Graph BMR (Included in the paper)

preserve 
use "$otherdata/202010308_CK_lr_covar_reduction_time_BMR.dta", clear

gen CKdyn_t_lr_dem = CKdyn_b_lr_dem/CKdyn_se_lr_dem
gen CKdyn_t_lr_dem_mg = CKdyn_b_dem_mg/CKdyn_se_dem_mg
gen CKdyn_t_lr_dem_cov = CKdyn_b_lr_dem_cov/CKdyn_se_lr_dem_cov
gen CKdyn_t_lr_dem_covmg = CKdyn_b_lr_dem_covmg/CKdyn_se_lr_dem_covmg

twoway ///	
	(connected CKdyn_b_lr_dem_covmg CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(navy) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(dash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_covmg CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_mg)>1.645, ///
		msymbol(O) lcolor(navy) msize(medlarge) mfcolor(navy) mlcolor(black) mlwidth(thin)) ///
	(connected CKdyn_b_lr_dem CKdyn_endyear, ///
		lw(.6) msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(shortdash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem CKdyn_endyear if abs(CKdyn_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(midblue) mlcolor(black) mlwidth(thin)) ///		
	(connected CKdyn_b_lr_dem_cov CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(solid) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_cov)>1.645, ///
		msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(pink*1.25) mlcolor(black) mlwidth(thin)) ///	
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov, yaxis(2) msymbol(none))  ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Sample End Year", size(medsmall)) ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(2(2)12, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small))  ///
	yscale(range(1 12.1)) ///
	yscale(range(1 12.1) axis(2)) ymtick(##2) ymtick(##2, axis(2)) ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)  axis(2))  ///
	ylabel(2(2)12, labsize(small) grid glc(black) glw(.2) glp(shortdash) angle(0) axis(2)) /// 
	xlabel(1992(1)2007) xscale(range(1991.5 2007.5) reverse) xmtick(##1)   ///
	legend(rows(2) order(- "Average LR Democracy Effect (BMR):" 3 1 5 - "Statistically significant (10% level):" 4 2  6) ///
	label(1 "MG w/ covariates") ///
	label(2 "") ///
	label(3 "C&K MG") ///
	label(4 "") ///
	label(5 "C&K MG w/ covariates") ///
	label(6 "") ///
	size(medsmall)) 
graph export "$texfolder/20220421_CKdyn_evolution_time_BMR.eps", as(eps) replace

list CKdyn_endyear_cov CKdyn_b_lr_dem_cov CKdyn_se_lr_dem_cov CKdyn_t_lr_dem_cov if (CKdyn_endyear_cov==1999  ///
	| CKdyn_endyear_cov==1992), noobs separator(0)
	
restore 


***
*** Cheibub et al (CGV) - Panel (b)
***

use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 
xtreg y demCGV yy* tradewb ginv, fe
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum
sort wbcode
sort wbcode2 year
gen demevent2 = 0
by wbcode2: replace demevent2 = 1 if (l.demCGV==0 & demCGV==1)
gen revevent2 = 0
by wbcode2: replace revevent2 = 1 if (l.demCGV==1 & demCGV==0)
by wbcode2: egen democ_sample = sum(demevent2) if e(sample)
by wbcode2: egen autoc_sample = sum(revevent2) if e(sample)
by wbcode2: egen always_dem = sum(demCGV) if e(sample)
replace always_dem = always_dem/csum
gen never_dem = 1 if always_dem==0 & e(sample)
replace never_dem = 0 if missing(never_dem) & e(sample)
replace always_dem = 0 if always_dem!=1 & e(sample)
* always_dem indicates the countries which are democracies throughout the sample
gen dem_sample = 0 if e(sample)
replace dem_sample = 1 if (democ_sample!=0 | autoc_sample!=0) & e(sample)
replace dem_sample = 0 if (wbcode=="BTN" | wbcode=="MRT")
* if a country has no democratisation or reversal it is either always a demo or never
xtreg y demCGV yy* tradewb ginv if dem_sample==1, fe 
* dem_sample indicates countries which did *switch* into democracy or autocracy during the sample
drop c csum

label var never_dem "Country never a democracy during sample period"
label var always_dem "Country always a democracy during sample period"
label var dem_sample "Country transitioned into or from democracy during sample period"

sort year
by year: egen yT = mean(y) if never_dem==1
by year: egen ginvT = mean(ginv) if never_dem==1
by year: egen tradewbT = mean(tradewb) if never_dem==1

by year: egen yNT = mean(yT)
by year: egen ginvNT = mean(ginvT)
by year: egen tradewbNT = mean(tradewbT)

sort wbcode2 year
tsset wbcode2 year
gen dginv = d.ginv
gen dtradewb = d.tradewb
gen dy = d.y
gen ddemCGV = d.demCGV

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddemCGV dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	


xtmg y demCGV ginv tradewb ///
		l0ddemCGV l1ddemCGV l2ddemCGV  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1, robust
drop if !e(sample)
tab year if e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
mat list e(betas)
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum

run "$xtmgpath/xtmg.ado"

matrix CKdyn=J(16,14,.)

local i = 1
cls
tsset wbcode2 year
forvalues k=2008(-1)1993{
	**** Plain vanilla model
	qui: xtmg y demCGV  ///
		l0ddemCGV l1ddemCGV l2ddemCGV ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
		if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',1]=`k'
	* Total number of observations
	mat CKdyn[`i',2]=e(N)
	* Number of Countries
	mat CKdyn[`i',3]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (plain vanilla)
	mat CKdyn[`i',4]=b_1
	mat CKdyn[`i',5]=se_1	

	**** Plain vanilla model (MG)
	qui: xtmg y demCGV  ///
		l0ddemCGV l1ddemCGV l2ddemCGV ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (MG)
	mat CKdyn[`i',6]=b_1
	mat CKdyn[`i',7]=se_1	
	
	**** Model with covariates
	qui: xtmg y demCGV ginv tradewb ///
		l0ddemCGV l1ddemCGV l2ddemCGV ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
	if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',8]=`k'
	* Total number of observations
	mat CKdyn[`i',9]=e(N)
	* Number of Countries
	mat CKdyn[`i',10]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (add covariates)
	mat CKdyn[`i',11]=b_1
	mat CKdyn[`i',12]=se_1	

	**** Model with covariates (MG)
	qui: xtmg y demCGV ginv tradewb ///
		l0ddemCGV l1ddemCGV l2ddemCGV ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (MG w/ covariates)
	mat CKdyn[`i',13]=b_1
	mat CKdyn[`i',14]=se_1	

	di in gr "Result #" in ye `i' in gr ": Completed results for end year " in ye `k' in gr ". N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list CKdyn
svmat CKdyn, names(CKdyn_)
rename CKdyn_1 CKdyn_endyear
rename CKdyn_2 CKdyn_obs
rename CKdyn_3 CKdyn_N
rename CKdyn_4 CKdyn_b_lr_dem
rename CKdyn_5 CKdyn_se_lr_dem
rename CKdyn_6 CKdyn_b_dem_mg
rename CKdyn_7 CKdyn_se_dem_mg
rename CKdyn_8 CKdyn_endyear_cov
rename CKdyn_9 CKdyn_obs_cov
rename CKdyn_10 CKdyn_N_cov
rename CKdyn_11 CKdyn_b_lr_dem_cov
rename CKdyn_12 CKdyn_se_lr_dem_cov
rename CKdyn_13 CKdyn_b_lr_dem_covmg
rename CKdyn_14 CKdyn_se_lr_dem_covmg

preserve
keep CKdyn_*
gen model = _n if !missing(CKdyn_endyear)
drop if _n>30
order model
sort CKdyn_endyear
save  "$otherdata/202010308_CK_lr_covar_reduction_time_CGV.dta", replace
restore

*** Graph CGV (Included in the paper)

preserve 
use "$otherdata/202010308_CK_lr_covar_reduction_time_CGV.dta", clear

gen CKdyn_t_lr_dem = CKdyn_b_lr_dem/CKdyn_se_lr_dem
gen CKdyn_t_lr_dem_mg = CKdyn_b_dem_mg/CKdyn_se_dem_mg
gen CKdyn_t_lr_dem_cov = CKdyn_b_lr_dem_cov/CKdyn_se_lr_dem_cov
gen CKdyn_t_lr_dem_covmg = CKdyn_b_lr_dem_covmg/CKdyn_se_lr_dem_covmg

twoway ///	
	(connected CKdyn_b_lr_dem_covmg CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(navy) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(dash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_covmg CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_mg)>1.645, ///
		msymbol(O) lcolor(navy) msize(medlarge) mfcolor(navy) mlcolor(black) mlwidth(thin)) ///
	(connected CKdyn_b_lr_dem CKdyn_endyear, ///
		lw(.6) msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(shortdash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem CKdyn_endyear if abs(CKdyn_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(midblue) mlcolor(black) mlwidth(thin)) ///		
	(connected CKdyn_b_lr_dem_cov CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(solid) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_cov)>1.645, ///
		msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(pink*1.25) mlcolor(black) mlwidth(thin)) ///	
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov, yaxis(2) msymbol(none))  ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Sample End Year", size(medsmall)) ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(2(2)10, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small))  ///
	yscale(range(1.9 10.1)) ///
	yscale(range(1.9 10.1) axis(2)) ymtick(##2) ymtick(##2, axis(2))  ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)  axis(2))  ///
	ylabel(2(2)10, labsize(small) grid glc(black) glw(.2) glp(shortdash) angle(0) axis(2)) /// 
	xlabel(1993(1)2008) xscale(range(1992.5 2008.5) reverse) xmtick(##1) yscale(range(.5 10.5) axis(2)) yscale(range(.5 10.5)) ///
	legend(rows(2) order(- "Average LR Democracy Effect (CGV):" 3 1 5 - "Statistically significant (10% level):" 4 2  6) ///
	label(1 "MG w/ covariates") ///
	label(2 "") ///
	label(3 "C&K MG") ///
	label(4 "") ///
	label(5 "C&K MG w/ covariates") ///
	label(6 "") ///
	size(medsmall)) 
graph export "$texfolder/20220421_CKdyn_evolution_time_CGV.eps", as(eps) replace

list CKdyn_endyear_cov CKdyn_b_lr_dem_cov CKdyn_se_lr_dem_cov CKdyn_t_lr_dem_cov if (CKdyn_endyear_cov==1993  ///
	| CKdyn_endyear_cov==2007), noobs separator(0)
	
restore 


***
*** P&S - Panel (c)
***

use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 
xtreg y demPS yy* tradewb ginv, fe
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum
sort wbcode
gen demevent2 = 0
by wbcode: replace demevent2 = 1 if demPS[_n-1]==0 & demPS==1
gen revevent2 = 0
by wbcode: replace revevent2 = 1 if demPS[_n-1]==1 & demPS==0
sort wbcode
by wbcode: egen democ_sample = sum(demevent2) if e(sample)
by wbcode: egen autoc_sample = sum(revevent2) if e(sample)
by wbcode: egen always_dem = sum(demPS) if e(sample)
replace always_dem = always_dem/csum
gen never_dem = 1 if always_dem==0 & e(sample)
replace never_dem = 0 if missing(never_dem) & e(sample)
replace always_dem = 0 if always_dem!=1 & e(sample)
* always_dem indicates the countries which are democracies throughout the sample
gen dem_sample = 0 if e(sample)
replace dem_sample = 1 if (democ_sample!=0 | autoc_sample!=0) & e(sample)
*replace dem_sample = 0 if (wbcode=="POL" | wbcode=="BDI" | wbcode=="BFA" | wbcode=="CAF" | wbcode=="CIV" | wbcode=="COG" | wbcode=="FJI" | wbcode=="GMB" | wbcode=="GNB" | wbcode=="KEN" | wbcode=="NER" | wbcode=="NPL" | wbcode=="PAK" | wbcode=="SDN")
*replace never_dem = 1 if (wbcode=="POL" | wbcode=="BDI" | wbcode=="BFA" | wbcode=="CAF" | wbcode=="CIV" | wbcode=="COG" | wbcode=="FJI" | wbcode=="GMB" | wbcode=="GNB" | wbcode=="KEN" | wbcode=="NER" | wbcode=="NPL" | wbcode=="PAK" | wbcode=="SDN")
tab wbcode2 demevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 demPS if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)

* if a country has no democratisation or reversal it is either always a demo or never
xtreg y demPS yy* tradewb ginv if dem_sample==1, fe 
* dem_sample indicates countries which did *switch* into democracy or autocracy during the sample
drop c csum

label var never_dem "Country never a democracy during sample period"
label var always_dem "Country always a democracy during sample period"
label var dem_sample "Country transitioned into or from democracy during sample period"

sort year
by year: egen yT = mean(y) if never_dem==1
by year: egen ginvT = mean(ginv) if never_dem==1
by year: egen tradewbT = mean(tradewb) if never_dem==1
by year: egen lmortT = mean(lmort) if never_dem==1 
by year: egen lgovT = mean(lgov) if never_dem==1 

by year: egen yNT = mean(yT)
by year: egen ginvNT = mean(ginvT)
by year: egen tradewbNT = mean(tradewbT)
by year: egen lmortNT = mean(lmortT)
by year: egen lgovNT = mean(lgovT)

sort wbcode2 year

tsset wbcode2 year
gen dginv = d.ginv
gen dtradewb = d.tradewb
gen dy = d.y
gen ddemPS = d.demPS

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddemPS dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	


xtmg y demPS ginv tradewb ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 
drop if !e(sample)
tab year if e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
mat list e(betas)
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab year if e(sample)

run "$xtmgpath/xtmg.ado"

matrix CKdyn=J(16,14,.)

local i = 1
cls
tsset wbcode2 year
forvalues k=2010(-1)1995{
	**** Plain vanilla model
	qui: xtmg y demPS  ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
		if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',1]=`k'
	* Total number of observations
	mat CKdyn[`i',2]=e(N)
	* Number of Countries
	mat CKdyn[`i',3]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (plain vanilla)
	mat CKdyn[`i',4]=b_1
	mat CKdyn[`i',5]=se_1	

	**** Plain vanilla model (MG)
	qui: xtmg y demPS  ///
		l0ddemPS l1ddemPS l2ddemPS  ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* LR Dem coefficient (MG)
	mat CKdyn[`i',6]=b_1
	mat CKdyn[`i',7]=se_1	
	
	**** Model with covariates
	qui: xtmg y demPS ginv tradewb ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
	if dem_sample==1 & year<=`k' & c==1, robust
	* Min number of observations
	mat CKdyn[`i',8]=`k'
	* Total number of observations
	mat CKdyn[`i',9]=e(N)
	* Number of Countries
	mat CKdyn[`i',10]=e(N_g)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (add covariates)
	mat CKdyn[`i',11]=b_1
	mat CKdyn[`i',12]=se_1	

	**** Model with covariates (MG)
	qui: xtmg y demPS ginv tradewb ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
	if e(sample), robust
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	* Dem coefficient (MG w/ covariates)
	mat CKdyn[`i',13]=b_1
	mat CKdyn[`i',14]=se_1	

	di in gr "Result #" in ye `i' in gr ": Completed results for end year " in ye `k' in gr ". N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list CKdyn
svmat CKdyn, names(CKdyn_)
rename CKdyn_1 CKdyn_endyear
rename CKdyn_2 CKdyn_obs
rename CKdyn_3 CKdyn_N
rename CKdyn_4 CKdyn_b_lr_dem
rename CKdyn_5 CKdyn_se_lr_dem
rename CKdyn_6 CKdyn_b_dem_mg
rename CKdyn_7 CKdyn_se_dem_mg
rename CKdyn_8 CKdyn_endyear_cov
rename CKdyn_9 CKdyn_obs_cov
rename CKdyn_10 CKdyn_N_cov
rename CKdyn_11 CKdyn_b_lr_dem_cov
rename CKdyn_12 CKdyn_se_lr_dem_cov
rename CKdyn_13 CKdyn_b_lr_dem_covmg
rename CKdyn_14 CKdyn_se_lr_dem_covmg

preserve
keep CKdyn_*
gen model = _n if !missing(CKdyn_endyear)
drop if _n>30
order model
sort CKdyn_endyear
save  "$otherdata/202010308_CK_lr_covar_reduction_time_PS.dta", replace
restore

*** Graph PS (Included in the paper)

preserve 
use "$otherdata/202010308_CK_lr_covar_reduction_time_PS.dta", clear

gen CKdyn_t_lr_dem = CKdyn_b_lr_dem/CKdyn_se_lr_dem
gen CKdyn_t_lr_dem_mg = CKdyn_b_dem_mg/CKdyn_se_dem_mg
gen CKdyn_t_lr_dem_cov = CKdyn_b_lr_dem_cov/CKdyn_se_lr_dem_cov
gen CKdyn_t_lr_dem_covmg = CKdyn_b_lr_dem_covmg/CKdyn_se_lr_dem_covmg

twoway ///	
	(connected CKdyn_b_lr_dem_covmg CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(navy) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(dash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_covmg CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_mg)>1.645, ///
		msymbol(O) lcolor(navy) msize(medlarge) mfcolor(navy) mlcolor(black) mlwidth(thin)) ///
	(connected CKdyn_b_lr_dem CKdyn_endyear, ///
		lw(.6) msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(shortdash) lw(.8)) ///
	(scatter CKdyn_b_lr_dem CKdyn_endyear if abs(CKdyn_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(midblue) msize(medlarge) mfcolor(midblue) mlcolor(black) mlwidth(thin)) ///		
	(connected CKdyn_b_lr_dem_cov CKdyn_endyear_cov, ///
		lw(.6) msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin) lpattern(solid) lw(.8)) ///
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov if abs(CKdyn_t_lr_dem_cov)>1.645, ///
		msymbol(O) lcolor(pink*1.25) msize(medlarge) mfcolor(pink*1.25) mlcolor(black) mlwidth(thin)) ///	
	(scatter CKdyn_b_lr_dem_cov CKdyn_endyear_cov, yaxis(2) msymbol(none))  ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Sample End Year", size(medsmall)) ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(6(2)14, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small))  ///
	yscale(range(4.5 14.5)) ///
	yscale(range(4.5 14.5) axis(2)) ymtick(##2) ymtick(##2, axis(2)) ///
	ytitle("Average Long-run Coefficient on Democracy", size(medsmall)  axis(2))  ///
	ylabel(6(2)14, labsize(small) grid glc(black) glw(.2) glp(shortdash) angle(0) axis(2)) /// 
	xlabel(1995(1)2010) xscale(range(1994.5 2010.5) reverse) xmtick(##1)   ///
	legend(rows(2) order(- "Average LR Democracy Effect (PS):" 3 1 5 - "Statistically significant (10% level):" 4 2  6) ///
	label(1 "MG w/ covariates") ///
	label(2 "") ///
	label(3 "C&K MG") ///
	label(4 "") ///
	label(5 "C&K MG w/ covariates") ///
	label(6 "") ///
	size(medsmall)) 
graph export "$texfolder/20220421_CKdyn_evolution_time_PS.eps", as(eps) replace

list CKdyn_endyear_cov CKdyn_b_lr_dem_cov CKdyn_se_lr_dem_cov CKdyn_t_lr_dem_cov if (CKdyn_endyear_cov==2002), noobs separator(0)
	
restore 
