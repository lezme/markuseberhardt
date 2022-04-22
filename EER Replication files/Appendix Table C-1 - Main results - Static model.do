******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*          Appendix Table C-1 - Main results - Static specification 
******************************************************************************
*                       Created: 26th February 2021
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

dropbox, nocd

global path0 "/Users/lezme/Dropbox"
global path0 "D:/Dropbox"

global path "$path0/Curious GMM/Acemoglu Replication"

global dofolder "$path/do_files"
global otherdata "$path/Data"
global outputfolder "$path/Output"
global rubbishfolder "$path/Rubbish"
global texfolder "$path/Tex"

global xtmgpath = "/Users/lezme/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"
global xtmgpath = "D:/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"

use "$path/DDCGdata_final.dta", clear

***
*** ANRR - Panel (a)
***
preserve 
xtreg y dem yy*, fe
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum 
sort wbcode
by wbcode: egen democ_sample = sum(demevent) if e(sample)
by wbcode: egen autoc_sample = sum(revevent) if e(sample)
by wbcode: egen always_dem = sum(dem) if e(sample)
replace always_dem = always_dem/csum
gen never_dem = 1 if always_dem==0 & e(sample)
replace never_dem = 0 if missing(never_dem) & e(sample)
replace always_dem = 0 if always_dem!=1 & e(sample)
* always_dem indicates the countries which are democracies throughout the sample
gen dem_sample = 0 if e(sample)
replace dem_sample = 1 if (democ_sample!=0 | autoc_sample!=0) & e(sample)
replace dem_sample = 0 if (wbcode=="POL" | wbcode=="SLV" | wbcode=="EST" | wbcode=="MKD")
* always democracy in this sample OR (SLV, MKD) the x-vars are missing so no pre-event
* data or the event itself is not in the sample period
* VEN, UGA, KHM, GMB are reversals only (4 reversals only)
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

xtmg y dem ginv tradewb ///
		yNT ///
		ginvNT   ///
		tradewbNT  ///
if dem_sample==1 
drop if !e(sample)
tab year if e(sample)
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum demevent
tab wbcode demevent if csum!=0 & !missing(csum) & c==1
tab wbcode revevent if csum!=0 & !missing(csum) & c==1
tab wbcode dem if csum!=0 & !missing(csum) & c==1

mat list e(betas)

run "$xtmgpath/xtmg.ado"

*** Full sample results

*** 105 demevents and 58 reversals (ANRR) in the full 83 country-sample

estimates clear 

local i = 1

** Plain vanilla MG - Column (1)
eststo: xtmg y dem  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if dem==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries1=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok - Column (2)
eststo: xtmg y dem  ///
		yNT ///
		ginvNT  ///
		tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if dem==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries1=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** MG with covariates - Column (3)
eststo: xtmg y dem ginv tradewb ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if dem==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries1=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates - Column (4) 
eststo: xtmg y dem ginv tradewb ///
		yNT ///
		ginvNT   ///
		tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if dem==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries1=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem _cons)"
esttab model* using "$texfolder/20220421_Fullsample_Static_DID_ANRR.tex", b(3) se(3) abs `esttabformat' `esttabstats'  style(tex) replace 
restore 


***
*** Boix et al - Panel (b)
***
preserve 
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
replace dem_sample = 0 if (wbcode=="MDA" | wbcode=="MKD" | wbcode=="SVN")
tab wbcode2 demevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
* 84 and 47
tab wbcode2 demBMR if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)

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
		yNT ///
		ginvNT   ///
		tradewbNT   ///
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
	
*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG - Column (1)
eststo: xtmg y demBMR  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demBMR==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok - Column (2)
eststo: xtmg y demBMR  ///
		yNT ///
		ginvNT  ///
		tradewbNT   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demBMR==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** MG with covariates - Column (3)
eststo: xtmg y demBMR ginv tradewb ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demBMR==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates - Column (4)
eststo: xtmg y demBMR ginv tradewb ///
		yNT ///
		ginvNT   ///
		tradewbNT  ///
if dem_sample==1 & c==1, robust
mat bBMR = e(betas)
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demBMR==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

label var demBMR "Democracy (BMR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(demBMR _cons _n_1)"
esttab model* using "$texfolder/20220421_Fullsample_Static_DID_BMR.tex", b(3) se(3) abs `esttabformat' `esttabstats'  style(tex) replace 
restore 
	


***
*** Cheibub et al - Panel (c)
***
preserve 
drop demevent revevent 
sort wbcode2 year
gen demevent = 0
by wbcode2: replace demevent = 1 if l.demCGV==0 & demCGV==1
gen revevent = 0
by wbcode2: replace revevent = 1 if l.demCGV==1 & demCGV==0
xtreg y demCGV yy* tradewb ginv, fe
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum
sort wbcode
by wbcode: egen democ_sample = sum(demevent) if e(sample)
by wbcode: egen autoc_sample = sum(revevent) if e(sample)
by wbcode: egen always_dem = sum(demCGV) if e(sample)
replace always_dem = always_dem/csum
gen never_dem = 1 if always_dem==0 & e(sample)
replace never_dem = 0 if missing(never_dem) & e(sample)
replace always_dem = 0 if always_dem!=1 & e(sample)
* always_dem indicates the countries which are democracies throughout the sample
gen dem_sample = 0 if e(sample)
replace dem_sample = 1 if (democ_sample!=0 | autoc_sample!=0) & e(sample)
tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
* 95 and 50
tab wbcode2 demCGV if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)

replace dem_sample = 0 if (wbcode=="ARM" | wbcode=="LTU" | wbcode=="LVA" | wbcode=="MDA" | wbcode=="MKD" | wbcode=="SVN")
* if a country has no observations in autoc or democ (usually due to trade and ginv data!)
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
		yNT ///
		ginvNT   ///
		tradewbNT   ///
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
	
*** Full sample results

estimates clear 

local i = 1

** Plain vanilla MG - Column (1)
eststo: xtmg y demCGV  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok - Column (2)
eststo: xtmg y demCGV  ///
		yNT ///
		ginvNT ///
		tradewbNT   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** MG with covariates - Column (3)
eststo: xtmg y demCGV ginv tradewb ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates - Column (4) 
eststo: xtmg y demCGV ginv tradewb ///
		yNT ///
		ginvNT   ///
		tradewbNT    ///
if dem_sample==1 & c==1, robust
mat bCGV = e(betas)
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent==1 & e(sample)
scalar Events = r(N)
count if revevent==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

label var demCGV "Democracy (CGV)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(demCGV _cons _n_1)"
esttab model* using "$texfolder/20220421_Fullsample_Static_DID_CGV.tex", b(3) se(3) abs `esttabformat' `esttabstats'  style(tex) replace 
restore 


***
*** P&S - Panel (d)
***
preserve 
xtreg y demPS yy* tradewb ginv, fe
gen c= 1 if e(sample)
sort wbcode year 
by wbcode: egen csum=sum(c)
tab csum
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
xtreg y demPS yy* if dem_sample==1, fe 
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
		yNT ///
		ginvNT   ///
		tradewbNT  ///
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
	
*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG - Column (1)
eststo: xtmg y demPS  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demPS==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok - Column (2)
eststo: xtmg y demPS  ///
		yNT ///
		ginvNT  ///
		tradewbNT  ///
if dem_sample==1 &  c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demPS==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** MG with covariates - Column (3)
eststo: xtmg y demPS ginv tradewb ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demPS==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates - Column (4) 
eststo: xtmg y demPS ginv tradewb ///
		yNT ///
		ginvNT   ///
		tradewbNT    ///
if dem_sample==1 & c==1, robust
 
mat bPS = e(betas)
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demPS==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
svmat betas, name(b1_)
replace b1_1=. if b1_1==0
eststo: rreg b1_1
estadd scalar countries=e(N)
drop b1_*
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/e(N)
est store model`i'

label var demPS "Democracy (PS)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem _cons _n_1)"
esttab model* using "$texfolder/20220421_Fullsample_Static_DID_PS.tex", b(3) se(3) abs `esttabformat' `esttabstats'  style(tex) replace 
restore 

 
exit