******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                       Table 2 - Alternative Samples 
******************************************************************************
*                       Created: 3rd March 2021
*                       Revised: 21st April 2022
******************************************************************************
*                          Markus Eberhardt
******************************************************************************
*            School of Economics, University of Nottingham,
*            University Park, Nottingham NG7 2RD, England
*              email: markus.eberhardt@nottingham.ac.uk
*            web: https://sites.google.com/site/medevecon/
******************************************************************************

* Note: Results are provided as .tex files separately for each panel (a-e)

clear all 
set matsize 11000

global path0 "/Users/lezme/Dropbox"
global path0 "D:/Dropbox"

global pspath "$path0/Curious GMM/Papaioannou replication/Data"
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
*** ANRR definition of democracy (full sample): Panel (a)
***
preserve 
tsset wbcode2 year
xtreg y dem yy* tradewb ginv, fe
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

* if a country has no democratisation or reversal it is either always a demo or never
xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
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
gen ddem = d.dem

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddem dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  
drop if !e(sample)
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimate identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum dem
tabstat csum dem, by(wbcode)

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG (Column 1)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok (Column 2)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** MG with covariates (Column 3)
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
*mat betas = e(betas)
*svmat betas, name(b1_)
*replace b1_1=. if b1_1==0
*eststo: rreg b1_1
*drop b1_*
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates (Column 4) 
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem)"
esttab model* using "$texfolder/20220421_Fullsample_DID_ANRR.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
restore 



***
*** ANRR definition of democracy (Drop reversal-only countries): Panel (b)
***

preserve 
tsset wbcode2 year
xtreg y dem yy* tradewb ginv, fe
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
* Robustness check without reversal-only countries
replace dem_sample = 0 if (wbcode=="GMB" | wbcode=="VEN" | wbcode=="ZWE" | wbcode=="UGA")

* if a country has no democratisation or reversal it is either always a demo or never
xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
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
gen ddem = d.dem

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddem dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  
drop if !e(sample)
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimate identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum dem
tabstat csum dem, by(wbcode)

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG (Column 1)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok (Column 2)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** MG with covariates (Column 3)
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates (Column 4) 
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem)"
esttab model* using "$texfolder/20220421_Fullsample_norev_DID_ANRR.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
restore 

***
*** ANRR definition of democracy (Single demo and drop reversal-only countries): Panel (c)
***

preserve 
tsset wbcode2 year
xtreg y dem yy* tradewb ginv, fe
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
* if a country has no democratisation or reversal it is either always a demo or never
xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
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
gen ddem = d.dem

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddem dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  

sort wbcode2 year
by wbcode2: egen dem_sum = sum(demevent) if e(sample)
by wbcode2: egen rev_sum = sum(revevent) if e(sample)

* Robustness check: countries with single democratisation (and no reversal-only countries)
replace dem_sample=0 if (dem_sum>1 & !missing(dem_sum)) 
replace dem_sample=0 if (rev_sum>0 & !missing(rev_sum)) 

drop if !e(sample)
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimate identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum dem
tabstat csum dem, by(wbcode)

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG (Column 1)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok (Column 2)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** MG with covariates (Column 3)
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates (Column 4) 
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem)"
esttab model* using "$texfolder/20220421_Fullsample_singledemnorev_DID_ANRR.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
restore 


***
*** ANRR definition of democracy (single demo): Panel (d)
***

preserve 
tsset wbcode2 year
xtreg y dem yy* tradewb ginv, fe
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
* if a country has no democratisation or reversal it is either always a demo or never
xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
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
gen ddem = d.dem

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddem dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  

sort wbcode2 year
by wbcode2: egen dem_sum = sum(demevent) if e(sample)
by wbcode2: egen rev_sum = sum(revevent) if e(sample)

* Robustness check: countries with single democratisation
replace dem_sample=0 if (dem_sum>1 & !missing(dem_sum)) 
replace dem_sample=0 if (rev_sum!=0 & dem_sum==0)

drop if !e(sample)
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimate identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum dem
tabstat csum dem, by(wbcode)

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG (Column 1)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok (Column 2)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** MG with covariates (Column 3)
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates (Column 4) 
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem)"
esttab model* using "$texfolder/20220421_Fullsample_singledem_DID_ANRR.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
restore 


***
*** ANRR definition of democracy (multiple demo): Panel (e)
***

preserve 
tsset wbcode2 year
xtreg y dem yy* tradewb ginv, fe
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
* if a country has no democratisation or reversal it is either always a demo or never
xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
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
gen ddem = d.dem

* p-bar_x and p-bar_y
sort wbcode2 year 
tsset wbcode2 year
local z "yNT ginvNT tradewbNT"
foreach k of local z{
	forvalues l=0/5{
		gen l`l'`k'=l`l'.`k'
	}
}

local z "ddem dtradewb dginv dy"
foreach k of local z{
	sort wbcode2 year
	gen l0`k' = `k'
	gen l1`k' = l1.`k'
	gen l2`k' = l2.`k'
	gen l3`k' = l3.`k'
	gen l4`k' = l4.`k'
}	

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  

sort wbcode2 year
by wbcode2: egen dem_sum = sum(demevent) if e(sample)
by wbcode2: egen rev_sum = sum(revevent) if e(sample)

* Robustness check: countries with single democratisation
replace dem_sample=0 if (dem_sum<2 & !missing(dem_sum)) 
replace dem_sample=0 if (rev_sum!=0 & dem_sum==0)

drop if !e(sample)
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimate identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum dem
tabstat csum dem, by(wbcode)

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG (Column 1)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Plain vanilla Chan and Kwok (Column 2)
eststo: xtmg y dem  ///
		l0ddem l1ddem l2ddem  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** MG with covariates (Column 3)
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

local i = `i' + 1

** Chan and Kwok with covariates (Column 4) 
eststo: xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bANRR = e(betas)
scalar countries1=e(N_g)
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
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var dem "Democracy (ANRR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(dem)"
esttab model* using "$texfolder/20220421_Fullsample_multiple_DID_ANRR.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
restore 

exit 
	