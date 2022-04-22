******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                       Table 3 - Quartile TE 
******************************************************************************
*                       Created: 2nd November 2021
*                       Revised: 21st April 2022
******************************************************************************
*                          Markus Eberhardt
******************************************************************************
*            School of Economics, University of Nottingham,
*            University Park, Nottingham NG7 2RD, England
*              email: markus.eberhardt@nottingham.ac.uk
*            web: https://sites.google.com/site/medevecon/
******************************************************************************

* Note: Results are provided as .tex files separately for quartile within 
* 		each panel (a-d). The Robust mean estaimates are again in separate tex 
*       files for each panel.

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


***
*** ANRR definition of democracy - Panel (a)
***
use "$path/DDCGdata_final.dta", clear
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

** Plain vanilla MG
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
mat betas = e(betas)
preserve
	keep if dem_sample==1 & c==1
	* Cross the threshold
	tsset wbcode2 year
	gen count_dem = 1 if ((dem==0 & l.dem==1) | (dem==1 & l.dem==0))
	by wbcode2: egen dem_count = sum(count_dem)
	replace count_dem = 0 if (dem==0 & l.dem==1)
	* Cross the threshold from below
	by wbcode2: egen trudem_dum_count = sum(count_dem)
	by wbcode2: egen dem_exposure = sum(dem)
	by wbcode2: egen sample_obs = sum(c)
	* Countries
	tabstat wbcode2 if dem_sample==1 & c==1, by(wbcode) save
	tabstatmat wbcode2stat, nototal
	* Years
	tabstat year if dem_sample==1 & c==1, by(wbcode) stats(min max) save
	tabstatmat yearstat, nototal
	tabstat dem_count trudem_dum_count dem_exposure sample_obs, by(wbcode) stats(mean) save
	tabstatmat exposure, nototal
	drop _all
	svmat wbcode2stat, name(wbcode)
	rename wbcode1 wbcode2
	svmat yearstat, name(year)
	rename year1 startyear
	rename year2 endyear
	svmat exposure, name(count)
	rename count1 flips
	rename count2 democratisations
	rename count3 exposure
	rename count4 total_obs
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	gen wbcode=""
	replace wbcode="ALB" if wbcode2==4
	replace wbcode="ARG" if wbcode2==6
	replace wbcode="BDI" if wbcode2==13
	replace wbcode="BEN" if wbcode2==15
	replace wbcode="BFA" if wbcode2==16
	replace wbcode="BGD" if wbcode2==17
	replace wbcode="BGR" if wbcode2==18
	replace wbcode="BOL" if wbcode2==25
	replace wbcode="BRA" if wbcode2==26
	replace wbcode="BTN" if wbcode2==29
	replace wbcode="CAF" if wbcode2==31
	replace wbcode="CHL" if wbcode2==34
	replace wbcode="CIV" if wbcode2==36
	replace wbcode="COG" if wbcode2==38
	replace wbcode="COM" if wbcode2==40
	replace wbcode="CPV" if wbcode2==41
	replace wbcode="DOM" if wbcode2==50
	replace wbcode="ECU" if wbcode2==52
	replace wbcode="ESP" if wbcode2==55
	replace wbcode="ETH" if wbcode2==57
	replace wbcode="FJI" if wbcode2==59
	replace wbcode="GHA" if wbcode2==65
	replace wbcode="GMB" if wbcode2==69
	replace wbcode="GNB" if wbcode2==70
	replace wbcode="GRC" if wbcode2==72
	replace wbcode="GRD" if wbcode2==73
	replace wbcode="GTM" if wbcode2==75
	replace wbcode="GUY" if wbcode2==78
	replace wbcode="HND" if wbcode2==80
	replace wbcode="HUN" if wbcode2==83
	replace wbcode="IDN" if wbcode2==84
	replace wbcode="KEN" if wbcode2==96
	replace wbcode="KOR" if wbcode2==101
	replace wbcode="LSO" if wbcode2==109
	replace wbcode="MDG" if wbcode2==116
	replace wbcode="MEX" if wbcode2==118
	replace wbcode="MLI" if wbcode2==120
	replace wbcode="MNG" if wbcode2==123
	replace wbcode="MOZ" if wbcode2==124
	replace wbcode="MRT" if wbcode2==125
	replace wbcode="MWI" if wbcode2==128
	replace wbcode="NER" if wbcode2==133
	replace wbcode="NIC" if wbcode2==135
	replace wbcode="NPL" if wbcode2==138
	replace wbcode="PAK" if wbcode2==141
	replace wbcode="PAN" if wbcode2==142
	replace wbcode="PER" if wbcode2==143
	replace wbcode="PHL" if wbcode2==144
	replace wbcode="PRT" if wbcode2==149
	replace wbcode="SDN" if wbcode2==158
	replace wbcode="SEN" if wbcode2==159
	replace wbcode="SLE" if wbcode2==163
	replace wbcode="SUR" if wbcode2==167
	replace wbcode="THA" if wbcode2==177
	replace wbcode="TUR" if wbcode2==183
	replace wbcode="UGA" if wbcode2==186
	replace wbcode="URY" if wbcode2==188
	replace wbcode="VEN" if wbcode2==192
	replace wbcode="ZAF" if wbcode2==199
	replace wbcode="ZMB" if wbcode2==201
	replace wbcode="ZWE" if wbcode2==202
	order wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	keep wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	save "$outputfolder/PlainVanilla_MG_country.dta", replace
restore
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

** Plain vanilla Chan and Kwok
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
mat betas = e(betas)
preserve
	keep if dem_sample==1 & c==1
	* Cross the threshold
	tsset wbcode2 year
	gen count_dem = 1 if ((dem==0 & l.dem==1) | (dem==1 & l.dem==0))
	by wbcode2: egen dem_count = sum(count_dem)
	replace count_dem = 0 if (dem==0 & l.dem==1)
	* Cross the threshold from below
	by wbcode2: egen trudem_dum_count = sum(count_dem)
	by wbcode2: egen dem_exposure = sum(dem)
	by wbcode2: egen sample_obs = sum(c)
	* Countries
	tabstat wbcode2 if dem_sample==1 & c==1, by(wbcode) save
	tabstatmat wbcode2stat, nototal
	* Years
	tabstat year if dem_sample==1 & c==1, by(wbcode) stats(min max) save
	tabstatmat yearstat, nototal
	tabstat dem_count trudem_dum_count dem_exposure sample_obs, by(wbcode) stats(mean) save
	tabstatmat exposure, nototal
	drop _all
	svmat wbcode2stat, name(wbcode)
	rename wbcode1 wbcode2
	svmat yearstat, name(year)
	rename year1 startyear
	rename year2 endyear
	svmat exposure, name(count)
	rename count1 flips
	rename count2 democratisations
	rename count3 exposure
	rename count4 total_obs
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	gen wbcode=""
	replace wbcode="ALB" if wbcode2==4
	replace wbcode="ARG" if wbcode2==6
	replace wbcode="BDI" if wbcode2==13
	replace wbcode="BEN" if wbcode2==15
	replace wbcode="BFA" if wbcode2==16
	replace wbcode="BGD" if wbcode2==17
	replace wbcode="BGR" if wbcode2==18
	replace wbcode="BOL" if wbcode2==25
	replace wbcode="BRA" if wbcode2==26
	replace wbcode="BTN" if wbcode2==29
	replace wbcode="CAF" if wbcode2==31
	replace wbcode="CHL" if wbcode2==34
	replace wbcode="CIV" if wbcode2==36
	replace wbcode="COG" if wbcode2==38
	replace wbcode="COM" if wbcode2==40
	replace wbcode="CPV" if wbcode2==41
	replace wbcode="DOM" if wbcode2==50
	replace wbcode="ECU" if wbcode2==52
	replace wbcode="ESP" if wbcode2==55
	replace wbcode="ETH" if wbcode2==57
	replace wbcode="FJI" if wbcode2==59
	replace wbcode="GHA" if wbcode2==65
	replace wbcode="GMB" if wbcode2==69
	replace wbcode="GNB" if wbcode2==70
	replace wbcode="GRC" if wbcode2==72
	replace wbcode="GRD" if wbcode2==73
	replace wbcode="GTM" if wbcode2==75
	replace wbcode="GUY" if wbcode2==78
	replace wbcode="HND" if wbcode2==80
	replace wbcode="HUN" if wbcode2==83
	replace wbcode="IDN" if wbcode2==84
	replace wbcode="KEN" if wbcode2==96
	replace wbcode="KOR" if wbcode2==101
	replace wbcode="LSO" if wbcode2==109
	replace wbcode="MDG" if wbcode2==116
	replace wbcode="MEX" if wbcode2==118
	replace wbcode="MLI" if wbcode2==120
	replace wbcode="MNG" if wbcode2==123
	replace wbcode="MOZ" if wbcode2==124
	replace wbcode="MRT" if wbcode2==125
	replace wbcode="MWI" if wbcode2==128
	replace wbcode="NER" if wbcode2==133
	replace wbcode="NIC" if wbcode2==135
	replace wbcode="NPL" if wbcode2==138
	replace wbcode="PAK" if wbcode2==141
	replace wbcode="PAN" if wbcode2==142
	replace wbcode="PER" if wbcode2==143
	replace wbcode="PHL" if wbcode2==144
	replace wbcode="PRT" if wbcode2==149
	replace wbcode="SDN" if wbcode2==158
	replace wbcode="SEN" if wbcode2==159
	replace wbcode="SLE" if wbcode2==163
	replace wbcode="SUR" if wbcode2==167
	replace wbcode="THA" if wbcode2==177
	replace wbcode="TUR" if wbcode2==183
	replace wbcode="UGA" if wbcode2==186
	replace wbcode="URY" if wbcode2==188
	replace wbcode="VEN" if wbcode2==192
	replace wbcode="ZAF" if wbcode2==199
	replace wbcode="ZMB" if wbcode2==201
	replace wbcode="ZWE" if wbcode2==202
	order wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	keep wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	save "$outputfolder/PlainVanilla_CKMG_country.dta", replace
restore
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

** MG with covariates
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
mat betas = e(betas)
preserve
	keep if dem_sample==1 & c==1
	* Cross the threshold
	tsset wbcode2 year
	gen count_dem = 1 if ((dem==0 & l.dem==1) | (dem==1 & l.dem==0))
	by wbcode2: egen dem_count = sum(count_dem)
	replace count_dem = 0 if (dem==0 & l.dem==1)
	* Cross the threshold from below
	by wbcode2: egen trudem_dum_count = sum(count_dem)
	by wbcode2: egen dem_exposure = sum(dem)
	by wbcode2: egen sample_obs = sum(c)
	* Countries
	tabstat wbcode2 if dem_sample==1 & c==1, by(wbcode) save
	tabstatmat wbcode2stat, nototal
	* Years
	tabstat year if dem_sample==1 & c==1, by(wbcode) stats(min max) save
	tabstatmat yearstat, nototal
	tabstat dem_count trudem_dum_count dem_exposure sample_obs, by(wbcode) stats(mean) save
	tabstatmat exposure, nototal
	drop _all
	svmat wbcode2stat, name(wbcode)
	rename wbcode1 wbcode2
	svmat yearstat, name(year)
	rename year1 startyear
	rename year2 endyear
	svmat exposure, name(count)
	rename count1 flips
	rename count2 democratisations
	rename count3 exposure
	rename count4 total_obs
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	gen wbcode=""
	replace wbcode="ALB" if wbcode2==4
	replace wbcode="ARG" if wbcode2==6
	replace wbcode="BDI" if wbcode2==13
	replace wbcode="BEN" if wbcode2==15
	replace wbcode="BFA" if wbcode2==16
	replace wbcode="BGD" if wbcode2==17
	replace wbcode="BGR" if wbcode2==18
	replace wbcode="BOL" if wbcode2==25
	replace wbcode="BRA" if wbcode2==26
	replace wbcode="BTN" if wbcode2==29
	replace wbcode="CAF" if wbcode2==31
	replace wbcode="CHL" if wbcode2==34
	replace wbcode="CIV" if wbcode2==36
	replace wbcode="COG" if wbcode2==38
	replace wbcode="COM" if wbcode2==40
	replace wbcode="CPV" if wbcode2==41
	replace wbcode="DOM" if wbcode2==50
	replace wbcode="ECU" if wbcode2==52
	replace wbcode="ESP" if wbcode2==55
	replace wbcode="ETH" if wbcode2==57
	replace wbcode="FJI" if wbcode2==59
	replace wbcode="GHA" if wbcode2==65
	replace wbcode="GMB" if wbcode2==69
	replace wbcode="GNB" if wbcode2==70
	replace wbcode="GRC" if wbcode2==72
	replace wbcode="GRD" if wbcode2==73
	replace wbcode="GTM" if wbcode2==75
	replace wbcode="GUY" if wbcode2==78
	replace wbcode="HND" if wbcode2==80
	replace wbcode="HUN" if wbcode2==83
	replace wbcode="IDN" if wbcode2==84
	replace wbcode="KEN" if wbcode2==96
	replace wbcode="KOR" if wbcode2==101
	replace wbcode="LSO" if wbcode2==109
	replace wbcode="MDG" if wbcode2==116
	replace wbcode="MEX" if wbcode2==118
	replace wbcode="MLI" if wbcode2==120
	replace wbcode="MNG" if wbcode2==123
	replace wbcode="MOZ" if wbcode2==124
	replace wbcode="MRT" if wbcode2==125
	replace wbcode="MWI" if wbcode2==128
	replace wbcode="NER" if wbcode2==133
	replace wbcode="NIC" if wbcode2==135
	replace wbcode="NPL" if wbcode2==138
	replace wbcode="PAK" if wbcode2==141
	replace wbcode="PAN" if wbcode2==142
	replace wbcode="PER" if wbcode2==143
	replace wbcode="PHL" if wbcode2==144
	replace wbcode="PRT" if wbcode2==149
	replace wbcode="SDN" if wbcode2==158
	replace wbcode="SEN" if wbcode2==159
	replace wbcode="SLE" if wbcode2==163
	replace wbcode="SUR" if wbcode2==167
	replace wbcode="THA" if wbcode2==177
	replace wbcode="TUR" if wbcode2==183
	replace wbcode="UGA" if wbcode2==186
	replace wbcode="URY" if wbcode2==188
	replace wbcode="VEN" if wbcode2==192
	replace wbcode="ZAF" if wbcode2==199
	replace wbcode="ZMB" if wbcode2==201
	replace wbcode="ZWE" if wbcode2==202
	order wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	keep wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	save "$outputfolder/FullCov_MG_country.dta", replace
restore
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

** Chan and Kwok with covariates 
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
mat betas = e(betas)
preserve
	keep if dem_sample==1 & c==1
	* Cross the threshold
	tsset wbcode2 year
	gen count_dem = 1 if ((dem==0 & l.dem==1) | (dem==1 & l.dem==0))
	by wbcode2: egen dem_count = sum(count_dem)
	replace count_dem = 0 if (dem==0 & l.dem==1)
	* Cross the threshold from below
	by wbcode2: egen trudem_dum_count = sum(count_dem)
	by wbcode2: egen dem_exposure = sum(dem)
	by wbcode2: egen sample_obs = sum(c)
	* Countries
	tabstat wbcode2 if dem_sample==1 & c==1, by(wbcode) save
	tabstatmat wbcode2stat, nototal
	* Years
	tabstat year if dem_sample==1 & c==1, by(wbcode) stats(min max) save
	tabstatmat yearstat, nototal
	tabstat dem_count trudem_dum_count dem_exposure sample_obs, by(wbcode) stats(mean) save
	tabstatmat exposure, nototal
	drop _all
	svmat wbcode2stat, name(wbcode)
	rename wbcode1 wbcode2
	svmat yearstat, name(year)
	rename year1 startyear
	rename year2 endyear
	svmat exposure, name(count)
	rename count1 flips
	rename count2 democratisations
	rename count3 exposure
	rename count4 total_obs
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	gen wbcode=""
	replace wbcode="ALB" if wbcode2==4
	replace wbcode="ARG" if wbcode2==6
	replace wbcode="BDI" if wbcode2==13
	replace wbcode="BEN" if wbcode2==15
	replace wbcode="BFA" if wbcode2==16
	replace wbcode="BGD" if wbcode2==17
	replace wbcode="BGR" if wbcode2==18
	replace wbcode="BOL" if wbcode2==25
	replace wbcode="BRA" if wbcode2==26
	replace wbcode="BTN" if wbcode2==29
	replace wbcode="CAF" if wbcode2==31
	replace wbcode="CHL" if wbcode2==34
	replace wbcode="CIV" if wbcode2==36
	replace wbcode="COG" if wbcode2==38
	replace wbcode="COM" if wbcode2==40
	replace wbcode="CPV" if wbcode2==41
	replace wbcode="DOM" if wbcode2==50
	replace wbcode="ECU" if wbcode2==52
	replace wbcode="ESP" if wbcode2==55
	replace wbcode="ETH" if wbcode2==57
	replace wbcode="FJI" if wbcode2==59
	replace wbcode="GHA" if wbcode2==65
	replace wbcode="GMB" if wbcode2==69
	replace wbcode="GNB" if wbcode2==70
	replace wbcode="GRC" if wbcode2==72
	replace wbcode="GRD" if wbcode2==73
	replace wbcode="GTM" if wbcode2==75
	replace wbcode="GUY" if wbcode2==78
	replace wbcode="HND" if wbcode2==80
	replace wbcode="HUN" if wbcode2==83
	replace wbcode="IDN" if wbcode2==84
	replace wbcode="KEN" if wbcode2==96
	replace wbcode="KOR" if wbcode2==101
	replace wbcode="LSO" if wbcode2==109
	replace wbcode="MDG" if wbcode2==116
	replace wbcode="MEX" if wbcode2==118
	replace wbcode="MLI" if wbcode2==120
	replace wbcode="MNG" if wbcode2==123
	replace wbcode="MOZ" if wbcode2==124
	replace wbcode="MRT" if wbcode2==125
	replace wbcode="MWI" if wbcode2==128
	replace wbcode="NER" if wbcode2==133
	replace wbcode="NIC" if wbcode2==135
	replace wbcode="NPL" if wbcode2==138
	replace wbcode="PAK" if wbcode2==141
	replace wbcode="PAN" if wbcode2==142
	replace wbcode="PER" if wbcode2==143
	replace wbcode="PHL" if wbcode2==144
	replace wbcode="PRT" if wbcode2==149
	replace wbcode="SDN" if wbcode2==158
	replace wbcode="SEN" if wbcode2==159
	replace wbcode="SLE" if wbcode2==163
	replace wbcode="SUR" if wbcode2==167
	replace wbcode="THA" if wbcode2==177
	replace wbcode="TUR" if wbcode2==183
	replace wbcode="UGA" if wbcode2==186
	replace wbcode="URY" if wbcode2==188
	replace wbcode="VEN" if wbcode2==192
	replace wbcode="ZAF" if wbcode2==199
	replace wbcode="ZMB" if wbcode2==201
	replace wbcode="ZWE" if wbcode2==202
	order wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs 
	keep wbcode2 wbcode startyear endyear beta1 flips democratisations exposure total_obs
	save "$outputfolder/FullCov_CKMG_country.dta", replace
restore
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
esttab model* using "$texfolder/20220421_Fullsample_DID_ANRR_qte.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 
 

***
*** Boix et al
***
use "$path/DDCGdata_final.dta", clear 
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

tab wbcode2 demevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 demBMR if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG
eststo: xtmg y demBMR  ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
if dem_sample==1 & c==1, robust
scalar countries1=e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_MG_country_BMR.dta", replace
restore
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

** Plain vanilla Chan and Kwok
eststo: xtmg y demBMR  ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar countries1 = e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_CKMG_country_BMR.dta", replace
restore
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

** MG with covariates
eststo: xtmg y demBMR ginv tradewb ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar countries1 = e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_MG_country_BMR.dta", replace
restore
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

** Chan and Kwok with covariates 
eststo: xtmg y demBMR ginv tradewb ///
		l0ddemBMR l1ddemBMR l2ddemBMR  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bBMR = e(betas)
scalar countries1 = e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_CKMG_country_BMR.dta", replace
restore
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1

est store model`i'

label var demBMR "Democracy (BMR)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(demBMR)"
esttab model* using "$texfolder/20220421_Fullsample_DID_BMR_qte.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 

	

***
*** Cheibub et al - panel (c)
***
use "$path/DDCGdata_final.dta", clear 
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

tab wbcode2 demevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 demCGV if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
	
*** Full sample results

estimates clear 

local i = 1

** Plain vanilla MG
eststo: xtmg y demCGV  ///
		l0ddemCGV l1ddemCGV l2ddemCGV  ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_MG_country_CGV.dta", replace
restore
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

** Plain vanilla Chan and Kwok
eststo: xtmg y demCGV  ///
		l0ddemCGV l1ddemCGV l2ddemCGV  ///
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
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_CKMG_country_CGV.dta", replace
restore
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

** MG with covariates
eststo: xtmg y demCGV ginv tradewb ///
		l0ddemCGV l1ddemCGV l2ddemCGV  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar obs1=e(N)
scalar countries1=e(N_g)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_MG_country_CGV.dta", replace
restore
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

** Chan and Kwok with covariates 
eststo: xtmg y demCGV ginv tradewb ///
		l0ddemCGV l1ddemCGV l2ddemCGV  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bCGV = e(betas)
scalar obs1=e(N)
scalar countries1=e(N_g)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demCGV==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_CKMG_country_CGV.dta", replace
restore
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var demCGV "Democracy (CGV)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(demCGV)"
esttab model* using "$texfolder/20220421_Fullsample_DID_CGV_qte.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace  


***
*** P&S
***
use "$path/DDCGdata_final.dta", clear
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
tab csum
run "$xtmgpath/xtmg.ado"

tab wbcode2 demevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent2 if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 demPS if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
	
*** Full sample results
estimates clear 

local i = 1

** Plain vanilla MG
eststo: xtmg y demPS  ///
		l0ddemPS l1ddemPS l2ddemPS  ///
if dem_sample==1 & c==1, robust
scalar countries1=e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_MG_country_PS.dta", replace
restore
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

** Plain vanilla Chan and Kwok
eststo: xtmg y demPS  ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT  ///
if dem_sample==1 & c==1, robust
scalar countries1=e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/PlainVanilla_CKMG_country_PS.dta", replace
restore
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

** MG with covariates
eststo: xtmg y demPS ginv tradewb ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
if dem_sample==1 & c==1, robust
scalar countries1=e(N_g)
scalar obs1=e(N)
scalar MSE1 = e(sigma)
scalar RMSE1 = sqrt(MSE1)
scalar paramet1 = e(df_m)
scalar minT1 = e(g_min)
mat betas = e(betas)
count if demevent2==1 & e(sample)
scalar Events = r(N)
count if revevent2==1 & e(sample)
scalar Reversals = r(N)
count if demPS==1 & e(sample)
scalar Indem = r(N)
mat betas = e(betas)
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_MG_country_PS.dta", replace
restore
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

** Chan and Kwok with covariates 
eststo: xtmg y demPS ginv tradewb ///
		l0ddemPS l1ddemPS l2ddemPS  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust
mat bPS = e(betas)
scalar countries1=e(N_g)
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
preserve
	clear 
	svmat betas, name(beta)
	replace beta1=. if beta1==0
	save "$outputfolder/FullCov_CKMG_country_PS.dta", replace
restore
estadd scalar countries=countries1
estadd scalar obs=obs1
estadd scalar RMSE = RMSE1
estadd scalar paramet = paramet1
estadd scalar minT = minT1
estadd scalar events = Events
estadd scalar reversals = Reversals
estadd scalar indem = Indem/countries1
est store model`i'

label var demPS "Democracy (PS)"
label var y "GDP per capita"

local esttabformat  "lines label nogaps numbers star(* 0.10 ** 0.05 *** 0.01) staraux mlabels(, depvars) note("") nolegend"
local esttabstats "stats(obs countries events reversals indem RMSE paramet minT, fmt(0 0 0 0 1 3 0 0) label("Observations" "Countries" "Democratisations" "Reversals" "Avg Years in Dem" "RMSE" "Parameters estimated" "Minimum T"))"
local esttabkeep "keep(demPS)"
esttab model* using "$texfolder/20211102_Fullsample_DID_PS_qte.tex", b(3) se(3) abs `esttabkeep' `esttabformat' `esttabstats'  style(tex) replace 

exit 


***
*** Compute Quartile TE
***

***ANRR

local z "PlainVanilla_MG_country PlainVanilla_CKMG_country FullCov_MG_country FullCov_CKMG_country"
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.25)
	est store model`i'
	local i = `i' + 1
}

label var beta1 "Democracy (ANRR)"
esttab model1 model2 model3 model4 using "$texfolder/20211102_ANRR_25.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.5)
	est store model`i'
	local i = `i' + 1
}
label var beta1 "Democracy (ANRR)"
esttab model1 model2 model3 model4 using "$texfolder/20211102_ANRR_50.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.75)
	est store model`i'
	local i = `i' + 1
}
label var beta1 "Democracy (ANRR)"
esttab model1 model2 model3 model4 using "$texfolder/20211102_ANRR_75.tex", b(3) se(3) abs style(tex) replace 

***BMR 
local z "PlainVanilla_MG_country_BMR PlainVanilla_CKMG_country_BMR FullCov_MG_country_BMR FullCov_CKMG_country_BMR"
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.25)
	est store model`i'
	local i = `i' + 1
}

esttab model1 model2 model3 model4 using "$texfolder/20211102_BMR_25.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.5)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_BMR_50.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.75)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_BMR_75.tex", b(3) se(3) abs style(tex) replace

***CGV
local z "PlainVanilla_MG_country_CGV PlainVanilla_CKMG_country_CGV FullCov_MG_country_CGV FullCov_CKMG_country_CGV"
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.25)
	est store model`i'
	local i = `i' + 1
}

esttab model1 model2 model3 model4 using "$texfolder/20211102_CGV_25.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.5)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_CGV_50.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.75)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_CGV_75.tex", b(3) se(3) abs style(tex) replace 

***PS
local z "PlainVanilla_MG_country_PS PlainVanilla_CKMG_country_PS FullCov_MG_country_PS FullCov_CKMG_country_PS"
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.25)
	est store model`i'
	local i = `i' + 1
}

esttab model1 model2 model3 model4 using "$texfolder/20211102_PS_25.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.5)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_PS_50.tex", b(3) se(3) abs style(tex) replace 

est clear 
local i = 1
foreach k of local z{
	use "$outputfolder/`k'.dta", clear 
	eststo: qreg beta1, quantile(.75)
	est store model`i'
	local i = `i' + 1
}
esttab model1 model2 model3 model4 using "$texfolder/20211102_PS_75.tex", b(3) se(3) abs style(tex) replace 

