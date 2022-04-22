******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*        Appendix Table E-4 - Sample Reductions Madsen et al (2015) 
******************************************************************************
*                       Created: 8th February 2021
*                       Revised: 21st April 2022
******************************************************************************
*                          Markus Eberhardt
******************************************************************************
*            School of Economics, University of Nottingham,
*            University Park, Nottingham NG7 2RD, England
*              email: markus.eberhardt@nottingham.ac.uk
*            web: https://sites.google.com/site/medevecon/
******************************************************************************

* Note: EER-D-14-00083_data_PANEL_1820_2000.dta is the dataset used in 
*       Madsen et al (2015).

global path "D:/Dropbox/Curious GMM/Acemoglu Replication"
global path "/Users/lezme/Dropbox/Curious GMM/Acemoglu Replication"

global dofolder "$path/do_files"
global otherdata "$path/Data"
global outputfolder "$path/Output"
global rubbishfolder "$path/Rubbish"
global texfolder "$path/Tex"

use "$path/Madsen historical data/EER-D-14-00083_data_PANEL_1820_2000.dta", clear

*** MADSEN ET AL - IV regressions Table 4 column (1)

xi: xtivreg2 lny  l.lny  (l.Dem= l.DemF)    i.time  ,   fe  cl(countryid)
sum l.Dem if e(sample)
local ldemM = r(sd)
nlcom (100*(_b[l.Dem]/(1-_b[l.lny]))*`ldemM')
* 95%

gen c=1 if e(sample)
sort countrycode
by countrycode: egen csum=sum(c)
tab csum if e(sample)

sort csum
by csum: tab countrycode if e(sample)

matrix Madsen_Table4Column1=J(14,9,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset countryid time
local z "1 2 3 4 5 6 7 8 9 12 13 14 15 16"
foreach k of local z{
	quietly{
	xi: xtivreg2 lny  l.lny  (l.Dem= l.DemF)   i.time if csum>=`k' , fe cl(countryid) ffirst
	* Min number of observations
	mat Madsen_Table4Column1[`i',1]=`k'
	* Total number of observations
	mat Madsen_Table4Column1[`i',2]=e(N)
	* Number of Countries
	mat Madsen_Table4Column1[`i',3]=e(N_g)
	* Dem coefficient (short-run)
	mat Madsen_Table4Column1[`i',6]=_b[l.Dem]
	mat Madsen_Table4Column1[`i',7]=_se[l.Dem]	
	sum l.Dem if e(sample)
	local ldemM = r(sd)
	nlcom (100*(_b[l.Dem]/(1-_b[l.lny]))*`ldemM') (_b[l.lny]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat Madsen_Table4Column1[`i',4]=b_2
	mat Madsen_Table4Column1[`i',5]=se_2
	* Dem coefficient (long-run)
	mat Madsen_Table4Column1[`i',8]=b_1
	mat Madsen_Table4Column1[`i',9]=se_1	
	}
	di in gr "Table 4(1), Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more decadal observations. n=" in ye e(N) in gr " obs." 

	local i = `i'+1
}	


mat list Madsen_Table4Column1
svmat Madsen_Table4Column1, names(IV41_)
rename IV41_1 IV41_minobs
rename IV41_2 IV41_obs
rename IV41_3 IV41_N
rename IV41_4 IV41_b_lagy
rename IV41_5 IV41_se_lagy
rename IV41_6 IV41_b_dem
rename IV41_7 IV41_se_dem
rename IV41_8 IV41_b_lr_dem
rename IV41_9 IV41_se_lr_dem

preserve
keep IV41_*
gen model = _n if !missing(IV41_minobs)
drop if _n>14
order model
sort IV41_minobs
save  "$otherdata/Madsen_4_1_reduction.dta", replace
restore

use "$otherdata/Madsen_4_1_reduction.dta", clear
gen IV41_t_lr_dem = IV41_b_lr_dem/IV41_se_lr_dem

twoway ///
	(connected IV41_b_lr_dem IV41_minobs, ///
		lw(.4) msymbol(O) lcolor(emerald) msize(medium) mfcolor(none) mlcolor(black) mlwidth(thin)) ///
	(scatter IV41_b_lr_dem IV41_minobs if abs(IV41_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(emerald) msize(medium) mfcolor(emerald) mlcolor(black) mlwidth(thin)) ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count", size(medsmall)) ///
	ytitle("Long-run Effect of 1 SD Increase in Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small)) ylabel(0(10)100) yscale(range(0 101)) ///
	yline(0, lpattern(dash) lw(.2)) xlabel(0(1)17)  ///
	legend(rows(1) order(- "IV Table 4(1):" 1 2) label(1 "Long-run Estimate for Democracy") ///
	label(2 "Statistically significant at 10% level") textfirst size(medsmall)) 
graph export "$texfolder/Madsen41_evolution.eps", as(eps) replace 


use "$path/Madsen historical data/EER-D-14-00083_data_PANEL_1820_2000.dta", clear

xi: xtivreg2 lny  l.lny  (l.Dem=    l.DemF) l.Lit   i.time   ,   fe  ffirst cl(countryid)
sum l.Dem if e(sample)
local ldemM = r(sd)
nlcom (100*(_b[l.Dem]/(1-_b[l.lny]))*`ldemM')
* 78%


gen c=1 if e(sample)
sort countrycode
by countrycode: egen csum=sum(c)
tab csum if e(sample)

matrix Madsen_Table4Column6=J(14,9,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset countryid time
local z "1 2 3 4 5 6 7 8 9 12 13 14 15 16"
foreach k of local z{
	quietly{
	xi: xtivreg2 lny  l.lny  (l.Dem= l.DemF)  l.Lit  i.time if csum>=`k' , fe cl(countryid) ffirst
	* Min number of observations
	mat Madsen_Table4Column6[`i',1]=`k'
	* Total number of observations
	mat Madsen_Table4Column6[`i',2]=e(N)
	* Number of Countries
	mat Madsen_Table4Column6[`i',3]=e(N_g)
	* Dem coefficient (short-run)
	mat Madsen_Table4Column6[`i',6]=_b[l.Dem]
	mat Madsen_Table4Column6[`i',7]=_se[l.Dem]	
	sum l.Dem if e(sample)
	local ldemM = r(sd)
	nlcom (100*(_b[l.Dem]/(1-_b[l.lny]))*`ldemM') (_b[l.lny]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat Madsen_Table4Column6[`i',4]=b_2
	mat Madsen_Table4Column6[`i',5]=se_2
	* Dem coefficient (long-run)
	mat Madsen_Table4Column6[`i',8]=b_1
	mat Madsen_Table4Column6[`i',9]=se_1	
	}
	di in gr "Table 4(6), Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more decadal observations. n=" in ye e(N) in gr " obs." 

	local i = `i'+1
}	


mat list Madsen_Table4Column6
svmat Madsen_Table4Column6, names(IV46_)
rename IV46_1 IV46_minobs
rename IV46_2 IV46_obs
rename IV46_3 IV46_N
rename IV46_4 IV46_b_lagy
rename IV46_5 IV46_se_lagy
rename IV46_6 IV46_b_dem
rename IV46_7 IV46_se_dem
rename IV46_8 IV46_b_lr_dem
rename IV46_9 IV46_se_lr_dem

preserve
keep IV46_*
gen model = _n if !missing(IV46_minobs)
drop if _n>14
order model
sort IV46_minobs
save  "$otherdata/Madsen_4_6_reduction.dta", replace
restore

use "$otherdata/Madsen_4_6_reduction.dta", clear
gen IV46_t_lr_dem = IV46_b_lr_dem/IV46_se_lr_dem

twoway ///
	(connected IV46_b_lr_dem IV46_minobs, ///
		lw(.4) msymbol(O) lcolor(emerald) msize(medium) mfcolor(none) mlcolor(black) mlwidth(thin)) ///
	(scatter IV46_b_lr_dem IV46_minobs if abs(IV46_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(emerald) msize(medium) mfcolor(emerald) mlcolor(black) mlwidth(thin)) ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count", size(medsmall)) ///
	ytitle("Long-run Effect of 1 SD Increase in Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small)) ylabel(0(10)100) yscale(range(0 101)) ///
	yline(0, lpattern(dash) lw(.2)) xlabel(0(1)17)  ///
	legend(rows(1) order(- "IV Table 4(6):" 1 2) label(1 "Long-run Estimate for Democracy") ///
	label(2 "Statistically significant at 10% level") textfirst size(medsmall)) 
graph export "$texfolder/Madsen46_evolution.eps", as(eps) replace 




use "$path/Madsen historical data/EER-D-14-00083_data_PANEL_1820_2000.dta", clear

xi: xtivreg2 lny l.lny (Dem= DemF) i.time   ,   fe  ffirst cl(countryid)
sum Dem if e(sample)
local ldemM = r(sd)
nlcom (100*(_b[Dem]/(1-_b[l.lny]))*`ldemM') 
* 122%

gen c=1 if e(sample)
sort countrycode
by countrycode: egen csum=sum(c)
tab csum if e(sample)

matrix Madsen_Table5Column9=J(14,9,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset countryid time
local z "1 2 3 4 5 6 7 8 9 12 13 14 15 16"
foreach k of local z{
	quietly{
	xi: xtivreg2 lny l.lny (Dem= DemF)  i.time if csum>=`k' , fe cl(countryid) ffirst
	* Min number of observations
	mat Madsen_Table5Column9[`i',1]=`k'
	* Total number of observations
	mat Madsen_Table5Column9[`i',2]=e(N)
	* Number of Countries
	mat Madsen_Table5Column9[`i',3]=e(N_g)
	* Dem coefficient (short-run)
	mat Madsen_Table5Column9[`i',6]=_b[Dem]
	mat Madsen_Table5Column9[`i',7]=_se[Dem]	
	sum Dem if e(sample)
	local demM = r(sd)
	nlcom (100*(_b[Dem]/(1-_b[l.lny]))*`ldemM') (_b[l.lny]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat Madsen_Table5Column9[`i',4]=b_2
	mat Madsen_Table5Column9[`i',5]=se_2
	* Dem coefficient (long-run)
	mat Madsen_Table5Column9[`i',8]=b_1
	mat Madsen_Table5Column9[`i',9]=se_1	
	}
	di in gr "Table 5(9), Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more decadal observations. n=" in ye e(N) in gr " obs." 

	local i = `i'+1
}	


mat list Madsen_Table5Column9
svmat Madsen_Table5Column9, names(IV59_)
rename IV59_1 IV59_minobs
rename IV59_2 IV59_obs
rename IV59_3 IV59_N
rename IV59_4 IV59_b_lagy
rename IV59_5 IV59_se_lagy
rename IV59_6 IV59_b_dem
rename IV59_7 IV59_se_dem
rename IV59_8 IV59_b_lr_dem
rename IV59_9 IV59_se_lr_dem

preserve
keep IV59_*
gen model = _n if !missing(IV59_minobs)
drop if _n>14
order model
sort IV59_minobs
save  "$otherdata/Madsen_5_9_reduction.dta", replace
restore

use "$otherdata/Madsen_5_9_reduction.dta", clear
gen IV59_t_lr_dem = IV59_b_lr_dem/IV59_se_lr_dem

twoway ///
	(connected IV59_b_lr_dem IV59_minobs, ///
		lw(.4) msymbol(O) lcolor(emerald) msize(medium) mfcolor(none) mlcolor(black) mlwidth(thin)) ///
	(scatter IV59_b_lr_dem IV59_minobs if abs(IV59_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(emerald) msize(medium) mfcolor(emerald) mlcolor(black) mlwidth(thin)) ///
, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count", size(medsmall)) ///
	ytitle("Long-run Effect of 1 SD Increase in Democracy", size(medsmall)) bgcolor(white) ///
	ylabel(, labsize(small) glcolor(gs14) glwidth(medthin) angle(0)) ///
	xlabel(, labsize(small)) ylabel(0(10)100) yscale(range(0 101)) ///
	yline(0, lpattern(dash) lw(.2)) xlabel(0(1)17)  ///
	legend(rows(1) order(- "IV Table 5(9):" 1 2) label(1 "Long-run Estimate for Democracy") ///
	label(2 "Statistically significant at 10% level") textfirst size(medsmall)) 
graph export "$texfolder/Madsen59_evolution.eps", as(eps) replace 


exit 

use "$otherdata/Madsen_4_1_reduction.dta", clear
gen IV41_t_lr_dem = IV41_b_lr_dem/IV41_se_lr_dem
merge model using "$otherdata/Madsen_4_6_reduction.dta"  ///
"$otherdata/Madsen_5_9_reduction.dta", sort
gen IV46_t_lr_dem = IV46_b_lr_dem/IV46_se_lr_dem
gen IV59_t_lr_dem = IV59_b_lr_dem/IV59_se_lr_dem


twoway ///
	(connected IV41_b_lr_dem IV41_minobs, ///
		lw(.6) msymbol(O) lcolor(navy) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV41_b_lr_dem IV41_minobs if abs(IV41_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(navy) msize(medlarge) mfcolor(navy) mlcolor(black) mlwidth(thin)) ///
	(connected IV46_b_lr_dem IV46_minobs, ///
		lw(.4) msymbol(O) lcolor(orange) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV46_b_lr_dem IV46_minobs if abs(IV46_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(orange) msize(medlarge) mfcolor(orange) mlcolor(black) mlwidth(thin)) ///		
	(connected IV59_b_lr_dem IV59_minobs, ///
		lw(.4) msymbol(O) lcolor(teal) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV59_b_lr_dem IV59_minobs if abs(IV59_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(teal) msize(medlarge) mfcolor(teal) mlcolor(black) mlwidth(thin)) ///		
	(line IV41_b_lr_dem IV41_minobs, ///
		 lcolor(none) yaxis(2)) ///
	(line IV46_b_lr_dem IV46_minobs, ///
		 lcolor(none) yaxis(2)) ///
	(line IV59_b_lr_dem IV59_minobs, ///
		 lcolor(none) yaxis(2)) ///
	, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count (in decades)", size(medsmall)) ///
	ytitle("LR Effect of 1 SD Increase in Democracy", size(medsmall)) bgcolor(white) ///
	ytitle("LR Effect of 1 SD Increase in Democracy", size(medsmall) axis(2))  ///
	ylabel(-80(40)120, grid glcolor(gs6) glwidth(medthin) glpattern(shortdash) labsize(small) angle(0) axis(2)) ///
	ylabel(-80(40)120, nogrid glcolor(gs6) glwidth(medthin) labsize(small) angle(0)) ///
	xlabel(, labsize(small))  yscale(range(-80 121) axis(2)) yscale(range(-80 121)) ///
	yline(0, lpattern(dash) lw(.2)) xlabel(0(1)17)  ///
	legend(rows(2) order(- "Long-run Estimate on Democracy:" 1 3 5 - "Statistical significant at 10% level:" 2 4 6) ///
	label(1 "IV Table 4(1)") ///
	label(2 "") ///
	label(3 "IV Table 4(6)") ///
	label(4 "") ///
	label(5 "IV Table 5(9)") ///
	label(6 "") ///
	size(medsmall))
graph export "$texfolder/MadsenIV_evolution_three.eps", as(eps) replace 
	
	

twoway ///
	(connected IV41_b_lr_dem IV41_minobs, ///
		lw(.8) msymbol(O) lcolor(pink*1.5) msize(large) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV41_b_lr_dem IV41_minobs if abs(IV41_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*1.5) msize(large) mfcolor(pink*1.5) mlcolor(black) mlwidth(thin)) ///
	(connected IV46_b_lr_dem IV46_minobs, lpattern(dash) ///
		lw(.8) msymbol(O) lcolor(pink*1) msize(large) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV46_b_lr_dem IV46_minobs if abs(IV46_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*1) msize(large) mfcolor(pink*1) mlcolor(black) mlwidth(thin)) ///		
	(connected IV59_b_lr_dem IV59_minobs, lpattern(shortdash) ///
		lw(.8) msymbol(O) lcolor(pink*.5) msize(large) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV59_b_lr_dem IV59_minobs if abs(IV59_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*.5) msize(large) mfcolor(pink*.5) mlcolor(black) mlwidth(thin)) ///		
	(line IV41_b_lr_dem IV41_minobs, ///
		 lcolor(none) yaxis(2)) ///
	(line IV46_b_lr_dem IV46_minobs, ///
		 lcolor(none) yaxis(2)) ///
	(line IV59_b_lr_dem IV59_minobs, ///
		 lcolor(none) yaxis(2)) ///
	, xsize(9) ysize(4) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count (in decades)", size(medsmall)) ///
	ytitle("LR Effect of 1 SD Increase in Democracy", size(medsmall)) bgcolor(white) ///
	ytitle("LR Effect of 1 SD Increase in Democracy", size(medsmall) axis(2))  ///
	ylabel(-80(40)120, grid glcolor(gs6) glwidth(medthin) glpattern(shortdash) labsize(small) angle(0) axis(2)) ///
	ylabel(-80(40)120, nogrid glcolor(gs6) glwidth(medthin) labsize(small) angle(0)) ///
	xlabel(, labsize(small))  yscale(range(-80 121) axis(2)) yscale(range(-80 121)) ///
	yline(0, lpattern(dash) lw(.2)) xlabel(0(1)17)  ///
	legend(rows(2) order(- "Long-run Estimate on Democracy:" 1 3 5 - "Statistical significant at 10% level:" 2 4 6) ///
	label(1 "IV Table 4(1)") ///
	label(2 "") ///
	label(3 "IV Table 4(6)") ///
	label(4 "") ///
	label(5 "IV Table 5(9)") ///
	label(6 "") ///
	size(medsmall))
graph export "$texfolder/MadsenIV_evolution_three_pink.eps", as(eps) replace 	


****
**** Table
****

use "$otherdata/Madsen_4_1_reduction.dta", clear
gen IV41_t_lr_dem = IV41_b_lr_dem/IV41_se_lr_dem

use "$otherdata/Madsen_4_6_reduction.dta", clear
gen IV46_t_lr_dem = IV46_b_lr_dem/IV46_se_lr_dem

use "$otherdata/Madsen_5_9_reduction.dta", clear
gen IV59_t_lr_dem = IV59_b_lr_dem/IV59_se_lr_dem
