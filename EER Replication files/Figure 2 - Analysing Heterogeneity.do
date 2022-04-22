******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                   Figure 2 - Analysing heterogeneity
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

*Notes: treisman_ANRR_sample.dta contains the democracy and GDPpc (levels in logs 
* and growth rate) data from ANRR alongside the information on elite biased 
* democratisation (Albertus and Menaldo, 2018) and democracy by mistake (Treisman, 2020).
* The latter two were compiled by the author from information in AM 2018 and the 
* supportive data file to Treisman (2020); the usual disclaimers apply.
* polity_IV_2016.dta contains the PolityIV dataset. 
* literacy_1970.dta are the 1970 literacy rates from Madsen et al (2015).
* country_helper.dta are the base year values for GDP pc in logs from ANRR.

* Paths for merging these additional datasets (all provided as part of the  
* replication files) will have to be adjusted by the user.


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
sort wbcode year 
save "$path/DDCGdata_final.dta", replace 

***
*** ANRR: estimate l-r democracy TE
***

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
*tab csum demevent

tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)

run "$xtmgpath/xtmg.ado"

** Chan and Kwok with covariates 
xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust

tabstat wbcode2 if e(sample), by(wbcode) save
tabstatmat wbcodestat2, nototal
mat betas = e(betas)
clear 
svmat wbcodestat2, name(nwbcode)
gen wbcode=""
replace wbcode="ALB" if nwbcode==4
replace wbcode="ARG" if nwbcode==6
replace wbcode="BDI" if nwbcode==13
replace wbcode="BEN" if nwbcode==15
replace wbcode="BFA" if nwbcode==16
replace wbcode="BGD" if nwbcode==17
replace wbcode="BGR" if nwbcode==18
replace wbcode="BOL" if nwbcode==25
replace wbcode="BRA" if nwbcode==26
replace wbcode="BTN" if nwbcode==29
replace wbcode="CAF" if nwbcode==31
replace wbcode="CHL" if nwbcode==34
replace wbcode="CIV" if nwbcode==36
replace wbcode="COG" if nwbcode==38
replace wbcode="COM" if nwbcode==40
replace wbcode="CPV" if nwbcode==41
replace wbcode="DOM" if nwbcode==50
replace wbcode="ECU" if nwbcode==52
replace wbcode="ESP" if nwbcode==55
replace wbcode="ETH" if nwbcode==57
replace wbcode="FJI" if nwbcode==59
replace wbcode="GHA" if nwbcode==65
replace wbcode="GMB" if nwbcode==69
replace wbcode="GNB" if nwbcode==70
replace wbcode="GRC" if nwbcode==72
replace wbcode="GRD" if nwbcode==73
replace wbcode="GTM" if nwbcode==75
replace wbcode="GUY" if nwbcode==78
replace wbcode="HND" if nwbcode==80
replace wbcode="HUN" if nwbcode==83
replace wbcode="IDN" if nwbcode==84
replace wbcode="KEN" if nwbcode==96
replace wbcode="KOR" if nwbcode==101
replace wbcode="LSO" if nwbcode==109
replace wbcode="MDG" if nwbcode==116
replace wbcode="MEX" if nwbcode==118
replace wbcode="MLI" if nwbcode==120
replace wbcode="MNG" if nwbcode==123
replace wbcode="MOZ" if nwbcode==124
replace wbcode="MRT" if nwbcode==125
replace wbcode="MWI" if nwbcode==128
replace wbcode="NER" if nwbcode==133
replace wbcode="NIC" if nwbcode==135
replace wbcode="NPL" if nwbcode==138
replace wbcode="PAK" if nwbcode==141
replace wbcode="PAN" if nwbcode==142
replace wbcode="PER" if nwbcode==143
replace wbcode="PHL" if nwbcode==144
replace wbcode="PRT" if nwbcode==149
replace wbcode="SDN" if nwbcode==158
replace wbcode="SEN" if nwbcode==159
replace wbcode="SLE" if nwbcode==163
replace wbcode="SUR" if nwbcode==167
replace wbcode="THA" if nwbcode==177
replace wbcode="TUR" if nwbcode==183
replace wbcode="UGA" if nwbcode==186
replace wbcode="URY" if nwbcode==188
replace wbcode="VEN" if nwbcode==192
replace wbcode="ZAF" if nwbcode==199
replace wbcode="ZMB" if nwbcode==201
replace wbcode="ZWE" if nwbcode==202
svmat betas, name(b1_)
replace b1_1=. if b1_1==0

drop b1_4-b1_22
rename b1_1 beta_dem
rename b1_2 beta_ginv
rename b1_3 beta_trade

rreg beta_dem

sort wbcode
save "$path/CK_ANRR_Dynamic.dta", replace 

exit 




***
*** Elite-Biased democratisation and Democracy by Mistake - Panels (a) and (b)
***

use "$path0/Curious GMM/Acemoglu replication/Data/treisman_ANRR_sample.dta", clear 
collapse (mean) treisman albertus albertus_annul, by(wbcode)
merge wbcode using "$path/CK_ANRR_Dynamic.dta", sort 

*** Treisman plot
gen mistake = 0 if !missing(treisman, beta_dem) 
replace mistake = 1 if treisman>4.5 & !missing(treisman, beta_dem)

reg beta_dem mistake, vce(robust) level(90)
est store ls
qreg beta_dem mistake, vce(robust)  level(90)
est store qr
rreg beta_dem mistake,  level(90)
est store rr

coefplot ///
	(ls, mcolor(navy) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(navy)  ///
		lw(.6) recast(rcap)) transform(mistake = @+7.694868) ci(90))  ///
	(qr, mcolor(midblue*.75) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(midblue*.75)  ///
		lw(.6) recast(rcap)) transform(mistake = @+4.670631)  ci(90) yaxis(2)) ///
	(rr, mcolor(pink) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(pink)   ///
		lw(.6) recast(rcap)) transform(mistake = @+6.131984) ci(90)) ///
, vertical  ///
	graphregion(color(white)) plotregion(lcolor(black) lpattern(solid)) ///
    xlabel(, tposition(inside) labsize(medium)) ysize(5) xsize(7) ymtick(##5) ///
	ylabel(-10(10)40, glcolor(gs5) grid glwidth(.1) glp(shortdash) ticks tposition(inside) labsize(medsmall))  ///
	ylabel(-10(10)40, glcolor(gs5) ticks tposition(inside) labsize(medsmall) axis(2))  ///
	coeflabels(mistake = "By Mistake [N=27]" ///
		_cons = "Chosen [N=27]", notick labgap(2)) ymtick(##5, axis(2)) ///
	legend(row(1) order(- "Estimates (+90% CI):" 2 6 5) label(2 "LS") ///
		label(6 "Median") label(5 "Robust")) ytitle("LR Democracy Coefficient") 	///
		order(_cons treisman5) yscale(range(-10.1 40.1)) ytitle("LR Democracy Coefficient", axis(2)) 	
graph export "$texfolder/treisman5_coefficients.eps", as(eps) replace
* Figure 2, Panel (b)

ttest beta_dem, by(mistake) 

*** Albertus plot
gen elite_biased = 0 if !missing(albertus, beta_dem) 
replace elite_biased = 1 if albertus>.5 & !missing(albertus, beta_dem) & albertus_annul!=1

reg beta_dem elite_biased, vce(robust) level(90)
est store ls
qreg beta_dem elite_biased, vce(robust)  level(90)
est store qr
rreg beta_dem elite_biased,  level(90)
est store rr

coefplot ///
	(ls, mcolor(navy) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(navy)  ///
		lw(.6) recast(rcap)) transform(elite_biased = @+17.01394) ci(90))  ///
	(qr, mcolor(midblue*.75) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(midblue*.75) ///
		lw(.6) recast(rcap)) transform(elite_biased = @+14.52895)  ci(90) yaxis(2)) ///
	(rr, mcolor(pink) msymbol(Oh) mlw(thick) msize(large) ciopts(lcolor(pink)   ///
		lw(.6) recast(rcap)) transform(elite_biased = @+15.60362) ci(90)) ///
, vertical  ///
	graphregion(color(white)) plotregion(lcolor(black) lpattern(solid)) ///
    xlabel(, tposition(inside) labsize(medium)) ysize(5) xsize(7) ymtick(##5) ///
	ylabel(-10(10)30, glcolor(gs5) ticks tposition(inside) labsize(medsmall) axis(2))  ///
	ylabel(-10(10)30, glcolor(gs5) grid glwidth(.1) glp(shortdash) ticks tposition(inside) labsize(medsmall))  ///
	coeflabels(elite_biased = "Elite-Biased [N=18]" ///
		_cons = "Popular [N=25]", notick labgap(2))  ymtick(##5, axis(2)) ///
	legend(row(1) order(- "Estimates (+90% CI):"  2 6 5) label(2 "LS") ///
		label(6 "Median") label(5 "Robust")) ytitle("LR Democracy Coefficient") 	///
		order(_cons elite_biased) yscale(range(1.49 3.51)) ytitle("LR Democracy Coefficient", axis(2)) 	
graph export "$texfolder/albertus_coefficients.eps", as(eps) replace
* Figure 2, Panel (a)

ttest beta_dem, by(elite_biased)


***
*** Geography - Panel (c)
***

use "$path/CK_ANRR_Dynamic.dta", clear 
merge wbcode using "$otherdata/country_helper.dta", sort
drop _merge
tab region, gen(r)
rreg beta_dem r1 r5 

gen lypc_base = log(ypc_base)
gen lypc_median = log(ypc_median)

gen r3_sum = 1 if r1==1
replace r3_sum = 2 if r5==1
replace r3_sum = 3 if r1==0 & r5==0		
		
*tabstat share , statistic(N median) by(r3_sum)

twoway ///
	(fpfitci beta_dem lypc_base if lypc_base>4.6 & lypc_base<8.99, lcolor(midblue) lpattern(dash) lw(1.5) level(90) ///
		alcolor(midblue) alwidth(.2) fcolor(midblue*.25))  ///
	(fpfit beta_dem lypc_base if r1==1, lcolor(navy) lpattern(solid) lw(1.5)) ///
	(fpfit beta_dem lypc_base if r5==1, lcolor(cranberry*.5) lw(1.5) lpattern(dash)) ///
	(fpfit beta_dem lypc_base if r1==0 & r5==0, lcolor(orange*.75) lw(1.5) lpattern(shortdash)) ///
	(lfitci beta_dem lypc_base if lypc_base>4.6 & lypc_base<8.99, lcolor(none) alcolor(none) fcolor(none) level(90) yaxis(2)) ///
	(fpfit beta_dem lypc_base if r1==1, lcolor(none)  yaxis(2)) ///
	(fpfit beta_dem lypc_base if r5==1, lcolor(none)  yaxis(2)) ///
	(fpfit beta_dem lypc_base if r1==0 & r5==0, lcolor(none) yaxis(2)) ///
, xsize(7) ysize(4)  ymtick(##2) xmtick(##5) yline(0, lpattern(dash) lw(.2) lcolor(black)) ///
	ylabel(-10(10)30, labsize(medsmall) grid glc(gs8) glw(.2) glp(shortdash) angle(0)) scheme(s2mono)   ///
	ytitle("Long-run democracy coefficient") xtitle("Per capita GDP in the base year (log scale)") ///
	ylabel(-10(10)30, labsize(medsmall) angle(0) axis(2)) ymtick(##2, axis(2)) ///
	legend(order(2 3 4 5) row(1) label(2 "All countries") label(3 "Africa") ///
		label(4 "LAC") label(5 "Other regions")) graphregion(color(white)) ///
		ytitle("Long-run democracy coefficient", axis(2)) ///
		yscale(range(-11 31)) yscale(range(-11 31) axis(2))
graph export "$texfolder/CK_regions_w.eps", as(eps) replace


***
*** Sample years in democracy - panel (d)
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
tab wbcode2 demevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 revevent if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)
tab wbcode2 dem if csum!=0 & !missing(csum) & c==1 & dem_sample==1 & e(sample)

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
run "$xtmgpath/xtmg.ado"
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
gen cc= 1 if e(sample)
sort wbcode2
by wbcode2: egen ccsum=sum(cc)	
gen helper = 0 if e(sample)
replace helper = 1 if dem==1 & e(sample)
by wbcode2: egen helpersum = sum(helper)

run "$xtmgpath/xtmg.ado"

matrix CKcov_years=J(35,11,.)

local i = 1
cls
tsset wbcode2 year
forvalues k=1/33{
	sort wbcode2 year
	qui: xtmg y dem ginv tradewb ///
			l0ddem l1ddem l2ddem   ///
			l0dginv l1dginv l2dginv   ///
			l0dtradewb l1dtradewb l2dtradewb   ///
			l0yNT ///
			l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
			l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
			if helpersum>=`k' & helpersum<=`k'+10 & csum>=24 & c==1, robust res(e`k')

	* Minimum and maximum Years in democracy
	mat CKcov_years[`i',1]=`k'
	mat CKcov_years[`i',2]=`k'+10
	* Total number of observations
	mat CKcov_years[`i',3]=e(N)
	* Number of Countries
	mat CKcov_years[`i',4]=e(N_g)
	
	* Democracy estimate (rolling window of years)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	mat CKcov_years[`i',5]=b_1
	mat CKcov_years[`i',6]=se_1	
	
	qui: xtmg y dem ginv tradewb ///
			l0ddem l1ddem l2ddem   ///
			l0dginv l1dginv l2dginv   ///
			l0dtradewb l1dtradewb l2dtradewb   ///
			l0yNT ///
			l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
			l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
			if helpersum>=`k' & helpersum<=`k'+10 & csum>=24 & c==1, robust res(e_`k')
	
	* Total number of observations
	mat CKcov_years[`i',8]=e(N)
	* Number of Countries
	mat CKcov_years[`i',9]=e(N_g)
	
	* Democracy estimate (rolling window of years)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	mat CKcov_years[`i',10]=b_1
	mat CKcov_years[`i',11]=se_1	
	qui: drop e`k' e_`k'
	di in gr "Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more years in democracy. N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list CKcov_years
svmat CKcov_years, names(CKcovar_)
rename CKcovar_1 CKcovar_minyrs
rename CKcovar_2 CKcovar_maxyrs
rename CKcovar_3 CKcovar_obs
rename CKcovar_4 CKcovar_N
rename CKcovar_5 CKcovar_b_dem
rename CKcovar_6 CKcovar_se_dem
rename CKcovar_7 CKcovar_cd
rename CKcovar_8 CKcovar_obs_l
rename CKcovar_9 CKcovar_N_l
rename CKcovar_10 CKcovar_b_dem_l
rename CKcovar_11 CKcovar_se_dem_l

preserve
keep CKcovar_*
gen model = _n if !missing(CKcovar_minyrs)
drop if _n>33
order model
sort CKcovar_maxyrs
save  "$otherdata/CK_lr_covar_years_alt2.dta", replace
restore

use  "$otherdata/CK_lr_covar_years_alt2.dta", clear

gen min90 = CKcovar_b_dem - 1.645*CKcovar_se_dem
gen max90 = CKcovar_b_dem + 1.645*CKcovar_se_dem

twoway ///
	(rarea min90 max90 CKcovar_maxyrs, lcolor(midblue) lpattern(dash) fcolor(midblue*.25)) ///
	(line CKcovar_b_dem CKcovar_maxyrs, lcolor(midblue) lpattern(dash) lw(1)) ///
	(function y = 0, range(10 45) lcolor(black) lpattern(dash) lw(.3)) ///
	(line CKcovar_N CKcovar_maxyrs, lcolor(gs3) lw(.7) lpattern(solid) yaxis(2)) ///
	, ytitle("Long-run Average Coefficient on Democracy") scheme(s2mono)  ///
		   legend(order(2 4) size(medsmall) rows(1) label(4 "Sample countries (right axis)") ///
			label(2 "Democracy Effect (+ 90% CI, left axis)")) ///
		   ylabel(-10(10)40, labsize(medsmall) grid glc(gs8) glw(.2) glp(shortdash) angle(0)) xmtick(##5) ///
		   xtitle("Countries with {it:t-10} to {it:t} years spent in democracy")  ///
		   graphregion(color(white))  ///
		   xlabel(10(5)45) xscale(range(10 44)) yscale(range(-10 41)) ///
		   yscale(range(-1 30.5) axis(2)) ylabel(0(5)30, axis(2) labsize(medsmall) angle(0)) ymtick(##5, axis(2)) ///
		   ytitle("Number of Countries in the Sample", axis(2)) ymtick(##5) ///
	xsize(7) ysize(4)
graph export "$texfolder/CK_years_in_democracy.eps", as(eps) replace 


***
*** Democratic Legacy - panel (e)
***

use  "$path0/Facets of Democracy/Data/Historical data/polity_IV_2016.dta", clear
keep if year<1976
tab wbcode democracy
gen c=1 if !missing(polity2)
collapse (sum) demo_count=democracy observed=c (lastnm) since_indep, by(wbcode)
gen dem_history = demo_count/observed
sort wbcode
save  "$path/dem_history.dta", replace

use "$path/CK_ANRR_Dynamic.dta", clear 
merge wbcode using "$path/dem_history.dta", sort 
drop if _merge ==2

twoway ///
	(lfitci beta_dem  demo_count if demo_count<100, lcolor(midblue) lpattern(dash) lw(1.5) level(90)  ///
		alcolor(midblue) alwidth(.2) fcolor(midblue*.25))  ///
	(hist demo_count if demo_count<100, w(1) frequency lcolor(none) fcolor(gs11) yaxis(2) barw(.9)) ///
, xsize(7) ysize(4) ymtick(##2) xmtick(##5) yline(0, lpattern(dash) lw(.2) lcolor(black)) ///
	ylabel(-10(10)30, labsize(medsmall) grid glc(gs8) glw(.2) glp(shortdash) angle(0)) scheme(s2mono)  ///
	ylabel(0(5)20, labsize(medsmall) angle(0) axis(2)) yscale(range(-0.1 20.4) axis(2))  ///
	xtitle("Years of democratic history in 1975") ytitle("Long-run democracy coefficient") ///
	xlabel(0(10)90) ytitle("Frequency count: democratic history", axis(2)) ///
	ymtick(##5, axis(2)) graphregion(color(white)) yscale(range(-10.1 30.1)) ///
	legend(order(3 2) size(medsmall) rows(1) label(3 "Sample countries (right axis)") ///
	label(2 "Linear Regression Line (+90% CI)")) 
graph export "$texfolder/CK_dem_history.eps", as(eps) replace  


***
*** By literacy 1970 - panel (f)	
***
	
use "$path/CK_ANRR_Dynamic.dta", clear 
merge wbcode using "$otherdata/literacy_1970.dta", sort	
drop if _merge==2
drop _merge
	
twoway ///
	(qfitci beta_dem Lit, lcolor(midblue) lpattern(dash) lw(1) level(90)  ///
		alcolor(midblue) alwidth(.2) fcolor(midblue*.25))  ///
	(qfitci beta_dem Lit, color(none) yaxis(2) level(90))  ///
, xsize(7) ysize(4) ymtick(##5) xmtick(##10) yline(0, lpattern(dash) lw(.2) lcolor(black)) ///
	ylabel(-10(10)60, labsize(medsmall) grid glc(gs8) glw(.2) glp(shortdash) angle(0)) ytitle("Long-run democracy coefficient") ///
	ylabel(-10(10)60, labsize(medsmall) axis(2) angle(0)) ytitle("Long-run democracy coefficient", axis(2)) ///
	xtitle("Literacy Rate in 1970 (in %)") legend(order(2 1) label(1 "90% CI") label(2 "Quadratic Regression Line")) ///
	graphregion(color(white)) yscale(range(-10.1 60.1)) yscale(range(-10.1 60.1) axis(2)) ymtick(##5, axis(2)) 
graph export "$texfolder/CK_dem_Lit70_relation_qfit.eps", as(eps) replace 
