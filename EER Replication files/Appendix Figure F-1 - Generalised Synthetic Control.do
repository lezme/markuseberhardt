******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*            Appendix Figure F-1 - Generalized Synthetic Controls 
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

* Notes: The fect command is written by Xu Yiqing and available on his personal
*        website: https://yiqingxu.org/software/ The usual disclaimers apply.

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
global texfolder "$path0/Apps/Overleaf/Democracy Growth Heterogeneity and Robustness"

global xtmgpath = "/Users/lezme/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"
global xtmgpath = "D:/Dropbox/Literature/econometric issues/Nonstationary Panel Metrics/Stata"

use "$path/DDCGdata_final.dta", clear


***
*** ANRR
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
sort wbcode2 year
by wbcode2: egen dem_sum = sum(demevent) if e(sample)
by wbcode2: egen rev_sum = sum(revevent) if e(sample)
*** MAKE THE ADJUSTMENT HERE (unless full multiple treatments)
* Only one democratisation, no reversals
*replace dem_sample=0 if (dem_sum>1 & !missing(dem_sum)) 
*replace dem_sample=0 if (rev_sum>0 & !missing(rev_sum)) 
* Only one democratisation, but with reversals (just not 'only-reversal')
*replace dem_sample=0 if (dem_sum>1 & !missing(dem_sum)) 
*replace dem_sample=0 if (rev_sum!=0 & dem_sum==0)
* More than one democratisation 
*replace dem_sample=0 if (dem_sum<2 & !missing(dem_sum)) 
*replace dem_sample=0 if (rev_sum!=0 & dem_sum==0)

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  
*drop if !e(sample)
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)

tab wbcode if dem_sample==1 & c==1
* N=2,433, N=61 treated countries

gen Y = y/100

gen fectsample = 0
replace fectsample = 1 if (dem_sample==1) & c==1
replace fectsample = 1 if never_dem==1 
* 2,443 (treated) and 1,194 (control) observations

save "$path/DDCGdata_final_fect.dta", replace
	
*** 
*** Estimate the GSynth model
*** 

use "$path/DDCGdata_final_fect.dta", clear

* VEN, UGA and GMB infeasible
local z "ALB ARG BDI BEN BFA BGD BGR BOL BRA BTN CAF CHL CIV COG COM CPV DOM ECU ESP ETH FJI GHA GMB GNB GRC GRD GTM GUY HND HUN IDN KEN KOR LSO MDG MEX MLI MNG MOZ MRT MWI NER NIC NPL PAK PAN PER PHL PRT SDN SEN SLE SUR THA TUR UGA URY ZAF ZMB ZWE"
foreach k of local z{
preserve
	gen fectsample2 = 0
	replace fectsample2 = 1 if wbcode=="`k'"
	replace fectsample2 = 1 if never_dem==1 
	replace fectsample2 = 1 if always_dem==1 
	keep if fectsample2==1
	tsset wbcode2 year
	fect y, treat(dem) unit(wbcode2) time(year) cov(ginv tradewb) method("ife") force("two-way") r(4)
	mat b_`k' = e(ATTs)
restore 	
}
save "$path/DDCGdata_final_fect_i.dta", replace

*** 
*** Plot the results
*** 


use "$path/DDCGdata_final_fect_i.dta", clear

*** ALB
mat ife_iest = b_ALB[22..46,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ALB}") 
graph export "$texfolder/fect_ALB.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** ARG
mat ife_iest = b_ARG[3..41,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ARG}") 
graph export "$texfolder/fect_ARG.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BDI
mat ife_iest = b_BDI[33..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BDI}") 
graph export "$texfolder/fect_BDI.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BEN
mat ife_iest = b_BEN[21..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BEN}") 
graph export "$texfolder/fect_BEN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BFA
mat ife_iest = b_BFA[7..20,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BFA}") 
graph export "$texfolder/fect_BFA.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BGD
mat ife_iest = b_BGD[7..33,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BGD}") 
graph export "$texfolder/fect_BGD.eps", as(eps) replace 
drop event_time ifei_2 ifei_3



*** BGR
mat ife_iest = b_BGR[21..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BGR}") 
graph export "$texfolder/fect_BGR.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BOL
mat ife_iest = b_BOL[12..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BOL}") 
graph export "$texfolder/fect_BOL.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BRA
mat ife_iest = b_BRA[11..47,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BRA}") 
graph export "$texfolder/fect_BRA.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** BTN
mat ife_iest = b_BTN[38..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:BTN}") 
graph export "$texfolder/fect_BTN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** CAF
mat ife_iest = b_CAF[23..43,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:CAF}") 
graph export "$texfolder/fect_CAF.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** CHL
mat ife_iest = b_CHL[7..38,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:CHL}") yscale(range(0 60.1))
graph export "$texfolder/fect_CHL.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** CIV
mat ife_iest = b_CIV[30..42,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:CIV}") yscale(range(0 10.1))
graph export "$texfolder/fect_CIV.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** COG
mat ife_iest = b_COG[19..34,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:COG}") 
graph export "$texfolder/fect_COG.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** COM
mat ife_iest = b_COM[6..24,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (4 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:COM}") 
graph export "$texfolder/fect_COM.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** CPV
mat ife_iest = b_CPV[22..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:CPV}") 
graph export "$texfolder/fect_CPV.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** DOM
mat ife_iest = b_DOM[8..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:DOM}") 
graph export "$texfolder/fect_DOM.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** ECU
mat ife_iest = b_ECU[8..50,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ECU}") 
graph export "$texfolder/fect_ECU.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** ESP
mat ife_iest = b_ESP[8..51,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ESP}") 
graph export "$texfolder/fect_ESP.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** ETH
mat ife_iest = b_ETH[25..50,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ETH}") 
graph export "$texfolder/fect_ETH.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** FJI
mat ife_iest = b_FJI[1..27,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:FJI}") 
graph export "$texfolder/fect_FJI.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** GHA
mat ife_iest = b_GHA[5..30,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GHA}") 
graph export "$texfolder/fect_GHA.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** GNB
mat ife_iest = b_GNB[25..40,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GNB}") 
graph export "$texfolder/fect_GNB.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** GRC
mat ife_iest = b_GRC[1..44,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GRC}") yscale(range(0 20.1))
graph export "$texfolder/fect_GRC.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** GRD
mat ife_iest = b_GRD[10..41,1..3]
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GRD}") 
graph export "$texfolder/fect_GRD.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** GTM
mat ife_iest = b_GTM[2..37,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GTM}") 
graph export "$texfolder/fect_GTM.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** GUY
mat ife_iest = b_GUY[22..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:GUY}") 
graph export "$texfolder/fect_GUY.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** HND
mat ife_iest = b_HND[12..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:HND}") 
graph export "$texfolder/fect_HND.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** HUN
mat ife_iest = b_HUN[20..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:HUN}") 
graph export "$texfolder/fect_HUN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** IDN
mat ife_iest = b_IDN[29..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:IDN}") 
graph export "$texfolder/fect_IDN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** KEN
mat ife_iest = b_KEN[32..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:KEN}") 
graph export "$texfolder/fect_KEN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** KOR
mat ife_iest = b_KOR[17..50,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:KOR}") 
graph export "$texfolder/fect_KOR.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** LSO
mat ife_iest = b_LSO[23..45,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:LSO}") 
graph export "$texfolder/fect_LSO.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** MDG
mat ife_iest = b_MDG[23..49,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MDG}") yscale(range(0 30.1))
graph export "$texfolder/fect_MDG.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** MEX
mat ife_iest = b_MEX[27..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MEX}") 
graph export "$texfolder/fect_MEX.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** MLI
mat ife_iest = b_MLI[22..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MLI}") 
graph export "$texfolder/fect_MLI.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** MNG
mat ife_iest = b_MNG[23..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MNG}") 
graph export "$texfolder/fect_MNG.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** MOZ
mat ife_iest = b_MOZ[24..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MOZ}") 
graph export "$texfolder/fect_MOZ.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** MRT
mat ife_iest = b_MRT[37..48,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MRT}") 
graph export "$texfolder/fect_MRT.eps", as(eps) replace 
drop event_time ifei_2 ifei_3



*** MWI
mat ife_iest = b_MWI[24..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:MWI}") 
graph export "$texfolder/fect_MWI.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** NER
mat ife_iest = b_NER[21..41,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:NER}") 
graph export "$texfolder/fect_NER.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** NIC
mat ife_iest = b_NIC[20..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:NIC}") 
graph export "$texfolder/fect_NIC.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** NPL
mat ife_iest = b_NPL[21..42,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:NPL}") 
graph export "$texfolder/fect_NPL.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** PAK
mat ife_iest = b_PAK[2..23,1..3] 	// 3
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:PAK}") 
graph export "$texfolder/fect_PAK.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** PAN
mat ife_iest = b_PAN[16..43,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:PAN}") 
graph export "$texfolder/fect_PAN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** PER
mat ife_iest = b_PER[2..30,1..3] 	// 4
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (4 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:PER}") 
graph export "$texfolder/fect_PER.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** PHL
mat ife_iest = b_PHL[12..46,1..3] 	// 4
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:PHL}") 
graph export "$texfolder/fect_PHL.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** PRT
mat ife_iest = b_PRT[6..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:PRT}") 
graph export "$texfolder/fect_PRT.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** SDN
mat ife_iest = b_SDN[7..21,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:SDN}") 
graph export "$texfolder/fect_SDN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** SEN
mat ife_iest = b_SEN[30..51,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:SEN}") 
graph export "$texfolder/fect_SEN.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** SLE
mat ife_iest = b_SLE[19..39,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:SLE}") 
graph export "$texfolder/fect_SLE.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** SUR
mat ife_iest = b_SUR[8..34,1..3] 	// 3
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:SUR}") 
graph export "$texfolder/fect_SUR.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** THA
mat ife_iest = b_THA[4..28,1..3] 	// 4
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:THA}") yscale(range(0 10.1))
graph export "$texfolder/fect_THA.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** TUR
mat ife_iest = b_TUR[1..31,1..3] 	// 3
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (3 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:TUR}") 
graph export "$texfolder/fect_TUR.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** URY
mat ife_iest = b_URY[3..39,1..3] 	// 2
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(black*1) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:URY}") 
graph export "$texfolder/fect_URY.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** ZAF
mat ife_iest = b_ZAF[24..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ZAF}") 
graph export "$texfolder/fect_ZAF.eps", as(eps) replace 
drop event_time ifei_2 ifei_3

*** ZMB
mat ife_iest = b_ZMB[21..51,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ZMB}") 
graph export "$texfolder/fect_ZMB.eps", as(eps) replace 
drop event_time ifei_2 ifei_3


*** ZWE
mat ife_iest = b_ZWE[8..27,1..3] 	// 1
mat list ife_iest
svmat ife_iest, names(ifei_)
rename ifei_1 event_time 
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(pink*1.5) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5) xmtick(##5)  ///
	ylabel(, labsize(small) glcolor(gs8) grid glwidth(.2) glp(shortdash)  angle(0)) ///
	xlabel(, labsize(small)) graphregion(color(white)) ///
	yline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xline(0, lpattern(shortdash) lw(.4) lcolor(black))  ///
	xtitle("Years in Democracy (2 events)", size(small)) ///
	ytitle("Democracy Effect (GDPpc increase in %)", size(small)) ///
	title("{bf:ZWE}") 
graph export "$texfolder/fect_ZWE.eps", as(eps) replace 


*** Empty
twoway ///
	(line ifei_3 event_time, sort /// 2
		lw(.8) lcolor(none) lpattern(solid)) /// 
, scheme(s2mono) ysize(5) xsize(5)   ///
	ylabel(, labsize(small) tstyle(none) labcolor(none) glcolor(none) nogrid  angle(0)) ///
	xlabel(, labsize(small) tstyle(none) labcolor(none)  glcolor(none)) graphregion(color(white)) ///
	xtitle("", size(small)) ///
	ytitle("", size(small)) ///
	title("")  yscale(lcolor(none)) xscale(lcolor(none))
graph export "$texfolder/fect_blank.eps", as(eps) replace 
drop event_time ifei_2 ifei_3
