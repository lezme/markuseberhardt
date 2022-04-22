******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                          Figure 1 - Histograms ******************************************************************************
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

global path "/Users/lezme/Dropbox/Curious GMM/Acemoglu Replication"
global path "D:/Dropbox/Curious GMM/Acemoglu Replication"

global dofolder "$path/do_files"
global otherdata "$path/Data"
global outputfolder "$path/Output"
global rubbishfolder "$path/Rubbish"
global texfolder "$path/Tex"

use "$path/DDCGdata_final.dta", clear

*** The following user-written commands need to be installed
*ssc install xtmg 
*ssc install spell

***
*** Data preparation for Figure 1 (a)-(c) - Histograms
***

sort wbcode2 year
xtreg y  dem yy* , fe
gen allevent = 0
replace allevent = 1 if demevent==1 & e(sample)
replace allevent = 1 if revevent==1 & e(sample)

gen redemevent = 0
replace redemevent = 1 if wbcode=="ALB" & year==1997 & e(sample)
replace redemevent = 1 if wbcode=="ARG" & year==1983 & e(sample)
replace redemevent = 1 if wbcode=="BGD" & year==2009 & e(sample)
replace redemevent = 1 if wbcode=="COM" & (year==1996) & e(sample)
replace redemevent = 1 if wbcode=="GHA" & (year==1979) & e(sample)
replace redemevent = 1 if wbcode=="GNB" & (year==1999) & e(sample)
replace redemevent = 1 if wbcode=="GTM" & year==1986 & e(sample)
replace redemevent = 1 if wbcode=="KGZ" & year==2010 & e(sample)
replace redemevent = 1 if wbcode=="LSO" & year==1999 & e(sample)
replace redemevent = 1 if wbcode=="NER" & (year==1999) & e(sample)
replace redemevent = 1 if wbcode=="NGA" & year==1999 & e(sample)
replace redemevent = 1 if wbcode=="NPL" & year==2006 & e(sample)
replace redemevent = 1 if wbcode=="PAK" & (year==1988) & e(sample)
replace redemevent = 1 if wbcode=="PER" & (year==1980) & e(sample)
replace redemevent = 1 if wbcode=="SDN" & year==1986 & e(sample)
replace redemevent = 1 if wbcode=="SLE" & year==2001 & e(sample)
replace redemevent = 1 if wbcode=="SUR" & year==1991 & e(sample)
replace redemevent = 1 if wbcode=="THA" & (year==1978) & e(sample)
replace redemevent = 1 if wbcode=="TUR" & (year==1973) & e(sample)

gen re2demevent = 0
replace re2demevent = 1 if wbcode=="THA" & (year==1992) & e(sample)
replace re2demevent = 1 if wbcode=="TUR" & (year==1983) & e(sample)
replace re2demevent = 1 if wbcode=="PAK" & (year==2008) & e(sample)
replace re2demevent = 1 if wbcode=="PER" & (year==1993) & e(sample)
replace re2demevent = 1 if wbcode=="NER" & (year==2010) & e(sample)
replace re2demevent = 1 if wbcode=="COM" & (year==2002) & e(sample)
replace re2demevent = 1 if wbcode=="GHA" & (year==1996) & e(sample)
replace re2demevent = 1 if wbcode=="GNB" & (year==2005) & e(sample)

gen re3demevent = 0
replace re3demevent = 1 if wbcode=="THA" & (year==2008) & e(sample)

gen all2event = 0
replace all2event = 1 if demevent==1
replace all2event = 1 if wbcode=="BGD" & year==2007 & e(sample)
replace all2event = 1 if wbcode=="COM" & (year==1995) & e(sample)
replace all2event = 1 if wbcode=="FJI" & year==2006 & e(sample)
replace all2event = 1 if wbcode=="GHA" & year==1981 & e(sample)
replace all2event = 1 if wbcode=="GNB" & year==2003 & e(sample)
replace all2event = 1 if wbcode=="HTI" & year==2000 & e(sample)
replace all2event = 1 if wbcode=="NER" & year==2009 & e(sample)
replace all2event = 1 if wbcode=="NGA" & year==1984 & e(sample)
replace all2event = 1 if wbcode=="PAK" & year==1999 & e(sample)
replace all2event = 1 if wbcode=="PER" & (year==1968) & e(sample)
replace all2event = 1 if wbcode=="SDN" & year==1989 & e(sample)
replace all2event = 1 if wbcode=="SLE" & year==1997 & e(sample)
replace all2event = 1 if wbcode=="SUR" & year==1990 & e(sample)
replace all2event = 1 if wbcode=="THA" & (year==1991 | year==2006) & e(sample)
replace all2event = 1 if wbcode=="TUR" & year==1980 & e(sample)

gen all3event = 0
replace all3event = 1 if demevent==1
replace all3event = 1 if wbcode=="COM" & (year==1999) & e(sample)
replace all3event = 1 if wbcode=="PER" & (year==1992) & e(sample)
replace all3event = 1 if wbcode=="THA" & (year==2006) & e(sample)

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

xtreg y dem yy* tradewb ginv if dem_sample==1, fe 
* dem_sample indicates countries which did *switch* into democracy or autocracy during the sample
drop c csum

label var never_dem "Country never a democracy during sample period"
label var always_dem "Country always a democracy during sample period"
label var dem_sample "Country transitioned into or from democracy during sample period"

* Run a dynamic C&K MG regression (to mark the sample)
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

***
*** Figure 1 (a) - Histograms
***

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1  
tab year if e(sample)
* 1964-2010
tabstat wbcode2 if e(sample), by(wbcode) statistic(mean)
* All estimates identified
mat list e(betas) 
gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum demevent

xtmg y dem ginv tradewb ///
		l0ddem l1ddem l2ddem  ///
		l0dginv l1dginv l2dginv   ///
		l0dtradewb l1dtradewb l2dtradewb   ///
		l0yNT ///
		l0ginvNT l1ginvNT l2ginvNT l3ginvNT  ///
		l0tradewbNT l1tradewbNT l2tradewbNT l3tradewbNT   ///
if dem_sample==1 & c==1, robust

twoway ///
	(hist year if allevent==1 & e(sample), discrete frequency lcolor(black) fcolor(midblue*.5) lwidth(.2)) ///
	(hist year if demevent==1 & e(sample), discrete frequency lcolor(black) fcolor(midblue*1.75) lwidth(.2)) ///
	(hist year if allevent==1 & e(sample), discrete frequency lcolor(none) fcolor(none) yaxis(2)) ///
	(hist year if demevent==1 & e(sample), discrete frequency lcolor(none) fcolor(none) yaxis(2)) ///
, legend(order(- "Democratisations (Total 78)" 2 - "Reversals (Total 42)" 1) label(1 "")  label(2 "") ///
	rows(1) size(medsmall)) ymtick(##5, axis(2)) /// 
	scheme(s2mono) graphregion(color(white)) xtitle("Year") ///
	ytitle("Count (democratisations or reversals)", axis(2)) ///
	ytitle("Count (democratisations or reversals)") ylabel(0(2)10, grid glc(black) glw(.2) glp(shortdash)) ///
	xlabel(1960(5)2010) ymtick(##2)	tmtick(##5) ysize(4) xsize(10) ///
	xline(1970.5 1980.5 1990.5 2000.5, lpattern(shortdash) lw(.2) lcolor(black)) ///
	ttext(9.5 1965 "2+3 events", placement(center)) ///
	ttext(9 1965 "(Democratisations", placement(center)) ///
	ttext(8.5   1965 "   + Reversals)", placement(center)) ///
	ttext(9.5 1975.5 "13+11 events", placement(center)) ///
	ttext(9.5 1985.5 "18+6 events", placement(center)) ///
	ttext(9.5 1995.5 "35+12 events", placement(center)) ///
	ttext(9.5 2006 "10+10 events", placement(center)) 
graph export "$texfolder/demo_reversal_events_2021_simple_dyn.eps", as(eps)  replace
*** THIS PLOTS IS INCLUDED IN THE PAPER: Figure 1(a)

 
***
*** Figure 1 (c) - Years spent in democracy by country
***

sort wbcode2 year
by wbcode2: egen dem_sum = sum(dem) if e(sample)
tab dem_sum if e(sample)

preserve
collapse (mean)  dem_sum, by(wbcode2)
twoway ///
	(hist dem_sum, discrete frequency lcolor(black) fcolor(emerald*0.9) lwidth(.2)) ///
	(scatteri 5.1 19 -0.1 19, recast(line) lpattern(solid) lw(1.5) lcolor(white)) ///
	(scatteri 5.1 19 -0.1 19, recast(line) lpattern(dash) lw(.5) lcolor(midblue*1.5)) ///
	(hist dem_sum, discrete frequency lcolor(none) fcolor(none) lwidth(.2) yaxis(2)) ///
	(scatteri 5.1 19 -0.1 19, recast(line) lpattern(solid) lw(.5) lcolor(none) yaxis(2)) ///
, legend(order(- "Countries" 1 - "Total Years in Democracy" 3) label(1 "Count")   label(3 "Median")  ///
	rows(1) size(medsmall)) ymtick(##2, axis(2)) /// 
	scheme(s2mono) graphregion(color(white)) xtitle("Total Years in Democracy") ///
	ytitle("Frequency Count", axis(2)) ylabel(0(1)5, axis(2)) yscale(range(0 5.25)) yscale(range(0 5.25) axis(2)) ///
	ytitle("Frequency Count") ylabel(0(1)5, grid glc(black) glw(.2) glp(shortdash)) ///
	xlabel(5(5)45) ymtick(##2)	xmtick(##5) ysize(4) xsize(10) xscale(range(0.5 36.5)) ///
	text(4.5 5.25 "12 Countries", justification(center) placement(0)) ///
	text(4.5 15.25 "26", justification(center) placement(0)) ///
	text(4.5 25.25 "11", justification(center) placement(0)) ///
	text(4.5 35.25 "10", justification(center) placement(0)) ///
	text(4.5 43.25 "2", justification(center) placement(0)) ///
	xline(10.5 20.5 30.5 40.5, lpattern(shortdash) lw(.2) lcolor(black)) 
graph export "$texfolder/demo_reversal_all_years_in_democracy_dyn.eps", as(eps)  replace
sum dem_sum, de
restore 
*** THIS PLOTS IS INCLUDED IN THE PAPER: Figure 1(c)


***
*** Figure 1 (b) - Years spent in democracy by spell
***

sort wbcode2 year 
spell dem if e(sample), by(wbcode2) 
tab _seq if _end==1 & dem==1
sort wbcode2 year
by wbcode2: egen maxspell = max(_spell)

sum _seq if _end==1 & dem==1, de 
* median 12, mean 13.4 years in democracy 
sum _seq if _end==1 & dem==1   & _spell==maxspell, de 
* Lasting democracy: median 18, mean 18 years in democracy
sum _seq if _end==1 & dem==1   & _spell!=maxspell, de 
* Overturned democracy: median 5, mean 8.1 years in democracy
tab _seq if _end==1 & dem==1
* one lasting democracy of length 45 years is omitted for ease of illustration
tab _seq  if _end==1 & dem==1 

twoway ///
	(hist _seq if _end==1 & dem==1 & _seq<40, discrete frequency lcolor(black) fcolor(pink*.3) lwidth(.2)) ///
	(hist _seq if _end==1 & dem==1 & _spell==maxspell  & _seq<40, discrete frequency lcolor(black) fcolor(pink*1.25) lwidth(.2)) ///
	(scatteri 9.4 11.5 -0.4 11.5, recast(line) lpattern(solid) lw(1.5) lcolor(white)) ///
	(scatteri 9.4 11.5 -0.4 11.5, recast(line) lpattern(dash) lw(.5) lcolor(midblue*1.5)) ///
	(scatteri 9.4 17.5 -0.4 17.5, recast(line) lpattern(solid) lw(1.5) lcolor(white)) ///
	(scatteri 9.4 17.5 -0.4 17.5, recast(line) lpattern(solid) lw(.5) lcolor(pink*1.75)) ///
	(hist _seq if _end==1 & dem==1  & _seq<40, discrete frequency lcolor(none) fcolor(none) lwidth(.2) yaxis(2)) ///
	(scatteri 9.4 17.5 -0.4 17.5, recast(line) lpattern(solid) lw(.5) lcolor(none) yaxis(2)) ///
, legend(order(- "Democratisations" 1 2 - "Median" 4 6) label(2 "Lasting")  label(1 "Overturned")   label(4 "Full Sample")  label(6 "Lasting Sample") ///
	rows(1) size(medsmall)) ymtick(##2, axis(2)) /// 
	scheme(s2mono) graphregion(color(white)) xtitle("Years in Democracy per spell") ///
	ytitle("Frequency Count", axis(2)) ylabel(0(1)9, axis(2)) yscale(range(0 9.01)) yscale(range(0 9.01) axis(2)) ///
	ytitle("Frequency Count") ylabel(0(1)9, grid glc(black) glw(.2) glp(shortdash)) ///
	xlabel(5(5)35) ymtick(##2)	xmtick(##5) ysize(4) xsize(10) xscale(range(0.5 36.5)) ///
	xline(10.5 20.5 30.5, lpattern(shortdash) lw(.2) lcolor(black)) ///
	text(8.5 5.25 "29 Overturned + 11 Lasting", justification(center) placement(0)) ///
	text(8.5 15.25 "10 + 19", justification(center) placement(0)) ///
	text(8.5 25.25 "1 + 13", justification(center) placement(0)) ///
	text(8.5 34.25 "1* + 5", justification(center) placement(0)) 
graph export "$texfolder/demo_reversal_events_years_in_democracy_dyn.eps", as(eps)  replace
*** THIS PLOTS IS INCLUDED IN THE PAPER: Figure 1(b)
