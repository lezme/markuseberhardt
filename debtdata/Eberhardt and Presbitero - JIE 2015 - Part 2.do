******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)		     
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			                 Asymmetric ECM 				
*                    Figure 4 and Appendix Table TA6
******************************************************************************
*			   	       Created: 10th November 2014					
*			   	       Revised: 17th July 2015					
******************************************************************************
* 			                Markus Eberhardt						
******************************************************************************
*  	 		    School of Economics, University of Nottingham,			
*	   		    University Park, Nottingham NG7 2RD, England			
*			    email: markus.eberhardt@nottingham.ac.uk				
*			    web: https://sites.google.com/site/medevecon/			
******************************************************************************
*
* If you use this code, please reference our paper:
*
*    Eberhardt and Presbitero (2015) 'Public debt and growth: Heterogeneity  
*       and non-linearity.' Journal of International Economics, forthcoming.
*
* The usual disclaimers apply. 
*
******************************************************************************
*
* The following commands need to be installed:
* xtmg
* outreg
*
* You also require the two do-files which mark the sample for countries with
* sufficient observations in both regimes (above and below the threshold). 
* These are called comparison60new.do and comparison90new.do -- see lines 219ff
*
******************************************************************************
clear all
set mem 200m
set matsize 11000

global path "C:\Users\lezme\Dropbox"

global data "$path\Debt\Stata"
global graphfolder "$path\Debt\Tex"
global out "$path\Debt\Do"
global output "$path\Debt\Output"
global rubbish "$path\Debt\Rubbish"
global do "$path\TeX\stata"
global do2 "$path\Literature\econometric issues\Nonstationary Panel Metrics\Stata"

***
*** Data construction
***

use "$data\data_imf.dta", clear

***
*** Variable construction
***

label var k "Real Capital stock per worker (US$2000 values)"
label var Debts "Real debt stock (US$2005 values)"
label var debts "Real debt stock per worker (US$2000 values)"
gen lK=log(K)
label var lK "Log of real Capital stock (US$2005 values)"
label var lY "Log of Real GDP (US$2005 values)"
gen lDebts=log(Debts)
label var lDebts "Log of real debt stock (US$2005 values)"
tsset nwbcode year
gen dlY=d.lY
gen dlK=d.lK
gen lL=ln(pop)
label var lL "Log of population"
gen L=pop
gen dlDebts=d.lDebts
gen dly=d.ly
gen dlk=d.lk
gen dldebts=d.ldebts
sort nwbcode year
gen lly=l.ly 
gen llk=l.lk 
gen lldebts=l.ldebts 
gen ldly=l.d.ly 
gen l2dly=l2.d.ly 
gen ldlk=l.d.lk 
gen dlL=d.lL
gen l2dlk=l2.d.lk 
gen ldldebts=l.d.ldebts
gen l2dldebts=l2.d.ldebts 
sort nwbcode year
gen dlYp=dlY*100
gen dlLp=dlL*100
gen dlKp=dlK*100
gen dlDebtsp=dlDebts*100
gen dlyp=dly*100
gen dlkp=dlk*100
gen dldebtsp=dldebts*100
gen gfcf2=inv/gdp
gen Debt2p=100*(Debts/gdp)
gen Debt2=(Debts/gdp)

* N=118 IMF sample
tsset nwbcode year
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
sort year nwbcode
by year: egen dlyT1=mean(dly) if e(sample)
by year: egen dlyT=mean(dlyT1)
by year: egen llyT1=mean(lly) if e(sample)
by year: egen llyT=mean(llyT1)
by year: egen llkT1=mean(llk) if e(sample)
by year: egen llkT=mean(llkT1)
by year: egen lldebtsT1=mean(lldebts) if e(sample)
by year: egen lldebtsT=mean(lldebtsT1)
by year: egen dlkT1=mean(dlk) if e(sample)
by year: egen dlkT=mean(dlkT1)
by year: egen dldebtsT1=mean(dldebts) if e(sample)
by year: egen dldebtsT=mean(dldebtsT1)
by year: egen ldlyT1=mean(ldly) if e(sample)
by year: egen ldlyT=mean(ldlyT1)
by year: egen l2dlyT1=mean(l2dly) if e(sample)
by year: egen l2dlyT=mean(l2dlyT1)
by year: egen ldlkT1=mean(ldlk) if e(sample)
by year: egen ldlkT=mean(ldlkT1)
by year: egen l2dlkT1=mean(l2dlk) if e(sample)
by year: egen l2dlkT=mean(l2dlkT1)
by year: egen ldldebtsT1=mean(ldldebts) if e(sample)
by year: egen ldldebtsT=mean(ldldebtsT1)
by year: egen l2dldebtsT1=mean(l2dldebts) if e(sample)
by year: egen l2dldebtsT=mean(l2dldebtsT1)

drop l2dldebtsT1 ldldebtsT1 l2dlkT1 ldlkT1 l2dlyT1 ldlyT1 dldebtsT1 dlkT1 lldebtsT1 llkT1 llyT1 dlyT1


***
*** Prep for asymmetric ECM
***

* 90% debt threshold
sort nwbcode  year
by nwbcode: gen dldebts_p90=dldebts if dldebts!= . & Debt2p>90 
by nwbcode: gen dldebts_m90=dldebts if dldebts!= . & Debt2p<=90
replace dldebts_p90=0 if missing(dldebts_p90)
replace dldebts_m90=0 if missing(dldebts_m90)
gen ccplus=1 if !missing(ldebts) & Debt2p>90 
gen ccminus=1 if !missing(ldebts) & Debt2p<=90
replace ccplus=0 if missing(ccplus) & !missing(ldebts)
replace ccminus=0 if missing(ccminus) & !missing(ldebts)
by nwbcode: gen ldebtplus=sum(dldebts_p90)
by nwbcode: gen ldebtminus=sum(dldebts_m90)
replace ldebtplus=. if missing(ccplus,ccminus)
replace ldebtminus=. if missing(ccplus,ccminus)
gen lldebts_p90=ldebtplus
gen lldebts_m90=ldebtminus
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
replace lldebts_p90=. if !e(sample)
replace lldebts_m90=. if !e(sample)
drop ldebtplus ldebtminus ldebtmi ldebtpl
drop ccplus ccminus


* 60% debt/GDP threshold (median)
sort nwbcode  year
by nwbcode: gen dldebts_p52=dldebts if dldebts!= . & Debt2p>52
by nwbcode: gen dldebts_m52=dldebts if dldebts!= . & Debt2p<=52
replace dldebts_p52=0 if missing(dldebts_p52)
replace dldebts_m52=0 if missing(dldebts_m52)
gen ccplus=1 if !missing(ldebts) & Debt2p>52
gen ccminus=1 if !missing(ldebts) & Debt2p<=52
replace ccplus=0 if missing(ccplus) & !missing(ldebts)
replace ccminus=0 if missing(ccminus) & !missing(ldebts)
by nwbcode: gen ldebtplus=sum(dldebts_p52)
by nwbcode: gen ldebtminus=sum(dldebts_m52)
replace ldebtplus=. if missing(ccplus,ccminus)
replace ldebtminus=. if missing(ccplus,ccminus)
gen lldebts_p52=ldebtplus
gen lldebts_m52=ldebtminus
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
replace lldebts_p52=. if !e(sample)
replace lldebts_m52=. if !e(sample)
drop ldebtplus ldebtminus ldebtmi ldebtpl
drop ccplus ccminus

* Determining the share of countries with >25% obs above/below threshold

* Mean debt/GDP ratio
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
gen c=1 if e(sample)
gen e=1 if e(sample) & Debt2p>59.90 & !missing(Debt2p)
sort nwbcode year
by nwbcode: egen csum=sum(c)
by nwbcode: egen esum=sum(e)
gen share=esum/csum
list nwbcode csum nwbcode esum share if year==2000, noo separator(0)
drop esum csum share c e

* 90% debt/GDP ratio
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
gen c=1 if e(sample)
gen e=1 if e(sample) & Debt2p>90 & !missing(Debt2p)
sort nwbcode year
by nwbcode: egen csum=sum(c)
by nwbcode: egen esum=sum(e)
gen share=esum/csum
list nwbcode csum nwbcode esum share if year==2000, noo separator(0)
drop esum csum share c e


***
*** Asymmetric ECM: Threshold selection
***

*run "$out\comparison60new.do"
*rename comparison60 comparison
*rename lldebts_p52 lldebts_p 
*rename lldebts_m52 lldebts_m 
*rename dldebts_p52 dldebts_p 
*rename dldebts_m52 dldebts_m 

run "$out\comparison90new.do"
rename comparison90 comparison
rename lldebts_p90 lldebts_p 
rename lldebts_m90 lldebts_m 
rename dldebts_p90 dldebts_p 
rename dldebts_m90 dldebts_m 


***
*** Asymmetric ECM: Estimation
***

run "$out\outreg.ado"
run "$do\xtmg.ado"

* CA are constructed from ALL countries in the sample
xtreg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m, fe
sort year
by year: egen lldebts_pT=mean(lldebts_p) if e(sample)
by year: egen lldebts_mT=mean(lldebts_m) if e(sample)
by year: egen dldebts_pT=mean(dldebts_p) if e(sample)
by year: egen dldebts_mT=mean(dldebts_m) if e(sample)

local z "dldebts_pT dldebts_mT"
foreach k of local z{
	tsset nwbcode year
	gen l1`k'=l.`k'
	gen l2`k'=l.`k'
	gen l3`k'=l.`k'
}

cd "$output"

* Standard MG with trend
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m   if comparison==1, trend robust res(eMGt) 
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma replace

* Standard Asymmetric CMG
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT if comparison==1, robust  res(eAACMG)
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma append

* Trend
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT if comparison==1 , robust trend res(eAACMGt)
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma append

* Long-run asymmetry only
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT dldebtsT if comparison==1 , robust  res(eA0CMG)
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma append

* Standard Asymmetric CCEMG add 1 lag
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT ldlkT  l1dldebts_pT l1dldebts_mT if comparison==1, robust  res(el1AACMG)
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma append

* Standard Asymmetric CCEMG add 2 lag
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT l2dlkT  l2dldebts_pT l2dldebts_mT if comparison==1, robust  res(el2AACMG)
	outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma append

	
*********************
*** MG with trend ***
*********************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m   if comparison==1, trend robust
mat bMGt=e(betas)
mat tMGt=e(tbetas)
svmat bMGt
svmat tMGt
replace bMGt3=. if bMGt3==0
replace bMGt4=. if bMGt4==0
replace bMGt6=. if bMGt6==0
replace bMGt7=. if bMGt7==0
rreg bMGt3, genwt(bMGt3weight) 
rreg bMGt4, genwt(bMGt4weight)

* EC Coeff t-stat and t-bar statistic
*reg bMGt1
reg tMGt1

* LRA results for >cutoff and <cutoff
sureg (eq1: bMGt3) (eq2: bMGt1) [w = bMGt3weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: bMGt4) (eq2: bMGt1) [w = bMGt4weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff and <cutoff
gen ALRdebttp=bMGt3/-bMGt1
gen ALRdebttm=bMGt4/-bMGt1
rreg ALRdebttp
rreg ALRdebttm

xtcd eMGt


*************
*** CCEMG ***
*************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT    if comparison==1, robust 
mat bCMG=e(betas)
mat tCMG=e(tbetas)
svmat bCMG
svmat tCMG
replace bCMG3=. if bCMG3==0
replace bCMG4=. if bCMG4==0
replace bCMG6=. if bCMG6==0
replace bCMG7=. if bCMG7==0
rreg bCMG3, genwt(bCMG3weight)
rreg bCMG4, genwt(bCMG4weight)

* EC Coeff and t-bar statistic
*reg bCMG1
reg tCMG1

* LRA results for >cutoff  and <cutoff
sureg (eq1: bCMG3) (eq2: bCMG1) [w=bCMG3weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: bCMG4) (eq2: bCMG1) [w=bCMG4weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff  and <cutoff
gen ALRdebtpC=bCMG3/-bCMG1
gen ALRdebtmC=bCMG4/-bCMG1
rreg ALRdebtpC
rreg ALRdebtmC

* CD Test
xtcd eAACMG


*************************
*** CCEMG with trends ***
*************************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT if comparison==1, robust trend
mat bCMGt=e(betas)
mat tCMGt=e(tbetas)
svmat bCMGt
svmat tCMGt
replace bCMGt3=. if bCMGt3==0
replace bCMGt4=. if bCMGt4==0
replace bCMGt6=. if bCMGt6==0
replace bCMGt7=. if bCMGt7==0
rreg bCMGt3, genwt(bCMGt3weight)
rreg bCMGt4, genwt(bCMGt4weight)

* EC Coeff and t-bar statistic
reg bCMGt1
reg tCMGt1

* LRA results for >cutoff  and <cutoff
sureg (eq1: bCMGt3) (eq2: bCMGt1) [w=bCMGt3weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: bCMGt4) (eq2: bCMGt1) [w=bCMGt4weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff  and <cutoff
gen ALRtdebtpC=bCMGt3/-bCMGt1
gen ALRtdebtmC=bCMGt4/-bCMGt1
rreg ALRtdebtpC
rreg ALRtdebtmC

* CD Test
xtcd eAACMGt

******************************************
*** CCEMG with long-run asymmetry only ***
******************************************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT dldebtsT if comparison==1, robust 
mat b2CMG=e(betas)
mat t2CMG=e(tbetas)
svmat b2CMG
svmat t2CMG
replace b2CMG3=. if b2CMG3==0
replace b2CMG4=. if b2CMG4==0
replace b2CMG6=. if b2CMG6==0
rreg b2CMG3, genwt(bCMGA03weight)
rreg b2CMG4, genwt(bCMGA04weight)

* EC Coeff and t-bar statistic
reg b2CMG1
reg t2CMG1

* LRA results for >cutoff  and <cutoff
sureg (eq1: b2CMG3) (eq2: b2CMG1) [w=bCMGA03weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: b2CMG4) (eq2: b2CMG1) [w=bCMGA04weight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff  and <cutoff
gen ALRdebtpC2=b2CMG3/-b2CMG1
gen ALRdebtmC2=b2CMG4/-b2CMG1
rreg ALRdebtpC2
rreg ALRdebtmC2

* CD Test
xtcd eA0CMG


*******************************************
*** CCEMG with one additional lagged CA ***
*******************************************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT ldlkT  l1dldebts_pT l1dldebts_mT if comparison==1, robust
mat b2CMGL=e(betas)
mat t2CMGL=e(tbetas)
svmat b2CMGL
svmat t2CMGL
replace b2CMGL3=. if b2CMGL3==0
replace b2CMGL4=. if b2CMGL4==0
replace b2CMGL6=. if b2CMGL6==0
replace b2CMGL7=. if b2CMGL7==0
rreg b2CMGL3, genwt(bCMGA03Lweight)
rreg b2CMGL4, genwt(bCMGA04Lweight)

* EC Coeff and t-bar statistic
reg b2CMGL1
reg t2CMGL1

* LRA results for >cutoff  and <cutoff
sureg (eq1: b2CMGL3) (eq2: b2CMGL1) [w=bCMGA03Lweight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: b2CMGL4) (eq2: b2CMGL1) [w=bCMGA04Lweight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff  and <cutoff
gen ALRdebtpC2L=b2CMGL3/-b2CMGL1
gen ALRdebtmC2L=b2CMGL4/-b2CMGL1
rreg ALRdebtpC2L
rreg ALRdebtmC2L

* CD Test
xtcd el1AACMG




*******************************************
*** CCEMG with two additional lagged CA ***
*******************************************
xtmg dly lly llk lldebts_p lldebts_m dlk  dldebts_p dldebts_m  ///
		dlyT llyT llkT lldebts_pT lldebts_mT dlkT  dldebts_pT dldebts_mT l2dlkT  l2dldebts_pT l2dldebts_mT if comparison==1, robust
mat b2CMGLL=e(betas)
mat t2CMGLL=e(tbetas)
svmat b2CMGLL
svmat t2CMGLL
replace b2CMGLL3=. if b2CMGLL3==0
replace b2CMGLL4=. if b2CMGLL4==0
replace b2CMGLL6=. if b2CMGLL6==0
replace b2CMGLL7=. if b2CMGLL7==0
rreg b2CMGLL3, genwt(bCMGA03LLweight)
rreg b2CMGLL4, genwt(bCMGA04LLweight)

* EC Coeff and t-bar statistic
reg b2CMGLL1
reg t2CMGLL1

* LRA results for >cutoff  and <cutoff
sureg (eq1: b2CMGLL3) (eq2: b2CMGLL1) [w=bCMGA03LLweight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]
sureg (eq1: b2CMGLL4) (eq2: b2CMGLL1) [w=bCMGA04LLweight]
nlcom [eq1]_b[_cons]/-[eq2]_b[_cons]

* ALR results for >cutoff  and <cutoff
gen ALRdebtpC2LL=b2CMGLL3/-b2CMGLL1
gen ALRdebtmC2LL=b2CMGLL4/-b2CMGLL1
rreg ALRdebtpC2LL
rreg ALRdebtmC2LL

* CD Test
xtcd el2AACMG


**************
** PICTURES **
**************

* We pick two specifications for 59.9% cutoff: standard CMG and CMG with 2 lags
* ALRdebtpC ALRdebtmC ALRdebtpC2LL ALRdebtmC2LL
* We also pick the same two specifications for 90% cutoff

preserve
keep	ALRdebtpC ALRdebtmC ALRdebtpC2LL ALRdebtmC2LL 
gen list=_n
*save "$data\asymcoeff60new.dta", replace
save "$data\asymcoeff90new.dta", replace
restore

exit

gen ldebtplus=lldebts if Debt2p>90
gen ldebtminus=lldebts if Debt2p<90
gen totdebtplus=Debt2p if Debt2p>90
gen totdebtminus=Debt2p if Debt2p<90

*gen ldebtplus=lldebts if Debt2p>59.9
*gen ldebtminus=lldebts if Debt2p<59.9
*gen totdebtplus=Debt2p if Debt2p>59.9
*gen totdebtminus=Debt2p if Debt2p<59.9

drop if comparison==0
keep totdebtplus totdebtminus ldebtplus ldebtminus comparison nwbcode pop lk Debt2p debts gdppc gdp lldebts_p lldebts_m
collapse  (mean) totdebtplus totdebtminus ldebtplus ldebtminus comparison  pop lk Debt2p debts gdppc gdp lldebts_p lldebts_m ///
	(max)  debt_peak=Debt2p, by(nwbcode)
gen list=_n
decode nwbcode, gen(wbcode)
sort list
merge list using  "$data\asymcoeff90new.dta", sort
save "$data\asymcoeffdat90.dta", replace
keep if _merge==3
drop _merge

* 59.9%
*use "$data\asymcoeffdat60.dta", clear
* 90%
use "$data\asymcoeffdat90.dta", clear

order nwbcode wbcode list

gen lgdppc=log(gdppc)
gen ldebts=log(debts)
gen ltotdebtplus=log(totdebtplus)
gen ltotdebtminus=log(totdebtminus)
drop if missing(wbcode)

******************************************************************************************************************
* These are the arrow graphs for 60% cutoffs
******************************************************************************************************************
* Standard CMG
******************************************************************************************************************
twoway ///  
	(pcbarrow ALRdebtmC ltotdebtminus  ALRdebtpC ltotdebtplus if comparison==1 & wbcode!="RWA",  ///
		msize(large) mlw(.4) mcolor(gray) mfcolor(blue)  mlcolor(blue) lwidth(.4) ///
		lpattern(solid) lcolor(blue)) ///
	(scatter ALRdebtmC ltotdebtminus if wbcode=="ZAR" | wbcode=="JOR" | wbcode=="NGA", ///
		msymbol(none) mlab(wbcode) mlabposition(8)) ///
	, legend(off) ytitle("Long-run debt coefficients") xtitle("Average debt/GDP ratio (in logs)") ///
	 scheme(s2mono) title("60% debt/GDP threshold") subtitle("N=54 countries; positive slope in 28/54 countries") ///
	xscale(range(2.5 6)) xlabel(2.5(.5)6) yline(0, lcolor(black) lwidth(.2)) 	///
	text(-.85 4.5 "RWA omitted to aid illustration", place(e)  just(right)) ///
	text(1.4 4 "Coefficient Change: mean -.006 [t=0.06];", place(e)  just(left) size(medsmall)) ///
	text(1.25 4.25 "robust mean .016 [t=0.54]", place(e)  just(left) size(medsmall))
graph export "$graphfolder\LRcomparison60CMG.png", as(png) replace 
******************************************************************************************************************
* CMG with two lags
******************************************************************************************************************
twoway ///  
	(pcbarrow ALRdebtmC2LL ltotdebtminus  ALRdebtpC2LL ltotdebtplus if comparison==1 & wbcode!="CMR"   ///
		& wbcode!="RWA" & wbcode!="JOR" , msize(large) mlw(.4) mcolor(gray) mfcolor(blue)  mlcolor(blue) lwidth(.4) ///
		lpattern(solid) lcolor(blue)) ///
	(scatter ALRdebtmC2LL ltotdebtminus if wbcode=="URY" | wbcode=="USA" | wbcode=="NGA" , ///
		msymbol(none) mlab(wbcode) mlabposition(8)) ///
	, legend(off) ytitle("Long-run debt coefficients") xtitle("Average debt/GDP ratio (in logs)") ///
	 scheme(s2mono) title("60% debt/GDP threshold") subtitle("N=53 countries; positive slope in 24/53 countries") ///
	xscale(range(2.5 6)) xlabel(2.5(.5)6) yline(0, lcolor(black) lwidth(.2)) 	///
	text(-.4 2.35 "CMR, JOR & RWA omitted to aid illustration", place(e)  just(right)) ///
	text(.9 2.5 "Coefficient Change: mean -.130 [t=0.99]; robust mean .010 [t=0.40]", place(e) just(left) size(medsmall)) 
graph export "$graphfolder\LRcomparison60CMGLL.png", as(png) replace 
******************************************************************************************************************
gen plus=1 if ALRdebtmC<ALRdebtpC
replace plus=0 if missing(plus) & !missing(ALRdebtmC)
tab plus
drop plus
gen plus=1 if ALRdebtmC2LL<ALRdebtpC2LL
replace plus=0 if missing(plus) & !missing(ALRdebtmC2LL)
tab plus
******************************************************************************************************************
* CMG: 26 countries have a negative slope, 28 a positive one
* CMG2LL: 24 countries have a negative slope, 29 a positive one
******************************************************************************************************************
gen coeffdiff=ALRdebtmC-ALRdebtpC
reg coeffdiff if coeffdiff>0
reg coeffdiff if coeffdiff<0
rreg coeffdiff if coeffdiff>0
rreg coeffdiff if coeffdiff<0
reg coeffdiff
rreg coeffdiff
drop coeffdiff
gen coeffdiff=ALRdebtmC2LL-ALRdebtpC2LL
reg coeffdiff
rreg coeffdiff
drop coeffdiff


******************************************************************************************************************
* These are the arrow graphs for 90% cutoffs
******************************************************************************************************************
* Standard CMG
******************************************************************************************************************
twoway ///  
	(pcbarrow ALRdebtmC ltotdebtminus  ALRdebtpC ltotdebtplus if comparison==1 ///
		& wbcode!="TZA" & wbcode!="ZAR",  ///
		msize(large) mlw(.4) mcolor(gray) mfcolor(red)  mlcolor(red) lwidth(.4) ///
		lpattern(solid) lcolor(red)) ///
	(scatter ALRdebtpC ltotdebtplus if wbcode=="HUN" | wbcode=="NGA" | wbcode=="TGO"  | wbcode=="BDI",  ///
		msymbol(none) mlab(wbcode) mlabposition(3)) ///
	, legend(off) ytitle("Long-run debt coefficients") xtitle("Average debt/GDP ratio (in logs)") ///
	 scheme(s2mono) title("90% debt/GDP threshold ({&beta}{sub:CMG})") subtitle("N=28 countries; positive slope in 14/28 countries") ///
	xscale(range(3 6)) xlabel(3(.5)6) yline(0, lcolor(black) lwidth(.2)) 	///
	text(-.9 4.5 "TZA, ZAR omitted to aid illustration", place(e)  just(right)) ///
	text(-.6 3 "Coefficient Change: mean .150 [t=1.62];", place(e)  just(left) size(medsmall)) ///
	text(-.7 3.5 "robust mean .020 [t=0.32]", place(e)  just(left) size(medsmall))
graph export "$graphfolder\LRcomparison90CMG.png", as(png) replace 
******************************************************************************************************************
* CMG with two lags
******************************************************************************************************************
twoway ///  
	( pcbarrow ALRdebtmC2LL ltotdebtminus  ALRdebtpC2LL ltotdebtplus, ///
		msize(large) mlw(.4) mcolor(gray) mfcolor(red)  mlcolor(red) lwidth(.4) lpattern(solid) lcolor(red) ) ///
	( scatter ALRdebtmC2LL ltotdebtminus if wbcode=="NGA" | wbcode=="CIV" | wbcode=="TZA", ///
		msymbol(none) mlab(wbcode) mlabposition(8) ) ///
	( scatter ALRdebtpC2LL ltotdebtplus if wbcode=="JOR" | wbcode=="HUN", ///
		msymbol(none) mlab(wbcode) mlabposition(2) ) ///
	, legend(off) ytitle("Long-run debt coefficients") xtitle("Average debt/GDP ratio (in logs)") ///
	 scheme(s2mono) title("90% debt/GDP threshold ({&beta}{sub:CMG 2 lags})") subtitle("N=27 countries; positive slope in 14/27 countries") ///
	xscale(range(3 6)) xlabel(3(.5)6) yline(0, lcolor(black) lwidth(.2)) 	///
	text(-.7 3 "Coefficient Change: mean .008 [t=0.09]; ", place(e) just(left) size(medsmall))  ///
	text(-.85 3.25 "robust mean -.068 [t=-0.95]", place(e) just(left) size(medsmall)) 
graph export "$graphfolder\LRcomparison90CMGLL.png", as(png) replace 
******************************************************************************************************************
gen plus=1 if ALRdebtmC<ALRdebtpC
replace plus=0 if missing(plus) & !missing(ALRdebtmC)
tab plus
drop plus
gen plus=1 if ALRdebtmC2LL<ALRdebtpC2LL
replace plus=0 if missing(plus) & !missing(ALRdebtmC2LL)
tab plus
******************************************************************************************************************
* CMG: 14 countries have a negative slope, 14 a positive one
* CMG2LL: 13 countries have a negative slope, 14 a positive one
******************************************************************************************************************
gen coeffdiff=ALRdebtmC-ALRdebtpC
reg coeffdiff if coeffdiff>0
reg coeffdiff if coeffdiff<0
rreg coeffdiff if coeffdiff>0
rreg coeffdiff if coeffdiff<0
reg coeffdiff
rreg coeffdiff
drop coeffdiff
gen coeffdiff=ALRdebtmC2LL-ALRdebtpC2LL
reg coeffdiff
rreg coeffdiff
drop coeffdiff
