******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			                 Co-Summability
*                         Appendix Table TA-9				
******************************************************************************
*			   	       Created: 14th November 2014					
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
* You require the do-file cosummhelp118.do to do the subsampling.
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

/*

* You only need to run everything that's commented out ONCE: these are the static
* regressions with levels, squared and cubed debt terms, from which the residuals
* are collected. At the end of this section, around line 446, these residuals are 
* stored in a separate file, residuals.dta, which can then be used for the co-
* summability analysis.

***
*** Data construction
***

use "$data\data_imf.dta", clear


***
*** Variable construction
***

gen ldebts2=ldebts^2
gen ldebts3=ldebts^3
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
gen l2ly=l2.ly 
gen l2lk=l2.lk 
gen l2ldebts=l2.ldebts 
gen l2ldebts2=l2.ldebts2
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
gen lldebts2=l.ldebts2


* N=118 IMF sample
tsset nwbcode year
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
sort year nwbcode
by year: egen lyT1=mean(ly) if e(sample)
by year: egen lyT=mean(lyT1)
by year: egen dlyT1=mean(dly) if e(sample)
by year: egen dlyT=mean(dlyT1)
by year: egen l2lyT1=mean(l2ly) if e(sample)
by year: egen l2lyT=mean(l2lyT1)
by year: egen llyT1=mean(lly) if e(sample)
by year: egen llyT=mean(llyT1)
by year: egen lkT1=mean(lk) if e(sample)
by year: egen lkT=mean(lkT1)
by year: egen l2lkT1=mean(l2lk) if e(sample)
by year: egen l2lkT=mean(l2lkT1)
by year: egen llkT1=mean(llk) if e(sample)
by year: egen llkT=mean(llkT1)
by year: egen ldebtsT1=mean(ldebts) if e(sample)
by year: egen ldebtsT=mean(ldebtsT1)
by year: egen ldebts2T1=mean(ldebts2) if e(sample)
by year: egen ldebts2T=mean(ldebts2T1)
by year: egen lldebts2T1=mean(lldebts2) if e(sample)
by year: egen lldebts2T=mean(lldebts2T1)
by year: egen l2ldebts2T1=mean(l2ldebts2) if e(sample)
by year: egen l2ldebts2T=mean(l2ldebts2T1)
by year: egen lldebtsT1=mean(lldebts) if e(sample)
by year: egen lldebtsT=mean(lldebtsT1)
by year: egen l2ldebtsT1=mean(l2ldebts) if e(sample)
by year: egen l2ldebtsT=mean(l2ldebtsT1)
tsset nwbcode year

***
*** Empirical models
***

* (1) MG
xtmg ly lk ldebts, trend res(eMGt)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: btrend 	= testbeta[1..118,3]#J(53,1,1)
mata: bcons 		= testbeta[1..118,4]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_trend=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_trend",btrend)
mata: st_store(.,"b_cons",bcons)
qui gen eMG1=ly-(b_k*lk+b_debt*ldebts) if e(sample)
drop b_k b_debt b_trend b_cons


* (2) MG with squared term
xtmg ly lk ldebts ldebts2, trend res(eMG2t)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: bldebtsq 	= testbeta[1..118,3]#J(53,1,1)
mata: btrend 	= testbeta[1..118,4]#J(53,1,1)
mata: bcons 	= testbeta[1..118,5]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_debtsq=0
gen b_trend=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_debtsq",bldebtsq)
mata: st_store(.,"b_trend",btrend)
mata: st_store(.,"b_cons",bcons)
qui gen eMG2=ly-(b_k*lk+b_debt*ldebts+b_debtsq*ldebts2) if e(sample)
drop b_k b_debt b_debtsq b_trend b_cons



* (3) CMG (full averages)
xtmg ly lk ldebts lyT lkT ldebtsT, robust res(eCMG)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: blybar	= testbeta[1..118,3]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,4]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,5]#J(53,1,1)
mata: bcons 	= testbeta[1..118,6]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG1=ly-(b_k*lk+b_debt*ldebts+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT) if e(sample)
drop b_k b_debt  b_cons b_ybar b_kbar b_debtbar 



* (4) CMG with squared term (full averages)
xtmg ly lk ldebts ldebts2  lyT lkT ldebtsT ldebts2T, robust res(e2CMG)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: bldebtsq 	= testbeta[1..118,3]#J(53,1,1)
mata: blybar	= testbeta[1..118,4]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,5]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,6]#J(53,1,1)
mata: bldebt2bar= testbeta[1..118,7]#J(53,1,1)
mata: bcons 	= testbeta[1..118,8]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_debtsq=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_debt2bar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_debtsq",bldebtsq)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_debt2bar",bldebt2bar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG2=ly-(b_k*lk+b_debt*ldebts+b_debtsq*ldebts2+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT+b_debt2bar*ldebts2T) if e(sample)
drop b_k b_debt b_debtsq  b_cons b_ybar b_kbar b_debtbar b_debt2bar  



* (5) CMG with one additional lag of CAs - linear
xtmg ly lk ldebts lyT lkT ldebtsT llyT llkT lldebtsT, robust res(eCMG1_)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: blybar	= testbeta[1..118,3]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,4]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,5]#J(53,1,1)
mata: bl2ybar 	= testbeta[1..118,6]#J(53,1,1)
mata: bl2kbar	= testbeta[1..118,7]#J(53,1,1)
mata: bl2debtbar= testbeta[1..118,8]#J(53,1,1)
mata: bcons 	= testbeta[1..118,9]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_lybar=0
gen b_lkbar=0
gen b_ldebtbar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_lybar",bl2ybar)
mata: st_store(.,"b_lkbar",bl2kbar)
mata: st_store(.,"b_ldebtbar",bl2debtbar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG1c=ly-(b_k*lk+b_debt*ldebts+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT+b_lybar*llyT+b_lkbar*llkT+b_ldebtbar*lldebtsT) if e(sample)
drop b_k b_debt b_cons b_ybar b_kbar b_debtbar  b_ldebtbar b_lkbar b_lybar


* (6) CMG with one additional lag of CAs - squared
xtmg ly lk ldebts ldebts2 lyT lkT ldebtsT ldebts2T llyT llkT lldebtsT lldebts2T, robust res(e2CMG1_)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2 	= testbeta[1..118,3]#J(53,1,1)
mata: blybar	= testbeta[1..118,4]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,5]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,6]#J(53,1,1)
mata: bldebt2bar= testbeta[1..118,7]#J(53,1,1)
mata: bl2ybar 	= testbeta[1..118,8]#J(53,1,1)
mata: bl2kbar	= testbeta[1..118,9]#J(53,1,1)
mata: bl2debtbar= testbeta[1..118,10]#J(53,1,1)
mata: bl2debt2bar= testbeta[1..118,11]#J(53,1,1)
mata: bcons 	 = testbeta[1..118,12]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_debt2=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_debt2bar=0
gen b_lybar=0
gen b_lkbar=0
gen b_ldebtbar=0
gen b_ldebt2bar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_debt2",bldebt2)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_debt2bar",bldebt2bar)
mata: st_store(.,"b_lybar",bl2ybar)
mata: st_store(.,"b_lkbar",bl2kbar)
mata: st_store(.,"b_ldebtbar",bl2debtbar)
mata: st_store(.,"b_ldebt2bar",bl2debt2bar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG2c=ly-(b_k*lk+b_debt*ldebts+b_debt2*ldebts2+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT+b_debt2bar*ldebts2T+b_lybar*llyT+b_lkbar*llkT+b_ldebt2bar*lldebts2T+b_ldebt2bar*lldebts2T) if e(sample)
drop b_k b_debt b_debt2 b_cons b_ybar b_kbar b_debtbar b_debt2bar  b_ldebtbar b_ldebt2bar b_lkbar b_lybar



* (7) CMG with two additional lags of CAs - linear
xtmg ly lk ldebts lyT lkT ldebtsT llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT, robust res(eCMG2_)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: blybar	= testbeta[1..118,3]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,4]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,5]#J(53,1,1)
mata: bl2ybar 	= testbeta[1..118,6]#J(53,1,1)
mata: bl2kbar	= testbeta[1..118,7]#J(53,1,1)
mata: bl2debtbar= testbeta[1..118,8]#J(53,1,1)
mata: bl3ybar 	= testbeta[1..118,9]#J(53,1,1)
mata: bl3kbar	= testbeta[1..118,10]#J(53,1,1)
mata: bl3debtbar= testbeta[1..118,11]#J(53,1,1)
mata: bcons 	= testbeta[1..118,12]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_lybar=0
gen b_lkbar=0
gen b_ldebtbar=0
gen b_l2ybar=0
gen b_l2kbar=0
gen b_l2debtbar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_lybar",bl2ybar)
mata: st_store(.,"b_lkbar",bl2kbar)
mata: st_store(.,"b_ldebtbar",bl2debtbar)
mata: st_store(.,"b_l2ybar",bl3ybar)
mata: st_store(.,"b_l2kbar",bl3kbar)
mata: st_store(.,"b_l2debtbar",bl3debtbar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG1cc=ly-(b_k*lk+b_debt*ldebts+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT+b_lybar*llyT+b_lkbar*llkT+b_ldebtbar*lldebtsT+b_l2ybar*l2lyT+b_l2kbar*l2lkT+b_l2debtbar*l2ldebtsT) if e(sample)
drop b_k b_debt b_cons b_ybar b_kbar b_debtbar  b_ldebtbar b_lkbar b_lybar b_l2debtbar b_l2kbar b_l2ybar


* (8) CMG with two additional lags of CAs - squared
xtmg ly lk ldebts ldebts2 lyT lkT ldebtsT ldebts2T llyT llkT lldebtsT lldebts2T l2lyT l2lkT l2ldebtsT l2ldebts2T, robust res(e2CMG2_)
mat testb=e(betas)
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk	 	= testbeta[1..118,1]#J(53,1,1)
mata: bldebt 	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2 	= testbeta[1..118,3]#J(53,1,1)
mata: blybar	= testbeta[1..118,4]#J(53,1,1)
mata: blkbar 	= testbeta[1..118,5]#J(53,1,1)
mata: bldebtbar = testbeta[1..118,6]#J(53,1,1)
mata: bldebt2bar = testbeta[1..118,7]#J(53,1,1)
mata: bl2ybar 	= testbeta[1..118,8]#J(53,1,1)
mata: bl2kbar	= testbeta[1..118,9]#J(53,1,1)
mata: bl2debtbar= testbeta[1..118,10]#J(53,1,1)
mata: bl2debt2bar= testbeta[1..118,11]#J(53,1,1)
mata: bl3ybar 	= testbeta[1..118,12]#J(53,1,1)
mata: bl3kbar	= testbeta[1..118,13]#J(53,1,1)
mata: bl3debtbar= testbeta[1..118,14]#J(53,1,1)
mata: bl3debt2bar= testbeta[1..118,15]#J(53,1,1)
mata: bcons 	= testbeta[1..118,16]#J(53,1,1)
gen b_k=0
gen b_debt=0
gen b_debt2=0
gen b_ybar=0
gen b_kbar=0
gen b_debtbar=0
gen b_debt2bar=0
gen b_lybar=0
gen b_lkbar=0
gen b_ldebtbar=0
gen b_ldebt2bar=0
gen b_l2ybar=0
gen b_l2kbar=0
gen b_l2debtbar=0
gen b_l2debt2bar=0
gen b_cons=0
mata: st_store(.,"b_k",blk)
mata: st_store(.,"b_debt",bldebt)
mata: st_store(.,"b_debt2",bldebt2)
mata: st_store(.,"b_ybar",blybar)
mata: st_store(.,"b_kbar",blkbar)
mata: st_store(.,"b_debtbar",bldebtbar)
mata: st_store(.,"b_debt2bar",bldebt2bar)
mata: st_store(.,"b_lybar",bl2ybar)
mata: st_store(.,"b_lkbar",bl2kbar)
mata: st_store(.,"b_ldebtbar",bl2debtbar)
mata: st_store(.,"b_ldebt2bar",bl2debt2bar)
mata: st_store(.,"b_l2ybar",bl3ybar)
mata: st_store(.,"b_l2kbar",bl3kbar)
mata: st_store(.,"b_l2debtbar",bl3debtbar)
mata: st_store(.,"b_l2debt2bar",bl3debt2bar)
mata: st_store(.,"b_cons",bcons)
qui gen eCMG2cc=ly-(b_k*lk+b_debt*ldebts+b_debt2*ldebts2+b_ybar*lyT+b_kbar*lkT+b_debtbar*ldebtsT+b_debt2bar*ldebts2T+b_lybar*llyT+b_lkbar*llkT+b_ldebtbar*lldebtsT+b_ldebt2bar*lldebts2T+b_l2ybar*l2lyT+b_l2kbar*l2lkT+b_l2debtbar*l2ldebtsT+b_l2debt2bar*l2ldebts2T) if e(sample)
drop b_k b_debt b_debt2 b_cons b_ybar b_kbar b_debtbar b_debt2bar b_ldebtbar b_ldebt2bar b_lkbar b_lybar b_l2debtbar b_l2kbar b_l2ybar b_l2debt2bar


preserve
sort nwbcode year
keep nwbcode year eMG1 eMG2 eCMG1 eCMG2 eCMG1c eCMG2c eCMG1cc eCMG2cc ly lk ldebts
save "$data\residuals.dta", replace
restore

exit
*/

***
*** Co-Summability
***

use "$data\residuals.dta", clear
* See above for the origins of 'residuals.dta'

local z "eMG1 eMG2 eCMG1 eCMG2 eCMG1c eCMG2c eCMG1cc eCMG2cc"
foreach k of local z{
	gen d`k'=d.`k'
}
gen c=0
sort nwbcode year
* Make adjustment between deMG1,deMG2,deCMG1,deCMG2 on the one hand 
* and deCMG1c,deCMG2c or deCMG1cc,deCMG2cc on the other (latter have shorter panel)
by nwbcode: replace c=1 if !missing(deCMG2cc)
by nwbcode: gen csum=sum(c)
drop if csum==0
** this guarantees the first observation for each country is never missing
drop c
gen c=1 if csum==1
gen list=sum(c)
drop c csum
order list
gen c=1
by nwbcode: egen csum=sum(c)
by nwbcode: replace csum=. if _n!=1
mkmat csum, matrix(robs)

mata:
numeric matrix excludemissing(numeric matrix A)
{
    return(select(A, rowmissing(A):==0))
}
end



******************************************************************************************************************************
******************************************************************************************************************************
******************************************************************************************************************************

mata: betas = J(118, 1, 0)
mata: deltas = J(118, 1, 0)

***
*** Pick the residual series to be tested for co-summability
***
// mkmat eMG1, matrix(z0t)
// mkmat eMG2, matrix(z0t)
// mkmat eCMG1, matrix(z0t)
// mkmat eCMG2, matrix(z0t)
// mkmat eCMG1c, matrix(z0t)
// mkmat eCMG2c, matrix(z0t)
// mkmat eCMG1cc, matrix(z0t)
mkmat eCMG2cc, matrix(z0t)

***
*** Country-specific estimates for beta and delta: CONSTANT ONLY CASE
***

mata:
// determining the start and end values for each country i (includes missing values)
z0tm	= st_matrix("z0t")
for (i=1; i<=118; i++) {
	robsm = st_matrix("robs")
	robsm = excludemissing(robsm)
	if (i<2) {
		rstart = 1
		rend = robsm[1,1]
	}
	else {
		robsm2 = robsm[1::(i-1),1]
		rstart = colsum(robsm2)+1
		robsm3 = robsm[1::(i),1]
		rend = colsum(robsm3)
	}

	z0ti = z0tm[rstart::rend,1]

// partial demeaning to get rid of the constant
	L = rows(z0ti)
	zt = J(L, 1, 0)
	for (l=1; l<=L; l++) {
		zt[l] = z0ti[l] - (sum(z0ti[1::l])/l)
	}
	zt=zt[2..L]
	
// estimation
mata:
	T = rows(zt)
	b = ceil(sqrt(T))+ 1
	ones = J(T, 1, 1)
	xdtn = log(runningsum(ones))
	ydtn = log(runningsum(zt):^2)
	ydtn = ydtn:-ydtn[1]
// getting rid of missing data
	datamat = zt , ydtn , xdtn 
	datamat = excludemissing(datamat )
	ydtn = datamat[.,2]
	xdtn = datamat[.,3]
	bmcon = invsym(xdtn'*xdtn)*(xdtn'*ydtn)
	betas[i] = bmcon
	deltas[i] = (bmcon-1)/2
}
st_matrix("betas",betas)
st_matrix("deltas",deltas)
end

svmat betas, name(b)
svmat deltas, name(d)


***
*** Create matrix for countries included
***
sort nwbcode year
mat summ_cM=0
mat summ_cMD=0
set seed 1102012


***
*** Start the subsampling
***
drop c csum list
set more off
forvalues l=1/107{
	qui run "$out\cosummhelp118.do"
 display in gr _skip(1) `l'-1+1 _skip(1) _c
}
mat drop z0t robs


***
*** Results
***

* Subsampling: differences in beta --- summ_cM summ_cMD summ_ctM summ_ctMD
mat summ_cM = summ_cM[2..107,1]
mat summ_cMD = summ_cMD[2..107,1]
svmat summ_cM 
svmat summ_cMD


* Mean and Median beta-differences in the N=105 sample
egen betas1m=mean(b1)
egen betas1md=median(b1)


*** constant only case
* Mean and Median Estimates
sum b1, de
scalar b1m=r(mean)
scalar b1md=r(p50)
gen bmco=log(12)*(abs((summ_cM1)-(betas1m)))
sort bmco
cumul bmco, gen(cum_cM1)
gen confic=abs(cum_cM1-0.95)
sort confic
scalar b1ci=bmco[1]
drop confic bmco
* b1m and b1ci
gen bmco=log(12)*(abs((summ_cMD1)-(betas1md)))
sort bmco
cumul bmco, gen(cum_cMD1)
gen confic=abs(cum_cMD1-0.95)
sort confic
scalar b1dci=bmco[1]
drop confic bmco
* b1md and b1dci

** CI for mean (constant)
scalar lowM	=b1m-((1/log(118))*b1ci)
scalar upM	=b1m+((1/log(118))*b1ci)

** CI for median (constant)
scalar lowMD	=b1md-((1/log(118))*b1dci)
scalar upMD		=b1md+((1/log(118))*b1dci)


display in gr _newline _newline _newline _skip(1) "Mean Co-Summability (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Mean" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowM-1)/2  _col(56)  (b1m-1)/2    _col(71)  (upM-1)/2 _continue ///
	in gr _newline _newline _newline _skip(1) "Median Co-Summability (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Median" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMD-1)/2  _col(56)  (b1md-1)/2    _col(71)  (upMD-1)/2 _continue ///
	in gr _newline _newline _newline

sum b1, de



