******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			             Static Nonlinear Models
*                               Table 2				
******************************************************************************
*			   	       Created: 17th November 2014					
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
gen lldebts2=lldebts^2
gen l2ly=l2.ly 
gen l2lk=l2.lk 
gen l2ldebts=l2.ldebts 
gen l2ldebts2=l2ldebts^2
gen l3ly=l3.ly 
gen l3lk=l3.lk 
gen l3ldebts=l3.ldebts 
gen l3ldebts2=l3ldebts^2
sort nwbcode year
gen gfcf2=inv/gdp
gen Debt2p=100*(Debts/gdp)
gen Debt2=(Debts/gdp)


* N=118 IMF sample
tsset nwbcode year
qui xtreg d.ly l.ly l.lk l.ldebts d.lk d.ldebts, fe
sort year nwbcode
by year: egen llyT1=mean(lly) if e(sample)
by year: egen llyT=mean(llyT1)
by year: egen llkT1=mean(llk) if e(sample)
by year: egen llkT=mean(llkT1)
by year: egen lldebtsT1=mean(lldebts) if e(sample)
by year: egen lldebtsT=mean(lldebtsT1)
by year: egen lldebts2T1=mean(lldebts2) if e(sample)
by year: egen lldebts2T=mean(lldebts2T1)
by year: egen l2lyT1=mean(l2ly) if e(sample)
by year: egen l2lyT=mean(l2lyT1)
by year: egen l2lkT1=mean(l2lk) if e(sample)
by year: egen l2lkT=mean(l2lkT1)
by year: egen l2ldebtsT1=mean(l2ldebts) if e(sample)
by year: egen l2ldebtsT=mean(l2ldebtsT1)
by year: egen l2ldebts2T1=mean(l2ldebts2) if e(sample)
by year: egen l2ldebts2T=mean(l2ldebts2T1)
by year: egen l3lyT1=mean(l3ly) if e(sample)
by year: egen l3lyT=mean(l3lyT1)
by year: egen l3lkT1=mean(l3lk) if e(sample)
by year: egen l3lkT=mean(l3lkT1)
by year: egen l3ldebtsT1=mean(l3ldebts) if e(sample)
by year: egen l3ldebtsT=mean(l3ldebtsT1)
by year: egen l3ldebts2T1=mean(l3ldebts2) if e(sample)
by year: egen l3ldebts2T=mean(l3ldebts2T1)

drop llyT1 llkT1 lldebtsT1 lldebts2T1 l2lyT1 l2lkT1 l2ldebtsT1 l2ldebts2T1 l3ldebts2T1 l3ldebtsT1 l3lkT1 l3lyT1

***
*** Static Linear Models
***

run "$out\outreg.ado"
run "$do\xtmg.ado"
cd "$output"


**** Nonsense regression
reg lly llk lldebts lldebts2 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	replace

****
**** 2FE - linear debt
****
xtreg lly llk lldebts i.year, robust fe 
scalar rmse=e(rmse)
scalar tsig=-999
predict e2FE if e(sample), e
xtcd e2FE 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
qui xtreg lly llk lldebts i.year, robust fe 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all


****
**** 2FE - squared debt
****
xtreg lly llk lldebts lldebts2 i.year, robust fe 
scalar rmse=e(rmse)
scalar tsig=-999
predict e2FE2 if e(sample), e
xtcd e2FE2 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
qui xtreg lly llk lldebts lldebts2 i.year, robust fe 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all

****
**** CCEP - linear
****
set seed 1102012
sort nwbcode year
xi: reg lly llk lldebts   ///
	i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT i.nwbcode 
scalar tsig=-999
scalar rmse=e(rmse)
predict eCCEP if e(sample), res
xtcd eCCEP 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
set seed 1102012
xi: reg lly llk lldebts   ///
	i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT i.nwbcode,  vce(bootstrap)
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all


****
**** CCEP - squared debt
****
sort nwbcode year
xi: reg lly llk lldebts lldebts2  ///
	i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT  i.nwbcode|lldebts2T i.nwbcode 
scalar tsig=-999
scalar rmse=e(rmse)
predict eCCEP2 if e(sample), res
xtcd eCCEP2 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xi: reg lly llk lldebts lldebts2  ///
	i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT  i.nwbcode|lldebts2T i.nwbcode,  vce(bootstrap)
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all


****
**** Linear MG with trend
****
xtmg lly llk lldebts , robust trend res(eMGt) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=e(trend_sig)
xtcd eMGt
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts, robust trend  
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: btrend	= testbeta[1..118,3]#J(53,1,1)
mata: bcons 	= testbeta[1..118,4]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_trend=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_trend",btrend)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_trend b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_trend b_cons  ///
 	pop lk debts Debt2p gdppc gdp ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_MGt.dta", replace
restore
drop b_*
scalar drop _all


****
**** MG squared debt and trend
****
xtmg lly llk lldebts lldebts2 , robust trend res(eMGt2) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=e(trend_sig)
xtcd eMGt2
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts lldebts2 , robust trend  
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2	= testbeta[1..118,3]#J(53,1,1)
mata: btrend	= testbeta[1..118,4]#J(53,1,1)
mata: bcons 	= testbeta[1..118,5]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_ldebt2=0
gen b_trend=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_ldebt2",bldebt2)
mata: st_store(.,"b_trend",btrend)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_ldebt2 b_trend b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_ldebt2 b_trend b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_MGt2.dta", replace
restore
drop b_*
scalar drop _all

****
**** Linear CMG 
****
xtmg lly llk lldebts , cce robust res(eCMG) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMG
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts , cce robust 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bcons 	= testbeta[1..118,6]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMG.dta", replace
restore
drop b_*
scalar drop _all


****
**** CMG squared debt 
****
xtmg lly llk lldebts lldebts2 , robust cce res(eCMG2) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMG2
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts lldebts2 , robust cce  
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2	= testbeta[1..118,3]#J(53,1,1)
mata: bcons 	= testbeta[1..118,8]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_ldebt2=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_ldebt2",bldebt2)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_ldebt2 b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_ldebt2 b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMG2.dta", replace
restore
drop b_*
scalar drop _all



****
**** Linear CMG with one additional lag of CA 
****
xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT , robust res(eCMGc) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMGc
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT , robust  
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bcons 	= testbeta[1..118,9]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMGc.dta", replace
restore
drop b_*
scalar drop _all


****
**** CMG with one additional lag: squared debt 
****
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT l2ldebtsT l2ldebts2T, robust res(eCMG2c) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMG2c
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT l2ldebtsT l2ldebts2T, robust 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2	= testbeta[1..118,3]#J(53,1,1)
mata: bcons 	= testbeta[1..118,12]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_ldebt2=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_ldebt2",bldebt2)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_ldebt2 b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_ldebt2 b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMG2c.dta", replace
restore
drop b_*
scalar drop _all


****
**** Linear CMG with 2 additional lags of CA 
****
xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT l3lyT l3lkT l3ldebtsT, robust res(eCMGcc) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMGcc
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT l3lyT l3lkT l3ldebtsT, robust
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bcons 	= testbeta[1..118,12]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_cons ///
 	pop lk debts Debt2p gdppc gdp  
sort nwbcode 
collapse (mean) b_lk b_ldebt b_cons  ///
 	pop lk debts Debt2p gdppc gdp  ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMGcc.dta", replace
restore
drop b_*
scalar drop _all


****
**** CMG with twp additional lag: squared debt 
****
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT l2ldebtsT l2ldebts2T l3lyT l3lkT l3ldebtsT l3ldebts2T, robust res(eCMG2cc) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
xtcd eCMG2cc
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT l2ldebtsT l2ldebts2T l3lyT l3lkT l3ldebtsT l3ldebts2T, robust 
outreg using nonlineardebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(4,3,3,3) append
scalar drop _all
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blk 		= testbeta[1..118,1]#J(53,1,1)
mata: bldebt	= testbeta[1..118,2]#J(53,1,1)
mata: bldebt2	= testbeta[1..118,3]#J(53,1,1)
mata: bcons 	= testbeta[1..118,16]#J(53,1,1)
gen b_lk=0
gen b_ldebt=0
gen b_ldebt2=0
gen b_cons=0
mata: st_store(.,"b_lk",blk)
mata: st_store(.,"b_ldebt",bldebt)
mata: st_store(.,"b_ldebt2",bldebt2)
mata: st_store(.,"b_cons",bcons)
preserve
keep nwbcode b_lk b_ldebt b_ldebt2 b_cons ///
 	pop lk debts Debt2p gdppc gdp  dly 
sort nwbcode 
collapse (mean) b_lk b_ldebt b_ldebt2 b_cons  ///
 	pop lk debts Debt2p gdppc gdp dly ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-17Static_CMG2cc.dta", replace
restore
drop b_*
scalar drop _all

exit

****
**** Residual PURTs
**** 

multipurt e2FE e2FE2 , lags(3)
multipurt eCCEP eCCEP2 , lags(3)
multipurt eMGt eMGt2 , lags(3)
multipurt eCMG eCMG2 , lags(3)
multipurt eCMGc eCMG2c , lags(3)
multipurt eCMGcc eCMG2cc , lags(3)

/*
run "$out\comparison60new.do"
run "$out\comparison90new.do"

xtcd eCMGcc if comparison60==1
xtcd eCMG2cc if comparison60==1
xtcd eCMGcc if comparison90==1
xtcd eCMG2cc if comparison90==1

multipurt eCMGcc eCMG2cc if comparison60==1, lags(3)
multipurt eCMGcc eCMG2cc if comparison90==1, lags(3)

xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT l3lyT l3lkT ///
	l3ldebtsT if comparison60==1, robust 
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT ///
	l2ldebtsT l2ldebts2T l3lyT l3lkT l3ldebtsT l3ldebts2T if comparison60==1, robust
 
xtmg lly llk lldebts llyT llkT lldebtsT l2lyT l2lkT l2ldebtsT l3lyT l3lkT ///
	l3ldebtsT if comparison90==1, robust 
xtmg lly llk lldebts lldebts2 llyT llkT lldebtsT lldebts2T l2lyT l2lkT ///
	l2ldebtsT l2ldebts2T l3lyT l3lkT l3ldebtsT l3ldebts2T if comparison90==1, robust


*/

****
**** Country results: Post-estimation coefficient analysis
**** 

* Standard CMG results
*use "$data\2014-11-17Static_CMGcc.dta", clear
use "$data\2014-11-17Static_CMG2cc.dta", clear

order nwbcode wbcode b_*
gen ltotdebt=ln(Debt2p)
gen lgdppc=ln(gdppc)
gen ldebtpeak=ln(debt_peak)
gen ldebts=log(debts)
gen ltotdebt_sq=ltotdebt^2

* Creates dummies for countries with substantial number of obs in both regimes
run "$out\comparison60new.do"
run "$out\comparison90new.do"

list wbcode b_ldebt b_ldebt2 gdppc, noobs separator(0)

rreg b_lk if comparison60==1
rreg b_ldebt if comparison60==1
rreg b_ldebt2 if comparison60==1

rreg b_lk if comparison90==1
rreg b_ldebt if comparison90==1
rreg b_ldebt2 if comparison90==1


* Graph: rich versus poor
sum gdppc, de
* Cutoff 75%ile of GDP pc: 8828.031

rreg b_ldebt if gdppc>=8828.032
scalar levels=_b[_cons]
rreg b_ldebt2 if gdppc>=8828.032
scalar square=_b[_cons]
nlcom exp(-levels/(square*2))

display 2.77e+07/25146.31
display  638.9038/2299.875

* Graph: box plot for various groupings
list wbcode b_ldebt b_ldebt2 gdppc, noobs separator(0)
gen minmax=exp(-b_ldebt/(b_ldebt2*2))
gen minmaxp=100*(minmax/gdppc)
gen shape=""
replace shape="concave" if b_ldebt>0 & b_ldebt2<0
replace shape="convex" if b_ldebt<0 & b_ldebt2>0
replace shape="negative" if b_ldebt<0 & b_ldebt2<0

list wbcode gdppc b_ldebt b_ldebt2 shape  minmaxp debt_peak, noobs separator(0)
gen location=""
replace location="beyond" if b_ldebt>0 & b_ldebt2<0 & minmaxp<debt_peak
replace location="below" if b_ldebt>0 & b_ldebt2<0 & minmaxp>debt_peak
replace location="beyond" if b_ldebt<0 & b_ldebt2>0 & minmaxp<debt_peak
replace location="below" if b_ldebt<0 & b_ldebt2>0 & minmaxp>debt_peak

graph box lgdppc, over(location) over(shape)

graph box dly, over(location) over(shape)

graph box Debt2p, over(location) over(shape) noout



 