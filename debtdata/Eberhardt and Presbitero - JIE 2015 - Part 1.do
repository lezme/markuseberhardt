******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			Dynamic Linear Models (ECM) and Coefficient Analysis
*                         Table 1 and Figure 3 				
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

* The following commands need to be installed:
* xtmg
* outreg

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
*** Dynamic Linear Models
***

cd "$output"
run "$out\outreg.ado"
run "$do\xtmg.ado"

****
**** 2FE - Table 1, Model [1]
****
qui xtreg dly lly llk lldebts dlk  dldebts i.year, robust fe 
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly]))  , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtb=a[1,2]
scalar debtse=sqrt(b[2,2])
scalar tsig=-999
qui xtreg dly lly llk lldebts dlk  dldebts i.year, robust fe 
scalar rmse=e(rmse)
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtb, "se-debt", debtse) ///
	adec(4,4,4,4,4) replace
predict e2FE if e(sample), e
xtcd e2FE 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3)  nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all


****
**** CCEP - Table 1, Model [2] 
****
set seed 1102012
sort nwbcode year
xi: reg dly lly llk lldebts  dlk  dldebts   ///
	i.nwbcode|dlyT i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT ///
	i.nwbcode|dlkT i.nwbcode|dldebtsT i.nwbcode,  vce(bootstrap)
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly])) , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtb=a[1,2]
scalar debtse=sqrt(b[2,2])
scalar tsig=-999
set seed 1102012
xi: reg dly lly llk lldebts  dlk  dldebts   ///
	i.nwbcode|dlyT i.nwbcode|llyT i.nwbcode|llkT i.nwbcode|lldebtsT ///
	i.nwbcode|dlkT i.nwbcode|dldebtsT i.nwbcode,  vce(bootstrap)
scalar rmse=e(rmse)
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtb, "se-debt", debtse) ///
	adec(4,4,4,4,4) append
predict eCCEP if e(sample), res
xtcd eCCEP 
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3) nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all


****
**** MG with trend  - Table 1, Model [3]
****
xtmg dly lly llk lldebts dlk  dldebts, robust trend res(eMGt) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=e(trend_sig)
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly])) , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtbb=a[1,2]
scalar debtse=sqrt(b[2,2])
xtmg dly lly llk lldebts  dlk  dldebts, robust trend 
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtbb , "se-debt", debtse) ///
	adec(4,4,4,4,4) append
xtcd eMG
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3)  nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all
mat drop abscorr cdtest
gen obs=1 if !missing(eMG)
drop eMG
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blly	 	= testbeta[1..118,1]#J(53,1,1)
mata: bllk 		= testbeta[1..118,2]#J(53,1,1)
mata: blldebt	= testbeta[1..118,3]#J(53,1,1)
mata: bdlk 		= testbeta[1..118,4]#J(53,1,1)
mata: bdldebt	= testbeta[1..118,5]#J(53,1,1)
mata: btrend	= testbeta[1..118,6]#J(53,1,1)
mata: bcons 	= testbeta[1..118,7]#J(53,1,1)
gen b_lly=0
gen b_llk=0
gen b_lldebt=0
gen b_dlk=0
gen b_dldebt=0
gen b_trend=0
gen b_cons=0
mata: st_store(.,"b_lly",blly)
mata: st_store(.,"b_llk",bllk)
mata: st_store(.,"b_lldebt",blldebt)
mata: st_store(.,"b_dldebt",bdldebt)
mata: st_store(.,"b_trend",btrend)
mata: st_store(.,"b_cons",bcons)
qui gen feMG=dly-(b_lly*lly+b_llk*llk+b_lldebt*lldebts+b_dlk*dlk+b_dldebt*dldebts+b_cons) if e(sample)
qui gen uMG=ly-((-b_llk/b_lly)*lk+(-b_lldebt/b_lly)*ldebts) if e(sample)
preserve
keep nwbcode obs b_lly b_llk b_lldebt b_dlk b_dldebt  b_trend b_cons ///
 	pop lk debts Debt Debt2p gdppc gdp inflation
sort nwbcode 
collapse (sum) obs (mean) b_lly b_llk b_lldebt b_dlk b_dldebt   b_trend b_cons  ///
 	pop lk debts Debt2p gdppc gdp inflation ///
	(max)  debt_peak=Debt2p, by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-10ECM_MGt.dta", replace
restore
drop b_lly b_llk b_lldebt b_dlk b_dldebt  b_trend b_cons
scalar drop _all



****
**** CMG (standard)  - Table 1, Model [4]
**** 
xtmg dly lly llk lldebts  dlk  dldebts  ///
	dlyT llyT llkT lldebtsT dlkT dldebtsT , robust res(eCMG) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly]))  , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtb=a[1,2]
scalar debtse=sqrt(b[2,2])
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT , robust  
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtb, "se-debt", debtse) ///
	adec(4,4,4,4,4) append
xtcd eCMG
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3)  nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all
mat drop abscorr cdtest
drop eCMG
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blly	 	= testbeta[1..118,1]#J(53,1,1)
mata: bllk 		= testbeta[1..118,2]#J(53,1,1)
mata: blldebt	= testbeta[1..118,3]#J(53,1,1)
mata: bdlk 		= testbeta[1..118,4]#J(53,1,1)
mata: bdldebt	= testbeta[1..118,5]#J(53,1,1)
mata: bcons 	= testbeta[1..118,12]#J(53,1,1)
gen bc_lly=0
gen bc_llk=0
gen bc_lldebt=0
gen bc_dlk=0
gen bc_dldebt=0
gen bc_cons=0
mata: st_store(.,"bc_lly",blly)
mata: st_store(.,"bc_llk",bllk)
mata: st_store(.,"bc_lldebt",blldebt)
mata: st_store(.,"bc_dlk",bdlk)
mata: st_store(.,"bc_dldebt",bdldebt)
mata: st_store(.,"bc_cons",bcons)
qui gen feCMG=dly-(bc_lly*lly+bc_llk*llk+bc_lldebt*lldebts+bc_dlk*dlk+bc_dldebt*dldebts+bc_cons) if e(sample)
qui gen uCMG=ly-((-bc_llk/bc_lly)*lk+(-bc_lldebt/bc_lly)*ldebts) if e(sample)
preserve
keep nwbcode bc_lly bc_llk bc_lldebt bc_dlk bc_dldebt  bc_cons
sort nwbcode 
collapse (mean) bc_lly bc_llk bc_lldebt bc_dlk bc_dldebt  bc_cons , by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-10ECM_CMG1.dta", replace
restore
drop bc_lly bc_llk bc_lldebt bc_dlk bc_dldebt  bc_cons 
scalar drop _all



****
**** CMG with one additional lag  - Table 1, Model [5]
**** 
xtmg dly lly llk lldebts  dlk  dldebts  ///
	dlyT llyT llkT lldebtsT dlkT dldebtsT ldlyT ldlkT ldldebtsT, robust res(eCMGlag2) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly]))  , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtb=a[1,2]
scalar debtse=sqrt(b[2,2])
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT ldlyT ldlkT ldldebtsT, robust  
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtb, "se-debt", debtse) ///
	adec(4,4,4,4,4) append
xtcd eCMG
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3)  nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all
mat drop abscorr cdtest
drop eCMGlag2
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blly	 	= testbeta[1..118,1]#J(53,1,1)
mata: bllk 		= testbeta[1..118,2]#J(53,1,1)
mata: blldebt	= testbeta[1..118,3]#J(53,1,1)
mata: bdlk 		= testbeta[1..118,4]#J(53,1,1)
mata: bdldebt	= testbeta[1..118,5]#J(53,1,1)
mata: bcons 	= testbeta[1..118,15]#J(53,1,1)
gen bc2_lly=0
gen bc2_llk=0
gen bc2_lldebt=0
gen bc2_dlk=0
gen bc2_dldebt=0
gen bc2_cons=0
mata: st_store(.,"bc2_lly",blly)
mata: st_store(.,"bc2_llk",bllk)
mata: st_store(.,"bc2_lldebt",blldebt)
mata: st_store(.,"bc2_dlk",bdlk)
mata: st_store(.,"bc2_dldebt",bdldebt)
mata: st_store(.,"bc2_cons",bcons)
qui gen feCMG2=dly-(bc2_lly*lly+bc2_llk*llk+bc2_lldebt*lldebts+bc2_dlk*dlk+bc2_dldebt*dldebts+bc2_cons) if e(sample)
qui gen uCMG2=ly-((-bc2_llk/bc2_lly)*lk+(-bc2_lldebt/bc2_lly)*ldebts) if e(sample)
preserve
keep nwbcode bc2_lly bc2_llk bc2_lldebt bc2_dlk bc2_dldebt  bc2_cons
sort nwbcode 
collapse (mean) bc2_lly bc2_llk bc2_lldebt bc2_dlk bc2_dldebt  bc2_cons , by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-10ECM_CMG2.dta", replace
restore
drop bc2_lly bc2_llk bc2_lldebt bc2_dlk bc2_dldebt  bc2_cons 
scalar drop _all


****
**** CMG with two additional lags  - Table 1, Model [6]
**** 
xtmg dly lly llk lldebts  dlk  dldebts  ///
	dlyT llyT llkT lldebtsT dlkT dldebtsT ldlyT ldlkT ldldebtsT l2dlyT l2dlkT l2dldebtsT, robust res(eCMGlag3) 
mat testb=e(betas)
scalar mse=e(sigma)
scalar rmse=sqrt(mse)
scalar tsig=-999
nlcom (_b[llk]/(-_b[lly]))  (_b[lldebts]/(-_b[lly]))  , post 
mat a=e(b)
mat b=e(V)
scalar capb=a[1,1]
scalar capse=sqrt(b[1,1])
scalar debtb=a[1,2]
scalar debtse=sqrt(b[2,2])
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT ldlyT ldlkT ldldebtsT l2dlyT l2dlkT l2dldebtsT, robust  
outreg using heterogdebt.csv, bdec(3) se nolabel br 3aster comma ///
	addstat("RMSE", rmse, "lr-cap", capb, "se-cap", capse, ///
	"lr-debt", debtb, "se-debt", debtse) ///
	adec(4,4,4,4,4) append
xtcd eCMGlag3
mat cdtest=r(pesaran)
scalar cd=cdtest[1,1]
mat abscorr=r(abscorr)
scalar corr=abscorr[1,1]
outreg lly using heterogdebt.csv, bdec(3)  nolabel br 3aster comma ///
	addstat("trends", tsig, "CD test", cd, "Abs Corr", corr) ///
	adec(2,4,4) append
scalar drop _all
mat drop abscorr cdtest
drop eCMG
sort wbcode year
mata: testbeta	= st_matrix("testb")
mata: blly	 	= testbeta[1..118,1]#J(53,1,1)
mata: bllk 		= testbeta[1..118,2]#J(53,1,1)
mata: blldebt	= testbeta[1..118,3]#J(53,1,1)
mata: bdlk 		= testbeta[1..118,4]#J(53,1,1)
mata: bdldebt	= testbeta[1..118,5]#J(53,1,1)
mata: bcons 	= testbeta[1..118,18]#J(53,1,1)
gen bc3_lly=0
gen bc3_llk=0
gen bc3_lldebt=0
gen bc3_dlk=0
gen bc3_dldebt=0
gen bc3_cons=0
mata: st_store(.,"bc3_lly",blly)
mata: st_store(.,"bc3_llk",bllk)
mata: st_store(.,"bc3_lldebt",blldebt)
mata: st_store(.,"bc3_dlk",bdlk)
mata: st_store(.,"bc3_dldebt",bdldebt)
mata: st_store(.,"bc3_cons",bcons)
qui gen feCMG3=dly-(bc3_lly*lly+bc3_llk*llk+bc3_lldebt*lldebts+bc3_dlk*dlk+bc3_dldebt*dldebts+bc3_cons) if e(sample)
qui gen uCMG3=ly-((-bc3_llk/bc3_lly)*lk+(-bc3_lldebt/bc3_lly)*ldebts) if e(sample)
preserve
keep nwbcode bc3_lly bc3_llk bc3_lldebt bc3_dlk bc3_dldebt  bc3_cons
sort nwbcode 
collapse (mean) bc3_lly bc3_llk bc3_lldebt bc3_dlk bc3_dldebt  bc3_cons , by(nwbcode)
decode nwbcode, gen(wbcode)
sort wbcode
save "$data\2014-11-10ECM_CMG3.dta", replace
restore
drop bc3_lly bc3_llk bc3_lldebt bc3_dlk bc3_dldebt  bc3_cons 
scalar drop _all


exit

****
**** t-bar statistics
**** 

* MG - Table 1, Model [3]
xtmg dly lly llk lldebts  dlk  dldebts,  trend
mat tMG=e(tbetas)
mat tMG=tMG[1..118,1]
svmat tMG
reg tMG1
drop tMG1

* CMG - Table 1, Model [4]
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT   
mat tCMG=e(tbetas)
mat tCMG=tCMG[1..118,1]
svmat tCMG
reg tCMG1
drop tCMG1
	
* CMG - Table 1, Model [5]
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT  ldlyT	ldlkT ldldebtsT
mat tCMG=e(tbetas)
mat tCMG=tCMG[1..118,1]
svmat tCMG
reg tCMG1
drop tCMG1
	
* CMG - Table 1, Model [6]
xtmg dly lly llk lldebts  dlk  dldebts   ///
	dlyT llyT llkT lldebtsT  dlkT dldebtsT  ldlyT	ldlkT ldldebtsT l2dlyT	l2dlkT l2dldebtsT
mat tCMG=e(tbetas)
mat tCMG=tCMG[1..118,1]
svmat tCMG
reg tCMG1
drop tCMG1
	
drop _Inwbcode* _Inwb*

****
**** L-R residual and Factor diagnostics
**** 
rename feCMG2 feCMG2_
rename feCMG3 feCMG3_
* Factors
xtcd feMG feCMG feCMGt feCMG2_ feCMG3_ 
* LR model residuals
xtcd feCMG2t feCMG3t

rename uCMG2 uCMG2_
rename uCMG3 uCMG3_
* Factors
xtcd uMG uCMG uCMGt uCMG2_ uCMG3_
* LR model residuals
xtcd uCMG2t uCMG3t


****
**** Post-estimation coefficient analysis
**** 

use "$data\2014-11-10ECM_MGt.dta", clear
merge wbcode using "$data\2014-11-10ECM_CMG1.dta" "$data\2014-11-10ECM_CMG2.dta" "$data\2014-11-10ECM_CMG3.dta"  ,  sort
drop _merge*

order nwbcode wbcode obs b_* bc_*  bc2_* bc3_*
gen ltotdebt=ln(Debt2p)
gen lgdppc=ln(gdppc)
gen ldebtpeak=ln(debt_peak)
gen ldebts=log(debts)
gen ltotdebt_sq=ltotdebt^2
gen ltotdebt_cu=ltotdebt^3

* MG
gen lr_debt	= b_lldebt/-b_lly
rreg lr_debt
drop lr_debt

* CMG
gen lr_debt	= bc_lldebt/-bc_lly
rreg lr_debt
drop lr_debt

* CMG one additional lag
gen lr_debt	= bc2_lldebt/-bc2_lly
rreg lr_debt
drop lr_debt

* CMG two additional lags
gen lr_debt	= bc3_lldebt/-bc3_lly
rreg lr_debt
drop lr_debt



*Capital
gen lr_k	= b_llk/-b_lly
rreg lr_k
drop lr_k

gen lr_k	= bc_llk/-bc_lly
rreg lr_k
drop lr_k

gen lr_k	= bc2_llk/-bc2_lly
rreg lr_k
drop lr_k

gen lr_k	= bc3_llk/-bc3_lly
rreg lr_k
drop lr_k





***
*** Long-run debt analysis
***

gen lr_k	= bc3_llk/-bc3_lly
gen lr_debt	= bc3_lldebt/-bc3_lly
gen sr_k	= bc_dlk
gen sr_debt	= bc_dldebt

gen lr_debt0	= b_lldebt/-b_lly
gen lr_debt1	= bc_lldebt/-bc_lly
gen lr_debt3	= bc2_lldebt/-bc2_lly
gen lr_debt4	= bc3_lldebt/-bc3_lly

format lr_debt %2.1f

rreg lr_debt4, genwt(w_4)
sum w_4, de
* w_4>.174 (5% cut)
rreg lr_debt1, genwt(w_1)
sum w_1, de
* w_1>0 (<10% cut)
rreg lr_debt3, genwt(w_3)
sum w_3, de
* w_3>0 (<10% cut)

tab wbcode if w_3==0 | w_1==0 | w_4<=.174
/*
        BHS |          1       14.29       14.29
        OMN |          1       14.29       28.57
        SGP |          1       14.29       42.86
        TTO |          1       14.29       57.14
        TZA |          1       14.29       71.43
        ZAF |          1       14.29       85.71
        ZAR |          1       14.29      100.00
*/
***********************************
*** FIGURE 4 IN THE FINAL DRAFT ***
***********************************





***
*** 4a, 4b
***

*******************************************************************************************************************
*** Following two plots are included in the paper using CMG-3 results (Graphs a & b in the final draft)
*******************************************************************************************************************
twoway	(fpfitci lr_debt4 ltotdebt if w_4>.174, level(90)) ///
		(scatter lr_debt4 ltotdebt if w_4>.174, msymbol(none) mlab(wbcode) mlabposition(0)) ///
, legend(off) ytitle("Long-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25)) ///
	 scheme(s2mono) title("(a) Debt/GDP ratio and long-run debt impact") subtitle("CMG model with two additional lags") ///
	yline(0, lwidth(.25)) 
graph export "$graphfolder\LRdebtratio1_IMF.png", as(png) replace 
*******************************************************************************************************************
twoway	(fpfitci lr_debt4 ldebtpeak if w_4>.174, level(90)) ///
		(scatter lr_debt4 ldebtpeak if w_4>.174, msymbol(none) mlab(wbcode) mlabposition(0)) ///
, legend(off) ytitle("Long-run debt coefficient") xtitle("Peak Debt/GDP (in logs)") ///
	 scheme(s2mono) title("(b) Peak Debt/GDP ratio and long-run debt impact") subtitle("CMG model with two additional lags") ///
	 xline(4.5, lwidth(.25))  yline(0, lwidth(.25))
graph export "$graphfolder\LRdebtpeak1_IMF.png", as(png) replace 
		*******************************************************************************************************************
		*** Dropping countries where net and gross debt deviates a lot
		*******************************************************************************************************************
		twoway	(fpfitci lr_debt4 ltotdebt if w_4>.174 & wbcode!="JPN" & wbcode!="KOR" & ///
					wbcode!="SWE" & wbcode!="FIN" & wbcode!="ISL" & wbcode!="NOR" & wbcode!="LUX", level(90)) ///
				(scatter lr_debt4 ltotdebt if w_4>.174 & wbcode!="JPN" & wbcode!="KOR" & ///
				wbcode!="SWE" & wbcode!="FIN" & wbcode!="ISL" & wbcode!="NOR" & wbcode!="LUX", msymbol(none) mlab(wbcode) mlabposition(0)) ///
		, legend(off) ytitle("Long-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25)) ///
			 scheme(s2mono) title("(a) Debt/GDP ratio and long-run debt impact") subtitle("CMG model with two additional lags") ///
			yline(0, lwidth(.25)) 
		graph export "$graphfolder\LRdebtratio1_IMF_reduced.png", as(png) replace 
		*******************************************************************************************************************
	
	
***
*** 4c
***

*******************************************************************************************************************
*** Box plots (Graph c in the final draft)
*******************************************************************************************************************
centile ltotdebt, centile(20 40 60 80)
gen debtgroup=0
replace debtgroup=1 if ltotdebt<3.558634
replace debtgroup=2 if ltotdebt>3.558634 & ltotdebt<3.780657
replace debtgroup=3 if ltotdebt>3.780657 & ltotdebt<4.014937
replace debtgroup=4 if ltotdebt>4.014937 & ltotdebt<4.365157
replace debtgroup=5 if ltotdebt>4.365157
* Cutoffs:
*      Debt2p |     118         20      35.11822        27.04121    38.51785
*             |                 40      43.84554        39.61153    47.70243
*             |                 60      55.42028        48.44176    60.78948
*             |                 80      78.66767        62.16687    87.75612
graph box lr_debt4, over(debtgroup) noout ytitle("Long-run debt coefficient") note("") ///
	title("(c) Debt/GDP ratio (quintiles) and long-run debt impact") ///
	subtitle("CMG model with two additional lags") box(1, color(gray)) medtype(cline) medline(lwidth(.6)) ///
	yline(0, lwidth(.25) lcolor(black)) scheme(s2mono) 
graph export "$graphfolder\LRdebtquint_IMF.png", as(png) replace 
*******************************************************************************************************************

***
*** 4d
***

*******************************************************************************************************************
* Debt levels and impact of debt (by income group):  Graph d in the final draft
*******************************************************************************************************************
sum lgdppc, de

twoway	(fpfitci lr_debt4 ltotdebt if lgdppc>9.085688 & w_4>.174, level(90) lcolor(black) fc(gs14) alp(dash)) ///
		(fpfitci lr_debt4 ltotdebt if lgdppc<9.085688 & w_4>.174, lcolor(red) lpattern(solid) level(90) fc(none) alp(dash)) ///
, legend(off) ytitle("Long-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25))  yline(0, lwidth(.25)) ///
	 scheme(s2mono) title("(d) Debt/GDP ratio and long-run debt impact") subtitle("CMG with two additional lags, red=bottom 75% by average income") ///
	yscale(range(-0.3 0.3)) ylabel(-0.3(0.1)0.3)   
graph export "$graphfolder\LRcoeffRichPoor_IMF.png", as(png) replace 

* Rich countries on their own (27 countries only)
twoway	(fpfitci lr_debt4 ltotdebt if lgdppc>9.085688 & w_4>.174, level(66) lcolor(black) fc(gs14) alp(dash)) ///
		(scatter lr_debt4 ltotdebt if lgdppc>9.085688 & w_4>.174, mlab(wbcode) mlabposition(9))
, legend(off) ytitle("Long-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25))  yline(0, lwidth(.25)) ///
	 scheme(s2mono) title("(d) Debt/GDP ratio and long-run debt impact") subtitle("CMG with two additional lags, rich countries only") 
*******************************************************************************************************************	 


***
*** 4e and f
***

*******************************************************************************************************************
* All three CMG models using fpfit (Graph e in the final draft) 
*******************************************************************************************************************
twoway	///
	(fpfit lr_debt1 ltotdebt if w_1>0, lcolor(blue) lwidth(.7) lpattern(solid)) ///
	(fpfit lr_debt3 ltotdebt if w_3>0, lcolor(red)  lwidth(.7) lpattern(dash)) ///
	(fpfit lr_debt4 ltotdebt if w_4>.174, lcolor(black)  lwidth(.5) lpattern(dash_dot)) ///
, legend(off) ytitle("Long-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25)) ///
	 scheme(s2mono) title("(e) Debt/GDP ratio and long-run debt impact") subtitle("Various CMG models") ///
	yscale(range(-0.2 0.2)) ylabel(-0.2(0.1)0.2)   yline(0, lwidth(.25)) 
graph export "$graphfolder\LRdebtratio_IMFRALL.png", as(png) replace 
*******************************************************************************************************************
* LR with debt peak (Graph f in the final draft)
*******************************************************************************************************************
twoway	///
	(fpfit lr_debt1 ldebtpeak if w_1>0, lcolor(blue) lwidth(.7) lpattern(solid)) ///
	(fpfit lr_debt3 ldebtpeak if w_3>0, lcolor(red)  lwidth(.7) lpattern(dash)) ///
	(fpfit lr_debt4 ldebtpeak if w_4>.174, lcolor(black)  lwidth(.5) lpattern(dash_dot)) ///
, legend(off) ytitle("Long-run debt coefficient") xtitle("Peak Debt/GDP (in logs)") ///
	 scheme(s2mono) title("Peak Debt/GDP ratio and long-run debt impact") subtitle("Various CMG models") ///
	 xline(4.5, lwidth(.25)) yscale(range(-0.2 0.2)) ylabel(-0.2(0.1)0.2) yline(0, lwidth(.25))
graph export "$graphfolder\LRdebtpeak_IMFALL.png", as(png) replace 
*******************************************************************************************************************

gen sr_debt1	= bc_dldebt
gen sr_debt3	= bc2_dldebt
gen sr_debt4	= bc3_dldebt

twoway	///
	(fpfit sr_debt1 ltotdebt if w_1>0, lcolor(blue) lwidth(.7) lpattern(solid))  ///
	(fpfit sr_debt3 ltotdebt if w_3>0, lcolor(red)  lwidth(.7) lpattern(dash)) ///
	(fpfit sr_debt4 ltotdebt if w_4>.174, lcolor(black)  lwidth(.5) lpattern(dash_dot)) ///
, legend(off) ytitle("Short-run debt coefficient") xtitle("Average Debt/GDP (in logs)") xline(4.5, lwidth(.25)) ///
	 scheme(s2mono) title("(f) Debt/GDP ratio and short-run debt impact") subtitle("Various CMG models") ///
	yscale(range(-0.3 0.3)) ylabel(-0.3(0.1)0.3)   yline(0, lwidth(.25)) 
graph export "$graphfolder\SRdebtratio_IMFRALL.png", as(png) replace 

