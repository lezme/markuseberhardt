******************************************************************************
*
*                           Replication files for:
*           "Democracy, Growth, Heterogeneity and Robustness"
*    (conditionally accepted for publication in the European Economic Review)
*
******************************************************************************
*                 Figure 4 - Sample reduction exercises - Panel (a) 
******************************************************************************
*                       Created: 12th November 2018
*                       Revised: 21st April 2022
******************************************************************************
*                          Markus Eberhardt
******************************************************************************
*            School of Economics, University of Nottingham,
*            University Park, Nottingham NG7 2RD, England
*              email: markus.eberhardt@nottingham.ac.uk
*            web: https://sites.google.com/site/medevecon/
******************************************************************************

* Note: The HHK estimation is quite time-consuming due to the bootstrap.

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


* SPMAT command needs to be installed manually
* search spmat 

cap ssc install xtabond2 
cap ssc install xtivreg2 


***
*** ANRR FE regressions (with 4 lags of the dependent variable)
***

use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 
xtreg y l(1/4).y dem yy*, fe r cluster(wbcode2)
nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y]))

gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum if e(sample)

matrix FETable2Column3=J(36,13,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset wbcode2 year
local z "6	7	9	10	12	13	14	15	16	17	18	19	20	21	22	23	25	26	27	28	29	30	31	32	35	36	37	38	39	40	41	42	44	45	46	47"
foreach k of local z{
	qui: xtreg y l(1/4).y dem yy* if csum>=`k', fe r cluster(wbcode2)
	* Min number of observations
	mat FETable2Column3[`i',1]=`k'
	* Total number of observations
	mat FETable2Column3[`i',2]=e(N)
	* Number of Countries
	mat FETable2Column3[`i',3]=e(N_g)
	qui nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y])) (_b[l.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat FETable2Column3[`i',4]=b_2
	mat FETable2Column3[`i',5]=se_2
	* Dem coefficient (long-run)
	mat FETable2Column3[`i',8]=b_1
	mat FETable2Column3[`i',9]=se_1	
	* OLS estimation
	qui reg y l(1/4).y dem yy*  if e(sample), robust
	qui nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y])) (_b[l.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	* Dem coefficient (long-run)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	mat FETable2Column3[`i',12]=b_1
	mat FETable2Column3[`i',13]=se_1	
	mat FETable2Column3[`i',10]=b_2
	mat FETable2Column3[`i',11]=se_2
	di in gr "Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more observations. N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list FETable2Column3
svmat FETable2Column3, names(FE23_)
rename FE23_1 FE23_minobs
rename FE23_2 FE23_obs
rename FE23_3 FE23_N
rename FE23_4 FE23_b_lagy
rename FE23_5 FE23_se_lagy
rename FE23_6 FE23_b_dem
rename FE23_7 FE23_se_dem
rename FE23_8 FE23_b_lr_dem
rename FE23_9 FE23_se_lr_dem
rename FE23_10 FE23_ls_b_lagy
rename FE23_11 FE23_ls_se_lagy
rename FE23_12 FE23_ls_lr_dem
rename FE23_13 FE23_ls_se_lr_dem

preserve
keep FE23_*
gen model = _n if !missing(FE23_minobs)
drop if _n>36
order model
sort FE23_minobs
save  "$otherdata/FE_2_3_reduction_2022.dta", replace
restore


***
*** ANRR AB regressions (with 4 lags of the dependent variable)
***

use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 

* Use AB sample to create available number of observations
xtabond2 y l(1/4).y dem yy* ,  gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan

gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum

matrix ABTable2Column7=J(36,13,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset wbcode2 year
local z "5	6	8	9	11	12	13	14	15	16	17	18	19	20	21	22	24	25	26	27	28	29	30	31	34	35	36	37	38	39	40	41	43	44	45	46"
foreach k of local z{
	qui: xtabond2 y l(1/4).y dem yy*   if csum>=`k',  gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan
	* Min number of observations
	mat ABTable2Column7[`i',1]=`k'
	* Total number of observations
	mat ABTable2Column7[`i',2]=e(N)
	* Number of Countries
	mat ABTable2Column7[`i',3]=e(N_g)
	* Dem coefficient (short-run)
	mat ABTable2Column7[`i',6]=_b[dem]
	mat ABTable2Column7[`i',7]=_se[dem]	
	qui nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y])) (_b[l.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat ABTable2Column7[`i',4]=b_2
	mat ABTable2Column7[`i',5]=se_2
	* Dem coefficient (long-run)
	mat ABTable2Column7[`i',8]=b_1
	mat ABTable2Column7[`i',9]=se_1	
	* OLS estimation
	qui reg y l(1/4).y dem yy*  if e(sample), robust
	qui nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y])) (_b[l.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	* Dem coefficient (long-run)
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	mat ABTable2Column7[`i',12]=b_1
	mat ABTable2Column7[`i',13]=se_1	
	mat ABTable2Column7[`i',10]=b_2
	mat ABTable2Column7[`i',11]=se_2
	di in gr "Table 2(7), Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more observations. N=" in ye e(N) in gr " obs."

	local i = `i'+1
}	

mat list ABTable2Column7
svmat ABTable2Column7, names(AB_2_7_)
rename AB_2_7_1 AB27_minobs
rename AB_2_7_2 AB27_obs
rename AB_2_7_3 AB27_N
rename AB_2_7_4 AB27_b_lagy
rename AB_2_7_5 AB27_se_lagy
rename AB_2_7_6 AB27_b_dem
rename AB_2_7_7 AB27_se_dem
rename AB_2_7_8 AB27_b_lr_dem
rename AB_2_7_9 AB27_se_lr_dem
rename AB_2_7_10 AB27_ls_b_lagy
rename AB_2_7_11 AB27_ls_se_lagy
rename AB_2_7_12 AB27_ls_lr_dem
rename AB_2_7_13 AB27_ls_se_lr_dem

preserve
keep AB27_*
gen model = _n
order model
sort AB27_minobs
drop if _n>36
save  "$otherdata/AB_2_7_reduction_2022.dta", replace
restore


***
*** ANRR HHK regressions (with 4 lags of the dependent variable)
***

* The following is taken verbatim from the code provided by ANRR
********************************************************************************

global limit=25                     /* Evaluate effects 25 years after transition */
local repsBS=100        /* Number of bootstrap repetitions            */

/*******************************************************************************************************************************/
/*******************************************************************************************************************************/
/*********************DEFINE REQUIRED PROGRAMS THAT WILL BE USED DURING THE EXECUTION OF THIS DO FILE **************************/
/*******************************************************************************************************************************/
/*******************************************************************************************************************************/

capture program drop vareffects
program define vareffects, eclass

quietly: nlcom (effect1: _b[shortrun]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post

quietly: nlcom (effect2: _b[effect1]*_b[lag1]+_b[shortrun]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post

quietly: nlcom (effect3: _b[effect2]*_b[lag1]+_b[effect1]*_b[lag2]+_b[shortrun]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
	  
quietly: nlcom (effect4: _b[effect3]*_b[lag1]+_b[effect2]*_b[lag2]+_b[effect1]*_b[lag3]+_b[shortrun]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post	  

forvalues j=5(1)$limit{	  
local j1=`j'-1
local j2=`j'-2
local j3=`j'-3
local j4=`j'-4

quietly: nlcom (effect`j': _b[effect`j1']*_b[lag1]+_b[effect`j2']*_b[lag2]+_b[effect`j3']*_b[lag3]+_b[effect`j4']*_b[lag4]+_b[shortrun]) ///
	  (effect`j1': _b[effect`j1']) ///
	  (effect`j2': _b[effect`j2']) ///
	  (effect`j3': _b[effect`j3']) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post	  	  

}

quietly: nlcom (effect$limit: _b[effect$limit]) ///
	  (longrun: _b[shortrun]/(1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4])) ///
      (shortrun: _b[shortrun]) ///
	  (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  , post
ereturn display
end


capture program drop vareffects8
program define vareffects8, eclass

quietly: nlcom (effect1: _b[shortrun]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post

quietly: nlcom (effect2: _b[effect1]*_b[lag1]+_b[shortrun]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post

quietly: nlcom (effect3: _b[effect2]*_b[lag1]+_b[effect1]*_b[lag2]+_b[shortrun]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post
	  
quietly: nlcom (effect4: _b[effect3]*_b[lag1]+_b[effect2]*_b[lag2]+_b[effect1]*_b[lag3]+_b[shortrun]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post
	  
quietly: nlcom (effect5: _b[effect4]*_b[lag1]+_b[effect3]*_b[lag2]+_b[effect2]*_b[lag3]+_b[effect1]*_b[lag4]+_b[shortrun]) ///
	  (effect4: _b[effect4]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post
	  
quietly: nlcom (effect6: _b[effect5]*_b[lag1]+_b[effect4]*_b[lag2]+_b[effect3]*_b[lag3]+_b[effect2]*_b[lag4]+_b[effect1]*_b[lag5]+_b[shortrun]) ///
	  (effect5: _b[effect5]) ///
	  (effect4: _b[effect4]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post	

quietly: nlcom (effect7: _b[effect6]*_b[lag1]+_b[effect5]*_b[lag2]+_b[effect4]*_b[lag3]+_b[effect3]*_b[lag4]+_b[effect2]*_b[lag5]+_b[effect1]*_b[lag6]+_b[shortrun]) ///
	  (effect6: _b[effect6]) ///
	  (effect5: _b[effect5]) ///
	  (effect4: _b[effect4]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post	

quietly: nlcom (effect8: _b[effect7]*_b[lag1]+_b[effect6]*_b[lag2]+_b[effect5]*_b[lag3]+_b[effect4]*_b[lag4]+_b[effect3]*_b[lag5]+_b[effect2]*_b[lag6]+_b[effect1]*_b[lag7]+_b[shortrun]) ///
	  (effect7: _b[effect7]) ///
	  (effect6: _b[effect6]) ///
	  (effect5: _b[effect5]) ///
	  (effect4: _b[effect4]) ///
	  (effect3: _b[effect3]) ///
	  (effect2: _b[effect2]) ///
	  (effect1: _b[effect1]) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post	
	  
	  
forvalues j=9(1)$limit{	  
local j1=`j'-1
local j2=`j'-2
local j3=`j'-3
local j4=`j'-4
local j5=`j'-5
local j6=`j'-6
local j7=`j'-7
local j8=`j'-8

quietly: nlcom (effect`j': _b[effect`j1']*_b[lag1]+_b[effect`j2']*_b[lag2]+_b[effect`j3']*_b[lag3]+_b[effect`j4']*_b[lag4]+_b[effect`j5']*_b[lag5]+_b[effect`j6']*_b[lag6]+_b[effect`j7']*_b[lag7]+_b[effect`j8']*_b[lag8]+_b[shortrun]) ///
	  (effect`j1': _b[effect`j1']) ///
	  (effect`j2': _b[effect`j2']) ///
	  (effect`j3': _b[effect`j3']) ///
	  (effect`j4': _b[effect`j4']) ///
	  (effect`j5': _b[effect`j5']) ///
	  (effect`j6': _b[effect`j6']) ///
	  (effect`j7': _b[effect`j7']) ///
	  (shortrun: _b[shortrun]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post	  	  

}

quietly: nlcom (effect$limit: _b[effect$limit]) ///
	  (longrun: _b[shortrun]/(1-_b[lag1]-_b[lag2]-_b[lag3]-_b[lag4]-_b[lag5]-_b[lag6]-_b[lag7]-_b[lag8])) ///
      (shortrun: _b[shortrun]) ///
	  (persistence: _b[lag1]+_b[lag2]+_b[lag3]+_b[lag4]+_b[lag5]+_b[lag6]+_b[lag7]+_b[lag8]) ///
	  (lag1: _b[lag1]) ///
	  (lag2: _b[lag2]) ///
	  (lag3: _b[lag3]) ///
	  (lag4: _b[lag4]) ///
	  (lag5: _b[lag5]) ///
	  (lag6: _b[lag6]) ///
	  (lag7: _b[lag7]) ///
	  (lag8: _b[lag8]) ///
	  , post
ereturn display
end


capture program drop helm
program define helm
*
* This program will do Helmert transformation for a list of variables
* NOTE:  must have variables named id, year   
* to use enter >> helm var1 var2...
* new variables will be names with h_ in front h_var1  and so on
*
qui while "`1'"~="" {
gsort id -year                /*sort years descending */
tempvar one sum n m w 
* capture drop h_`1'         /* IF the variable exist - it will remain and not generated again */
gen `one'=1 if `1'~=.             /*generate one if x is nonmissing */
qui by id: gen `sum'=sum(`1')-`1' /*running sum without current element */
qui by id: gen `n'=sum(`one')-1     /*number of obs included in the sum */
replace `n'=. if `n'<=0             /* n=0 for last observation and =-1 if
                                   last observation is missing*/
gen `m'=`sum'/`n'                 /* m is forward mean of variable x*/
gen `w'=sqrt(`n'/(`n'+1))         /* weight on mean difference */
capture gen h_`1'=`w'*(`1'-`m')             /* transformed variable */ 
sort id year
mac shift
}
end

capture program drop hhkBS
capture program hhkBS, eclass
syntax anything[, ydeep(integer 1960) ystart(integer 1964) yfinal(integer 2009) truncate(integer 4) depvarlags(integer 4)]
	local 0 `anything' 
	gettoken yvar 0 : 0 /*dependent variable*/
	gettoken seqexg   0 : 0, match(par) /*Sequentially exogenous variables*/
	gettoken gmminst  0 : 0, match(par) /*gmm style instruments*/
	gettoken gmmtrunc 0 : 0, match(par) /*gmm style instruments, truncated*/
	gettoken excov 0 : 0, match(par) /*exogenous covariates: diff coefficient for each equation*/
		
/*declares panel structure and defines estimation sample*/
quietly: tsset, clear
quietly: tsset newcl year
sort newcl year
quietly: xtreg `yvar' `seqexg' `excov', fe
quietly: gen tsample=e(sample)

/**************************************************************************************
*********Helmert transformations and partialing out covariates************************/
quietly: gen id=newcl
sort newcl year
quietly: reg `yvar' `excov' if tsample==1
quietly: predict `yvar'_res if tsample==1, resid
quietly: helm `yvar'_res
rename h_`yvar'_res h_`yvar'
drop `yvar'_res

local num_seqexg=0
local seqexg_helm
foreach var of local seqexg{
sort newcl year
local num_seqexg=`num_seqexg'+1
quietly: gen seqexg`num_seqexg'=`var'
quietly: reg seqexg`num_seqexg' `excov' if tsample==1
quietly: predict seqexg`num_seqexg'_res if tsample==1, resid
quietly: helm seqexg`num_seqexg'_res
rename h_seqexg`num_seqexg'_res h_seqexg`num_seqexg'
local seqexg_helm `seqexg_helm' h_seqexg`num_seqexg'
drop seqexg`num_seqexg' seqexg`num_seqexg'_res
}

/***************************************************************
***************Creation of GMM instruments**********************
***************************************************************/
local gmmlist
local num_gmm=0
foreach var of local gmminst{
sort newcl year
local num_gmm=`num_gmm'+1
quietly: gen abond`num_gmm'=`var'
quietly: reg abond`num_gmm' `excov' if tsample==1
quietly: predict abond`num_gmm'_res, resid
drop abond`num_gmm' 
rename abond`num_gmm'_res abond`num_gmm'

quietly: replace abond`num_gmm'=0 if abond`num_gmm'==.
local gmmlist `gmmlist' abond`num_gmm' 

}

local gmmlist_trunc
local num_gmm_trunc=0
foreach var of local gmmtrunc{
sort newcl year
local num_gmm_trunc=`num_gmm_trunc'+1
quietly: gen abond_trunc`num_gmm_trunc'=`var'
quietly: reg abond_trunc`num_gmm_trunc' `excov' if tsample==1
quietly: predict abond_trunc`num_gmm_trunc'_res, resid
drop abond_trunc`num_gmm_trunc' 
rename abond_trunc`num_gmm_trunc'_res abond_trunc`num_gmm_trunc'
quietly: replace abond_trunc`num_gmm_trunc'=0 if abond_trunc`num_gmm_trunc'==.
local gmmlist_trunc `gmmlist_trunc' abond_trunc`num_gmm_trunc' 
}

/***************************************************************
***************Estimator  for year j between maxy and  start***/
***************************************************************/
sort newcl year

*initialize objects*
local obs=0
quietly: gen samptemp=.
matrix def Num1=J(`num_seqexg', 1, 0)
matrix def Num2=J(`num_seqexg', `num_seqexg', 0) 
matrix def Den1=J(`num_seqexg', `num_seqexg', 0) 

forvalues maxyear=`ystart'(1)`yfinal'{

local maxinst=`maxyear'-`ydeep' /*deeper gmm lags until ydeep. Warning: many of these instruments may be zero if no data*/
local maxtrunc=min(`truncate',`maxinst')
local j=`maxyear'-1960+1

*ivreg to check for collinearities and obtain degree of overidentification*
cap ivreg2 h_`yvar' (`seqexg_helm'=l(1/`maxinst').(`gmmlist') l(1/`maxtrunc').(`gmmlist_trunc'))  if year==`maxyear',   noid  noconstant

if _rc==0{
/*Runs k-class estimator*/
local lambda=1+e(sargandf)/e(N)
cap ivreg2 h_`yvar' (`seqexg_helm'=l(1/`maxinst').(`gmmlist') l(1/`maxtrunc').(`gmmlist_trunc'))  if year==`maxyear',  k(`lambda') nocollin coviv  noid noconstant

if _rc==0{
quietly: replace samptemp=e(sample) if year==`maxyear'
local mobs=e(N)
local obs=`obs'+`mobs'

/*Construct locals with restrictions*/
local restriction
local mresults
forvalues m=1(1)`num_seqexg'{
local restriction `restriction' & _b[h_seqexg`m']!=0
local mresults `mresults' (seqexg`m': _b[h_seqexg`m'])
}

/*Extracts results and weights them by the adjusted variance*/
if _rc==0   `restriction'  {
quietly: nlcom `mresults', post
matrix b=e(b)
matrix V=e(V)
cap matrix Num1=Num1+`mobs'*inv(V)*b'
cap matrix Num2=Num2+`mobs'^2*inv(V)
cap matrix Den1=Den1+`mobs'*inv(V)
}
}
}
}

/******************Compiles results and post them******************/
matrix est2top=inv(Den1)*Num1
matrix var2top=inv(Den1)*Num2*inv(Den1)
matrix b=est2top'
matrix V=var2top
mat colnames b = `seqexg' 
mat colnames V = `seqexg' 
mat rownames V = `seqexg' 
/*Countries in sample*/
quietly: bysort newcl: egen newsamp=max(samptemp)
quietly: replace newsamp=0 if newsamp!=1
/*Post estimation results*/
ereturn post b V, obs(`obs') depname("Dep var") esample(newsamp)

drop tsample samptemp id  h_`yvar' `seqexg_helm'  `gmmlist' `gmmlist_trunc' 
end


use "$path/DDCGdata_final.dta", clear
tsset wbcode2 year 
gen instrument=demreg

*Creation of spatial lags*
forvalues j=1960(1)2010{
preserve
keep if year==`j'
keep if cen_lon!=.&cen_lat!=.&y!=.&dem!=.
keep wbcode2 cen_lon cen_lat y dem year 
sort wbcode2

spmat idistance dmat cen_lon cen_lat, id(wbcode2) dfunction(dhaversine) norm(row) replace
spmat lag  y_w dmat y
spmat lag  dem_w dmat dem

keep y_w dem_w wbcode2 year
tempfile m
save `m', replace
restore

merge 1:1 wbcode2 year using `m', update
drop _m

}

sort wbcode2 year
********************************************************************************
* End of verbatim from the code provided by ANRR

sort wbcode2 year
gen newcl=wbcode2

* Use AB sample to create available number of observations
xtabond2 y l(1/4).y dem yy* ,  gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan

gen c= 1 if e(sample)
sort wbcode
by wbcode: egen csum=sum(c)
tab csum
matrix HHKTable2Column11=J(36,9,.)
* 39 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem

local i = 1
cls
tsset wbcode2 year
local z "5	6	8	9	11	12	13	14	15	16	17	18	19	20	21	22	24	25	26	27	28	29	30	31	34	35	36	37	38	39	40	41	43	44	45	46"
foreach k of local z{
	preserve
	keep if csum>=`k'
	sort wbcode2 year
	replace newcl=wbcode2
	tsset newcl year
	bootstrap _b, seed(12345) reps(100) cluster(wbcode2) idcluster(newcl):  ///
		hhkBS y (dem L1.y L2.y L3.y L4.y) (y dem) () (yy*), ydeep(1960) ystart(1964) yfinal(2009) truncate(4)
	* Min number of observations
	mat HHKTable2Column11[`i',1]=`k'
	* Dem coefficient (short-run)
	mat HHKTable2Column11[`i',6]=_b[dem]
	mat HHKTable2Column11[`i',7]=_se[dem]	
	qui nlcom (_b[dem]/(1-_b[l.y]-_b[l2.y]-_b[l3.y]-_b[l4.y])) (_b[l.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat HHKTable2Column11[`i',4]=b_2
	mat HHKTable2Column11[`i',5]=se_2
	* Dem coefficient (long-run)
	mat HHKTable2Column11[`i',8]=b_1
	mat HHKTable2Column11[`i',9]=se_1
	* Determine number of obs and countries (same as AB; bug in ANRR code)
	qui xtabond2 y l(1/4).y dem yy*   ,  gmmstyle(y, laglimits(2 .)) gmmstyle(dem, laglimits(1 .)) ivstyle(yy*   , p) noleveleq robust nodiffsargan
	* Total number of observations
	mat HHKTable2Column11[`i',2]=e(N)
	* Number of Countries
	mat HHKTable2Column11[`i',3]=e(N_clust)
	di in gr "Result #" in ye `i' in gr " for four lags: Completed results for " in ye `k' in gr " or more observations. N=" in ye e(N) in gr " obs."

	local i = `i'+1
	restore
}	

mat list HHKTable2Column11
svmat HHKTable2Column11, names(HHK211_)
rename HHK211_1 HHK211_minobs
rename HHK211_2 HHK211_obs
rename HHK211_3 HHK211_N
rename HHK211_4 HHK211_b_lagy
rename HHK211_5 HHK211_se_lagy
rename HHK211_6 HHK211_b_dem
rename HHK211_7 HHK211_se_dem
rename HHK211_8 HHK211_b_lr_dem
rename HHK211_9 HHK211_se_lr_dem

preserve
keep HHK211_*
gen model = _n if !missing(HHK211_minobs)
drop if _n>36
order model
sort HHK211_minobs
save  "$otherdata/HHK_2_11_reduction_2022.dta", replace
restore

***
*** ANRR IV regressions (with 4 lags of the dependent variable)
***

* The following is taken verbatim from the code provided by ANRR
********************************************************************************
use "$path/DDCGdata_final.dta", clear
gen instrument=demreg

*Creation of spatial lags*
forvalues j=1960(1)2010{
cd "$path"
preserve
keep if year==`j'
keep if cen_lon!=.&cen_lat!=.&y!=.&dem!=.
keep wbcode2 cen_lon cen_lat y dem year 
sort wbcode2

spmat idistance dmat cen_lon cen_lat, id(wbcode2) dfunction(dhaversine) norm(row) replace
spmat lag  y_w dmat y
spmat lag  dem_w dmat dem

keep y_w dem_w wbcode2 year
tempfile m
save `m', replace
restore

merge 1:1 wbcode2 year using `m', update
drop _m

}
********************************************************************************
* End of verbatim

preserve 
matrix IVTable6Column2=J(36,9,.)
* 36 sample reductions; columns for: min obs, total obs, countries, b_lagy, se_lagy, b_dem, se_dem, b_lr_dem, se_lr_dem, ols_b_lagy, ols_se_lagy, ols_lr_dem, ols_se_lr_dem

local i = 1
cls
tsset wbcode2 year
local z "6	7	9	10	12	13	14	15	16	17	18	19	20	21	22	23	25	26	27	28	29	30	31	32	35	36	37	38	39	40	41	42	44	45	46	47"
foreach k of local z{
	qui: xtivreg2 y l(1/4).y (dem=l(1/4).instrument) yy* if csum>=`k', fe cluster(wbcode2) r partial(yy*)
	* Min number of observations
	mat IVTable6Column2[`i',1]=`k'
	* Total number of observations
	mat IVTable6Column2[`i',2]=e(N)
	* Number of Countries
	mat IVTable6Column2[`i',3]=e(N_g)
	* Dem coefficient (short-run)
	mat IVTable6Column2[`i',6]=_b[dem]
	mat IVTable6Column2[`i',7]=_se[dem]	
	qui nlcom (_b[dem]/(1-(_b[l1.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]))) (_b[l1.y]+_b[l2.y]+_b[l3.y]+_b[l4.y]), post
	mat b1 = e(b)
	mat v1 = e(V)
	scalar b_1 = b1[1,1]
	scalar se_1 = sqrt(v1[1,1])	
	scalar b_2 = b1[1,2]
	scalar se_2 = sqrt(v1[2,2])	
	* Lag y coefficient	
	mat IVTable6Column2[`i',4]=b_2
	mat IVTable6Column2[`i',5]=se_2
	* Dem coefficient (long-run)
	mat IVTable6Column2[`i',8]=b_1
	mat IVTable6Column2[`i',9]=se_1	
	di in gr "Table 6(2), Result #" in ye `i' in gr ": Completed results for " in ye `k' in gr " or more observations. n=" in ye e(N) in gr " obs and N=" in ye e(N_g) in gr " countries." 

	local i = `i'+1
}	

mat list IVTable6Column2
svmat IVTable6Column2, names(IV62_)
rename IV62_1 IV62_minobs
rename IV62_2 IV62_obs
rename IV62_3 IV62_N
rename IV62_4 IV62_b_lagy
rename IV62_5 IV62_se_lagy
rename IV62_6 IV62_b_dem
rename IV62_7 IV62_se_dem
rename IV62_8 IV62_b_lr_dem
rename IV62_9 IV62_se_lr_dem
keep IV62_*
gen model = _n if !missing(IV62_minobs)
drop if _n>36
order model
sort IV62_minobs
save  "$otherdata/IV_6_2_reduction_2022.dta", replace
restore


***
*** Merge estimates from the four empirical models
***


cls
use "$otherdata/FE_2_3_reduction_2022.dta", clear
gen FE23_t_lr_dem = FE23_b_lr_dem/FE23_se_lr_dem

merge model using "$otherdata/AB_2_7_reduction_2022.dta", sort
gen AB27_t_lr_dem = AB27_b_lr_dem/AB27_se_lr_dem
drop _merge*

merge model using "$otherdata/HHK_2_11_reduction_2022.dta", sort
gen HHK211_t_lr_dem = HHK211_b_lr_dem/HHK211_se_lr_dem
drop _merge*

merge model using "$otherdata/IV_6_2_reduction_2022.dta", sort
gen IV62_t_lr_dem = IV62_b_lr_dem/IV62_se_lr_dem



***
*** Democracy estimates - Panel (a)
***

twoway ///
	(connected FE23_b_lr_dem FE23_minobs, lpattern(solid) ///
		lw(.8) msymbol(O) lcolor(pink*.75) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter FE23_b_lr_dem FE23_minobs if abs(FE23_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*.75) msize(medlarge) mfcolor(pink*.75) mlcolor(black) mlwidth(thin)) ///
	(connected AB27_b_lr_dem AB27_minobs, lpattern(dash) ///
		lw(.8) msymbol(O) lcolor(pink*1) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter AB27_b_lr_dem AB27_minobs if abs(AB27_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*1) msize(medlarge) mfcolor(pink*1) mlcolor(black) mlwidth(thin)) ///		
	(connected HHK211_b_lr_dem HHK211_minobs, lpattern(shortdash) ///
		lw(.8) msymbol(O) lcolor(pink*1.5) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter HHK211_b_lr_dem HHK211_minobs if abs(HHK211_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*1.5) msize(medlarge) mfcolor(pink*1.5) mlcolor(black) mlwidth(thin)) ///		
	(connected IV62_b_lr_dem IV62_minobs, ///
		lw(.8) msymbol(O) lcolor(pink*2) msize(medlarge) mfcolor(white) mlcolor(black) mlwidth(thin)) ///
	(scatter IV62_b_lr_dem IV62_minobs if abs(IV62_t_lr_dem)>1.645, ///
		msymbol(O) lcolor(pink*2) msize(medlarge) mfcolor(pink*2) mlcolor(black) mlwidth(thin)) ///	
	(scatter IV62_b_lr_dem IV62_minobs, yaxis(2) msymbol(none))  ///
, xsize(10) ysize(5) scheme(s2mono) graphregion(color(white)) plotregion(color(white)) ///
	xtitle("Minimum observation count", size(medsmall)) ///
	ytitle("Long-run Coefficient on Democracy", size(medsmall)) yscale(range(-14.5 33)) ///
	ylabel(-10(10)30, labsize(small) glcolor(gs14) glwidth(medthin) grid) ///
	yscale(range(-14.5 33) axis(2)) ///
	ytitle("Long-run Coefficient on Democracy", size(medsmall) axis(2))  ///
	ylabel(-10(10)30, labsize(small) grid glc(black) ///
		glw(.2) glp(shortdash)  axis(2)) /// 
	xlabel(5(5)45) xmtick(##5) yline(0, lpattern(dash) lw(.2))  ///
	legend(rows(2) order(- "Long-run Estimate on Democracy (4-lags):" 1 3 5 7 - "Statistically significant at 10% level:" 2 4 6 8) ///
	label(1 "2FE") ///
	label(2 "") ///
	label(3 "AB") ///
	label(4 "") ///
	label(5 "HHK") ///
	label(6 "") ///
	label(7 "IV") ///
	label(8 "") ///
	size(medsmall)) 
graph export "$texfolder/Comparison_para_evolution_full.eps", as(eps) replace
