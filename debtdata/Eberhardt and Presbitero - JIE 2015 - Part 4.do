******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			                    Balance
*                         Appendix Table TA-8				
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
* You require the do-file balancehelp118.do to do the subsampling.
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


***
*** Balance
***

gen c=0
sort nwbcode year
by nwbcode: replace c=1 if !missing(dly,dldebts,dlk)
* Note that differences have to be used here since a few countries otherwise 
* have one obs at the very start of the their time series which leads to 
* problems when we construct the balance analysis.
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
// This matrix is crucial since it contains the number of observations per country
mkmat csum, matrix(robs)

mata:
numeric matrix excludemissing(numeric matrix A)
{
    return(select(A, rowmissing(A):==0))
}
end

mata: betas = J(118, 1, 0)

******************************************************************************************************************************
******************************************************************************************************************************
******************************************************************************************************************************

** EACH TIME RUN/DO THE ENTIRE DO-FILE AND CHANGE THE SPECIFICATION 
** IN balancehelp.do AND BELOW

mkmat ly, matrix(z0t1)
mkmat lk, matrix(z0t2)
mkmat ldebts, matrix(z0t3)
mkmat ldebts2, matrix(z0t4)

mkmat lyT, matrix(zT0t1)
mkmat lkT, matrix(zT0t2)
mkmat ldebtsT, matrix(zT0t3)
mkmat ldebts2T, matrix(zT0t4)

mkmat llyT, matrix(zT1t1)
mkmat llkT, matrix(zT1t2)
mkmat lldebtsT, matrix(zT1t3)
mkmat lldebts2T, matrix(zT1t4)

mkmat l2lyT, matrix(zT2t1)
mkmat l2lkT, matrix(zT2t2)
mkmat l2ldebtsT, matrix(zT2t3)
mkmat l2ldebts2T, matrix(zT2t4)


// *******************************************************************************	
// Change the sum of RHS variables to reflect the specification:
// Pick one of the following.
// *******************************************************************************	
// mat zedt = z0t2+z0t3
// mat zedt = z0t2+z0t3+z0t4
*******************************************************************************
// Including cross-section averages: 
// mat zedt = z0t2+z0t3				+zT0t1+zT0t2+zT0t3		
// mat zedt = z0t2+z0t3+z0t4		+zT0t1+zT0t2+zT0t3+zT0t4
*******************************************************************************
// Including cross-section averages and one lag of cross-section averages: 
// mat zedt = z0t2+z0t3				+zT0t1+zT0t2+zT0t3			+zT1t1+zT1t2+zT1t3
// mat zedt = z0t2+z0t3+z0t4		+zT0t1+zT0t2+zT0t3+zT0t4	+zT1t1+zT1t2+zT1t3+zT1t4
// *******************************************************************************	
// Including cross-section averages and two lags of cross-section averages: 
// mat zedt = z0t2+z0t3				+zT0t1+zT0t2+zT0t3			+zT1t1+zT1t2+zT1t3		+zT2t1+zT2t2+zT2t3	
mat zedt = z0t2+z0t3+z0t4		+zT0t1+zT0t2+zT0t3+zT0t4	+zT1t1+zT1t2+zT1t3+zT1t4+zT2t1+zT2t2+zT2t3+zT2t4
// *******************************************************************************	


// Data matrix (LHS, summed RHS)
mat z0t   = z0t1 , zedt


***
*** Country-specific estimates for beta and delta: CONSTANT ONLY CASE
***

mata:
z0tm  = st_matrix("z0t")

// Loop over countries starts here
for (i=1; i<=118; i++) {

// determining the start and end values for each country i (includes missing values)
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

	z0ti = z0tm[rstart::rend,1::2]

// partial demeaning to get rid of the constant
	L = rows(z0ti)
	z1t = J(L, 1, 0)
	z2t = J(L, 1, 0)
for (l=1; l<=L; l++) {
		z1t[l] = z0ti[l,1] - (sum(z0ti[1::l,1])/l)
		z2t[l] = z0ti[l,2] - (sum(z0ti[1::l,2])/l)
	}
	z1t=z1t[2::L]
	z2t=z2t[2::L]

mata:
// estimation
	T = rows(z1t)
	ones = J(T, 1, 1)
	xdtn = log(runningsum(ones))
	ydtn = log(runningsum(z1t):^2)
	ydtn = ydtn:-ydtn[1]
	yd2tn= log(runningsum(z2t):^2)
	yd2tn= yd2tn:-yd2tn[1]
// getting rid of missing data
	datamat = z1t , ydtn , yd2tn , xdtn 
	datamat = excludemissing(datamat)
	ydtn1 = datamat[.,2]
	ydtn2 = datamat[.,3]
	xdtn  = datamat[.,4]
	bmcon1 = invsym(xdtn'*xdtn)*(xdtn'*ydtn1)
	bmcon2 = invsym(xdtn'*xdtn)*(xdtn'*ydtn2)
	deltayz= (bmcon1-bmcon2)
	betas[i]= deltayz
}

st_matrix("betas",betas)
end

svmat betas, name(diffbetas)

* delta1 is difference in betas (model with constant only)


mat drop betas
mata: betas = J(118, 1, 0)


***
*** Country-specific estimates for beta and delta: LINEAR TREND CASE
***

// Data matrix (LHS, summed RHS)
mat z0t   = z0t1 , zedt

mata:
// determining the start and end values for each country i (includes missing values)
z0tm = st_matrix("z0t")

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

	z0ti = z0tm[rstart::rend,1::2]

// double partial demeaning to get rid of the constant and linear trend
	L = rows(z0ti)
	z1t  = J(L, 1, 0)
	z1bt = J(L, 1, 0)
	z2t  = J(L, 1, 0)
	z2bt = J(L, 1, 0)	
	
	for (l=1; l<=L; l++) {
		z1bt[l] = z0ti[l,1] - (sum(z0ti[1::l,1])/l)
		z1t[l]  = z0tm[l,1] - (sum(z0tm[1::l,1])/l) - ((2/l)*sum(z1bt[1::l]))
		z2bt[l] = z0ti[l,2] - (sum(z0ti[1::l,2])/l)
		z2t[l]  = z0tm[l,2] - (sum(z0tm[1::l,2])/l) - ((2/l)*sum(z2bt[1::l]))
	}
	
	z1t=z1t[2::L]
	z2t=z2t[2::L]

// estimation
mata:
	T = rows(z1t)
	ones = J(T, 1, 1)
	xdtn = log(runningsum(ones))
	ydtn = log(runningsum(z1t):^2)
	ydtn = ydtn:-ydtn[1]
	yd2tn= log(runningsum(z2t):^2)
	yd2tn= yd2tn:-yd2tn[1]
// getting rid of missing data
	datamat = z1t , ydtn , yd2tn , xdtn 
	datamat = excludemissing(datamat)
	ydtn1 = datamat[.,2]
	ydtn2 = datamat[.,3]
	xdtn  = datamat[.,4]
	bmcon1 = invsym(xdtn'*xdtn)*(xdtn'*ydtn1)
	bmcon2 = invsym(xdtn'*xdtn)*(xdtn'*ydtn2)
	deltayz = (bmcon1-bmcon2)
	betas[i]=deltayz
}
st_matrix("betas",betas)
end

svmat betas, name(diffbetast)

* deltat1 is difference in betas (model with constant and trend)


***
*** Create matrix for countries included
***
sort nwbcode year
mat summ_cM=0
mat summ_ctM=0
mat summ_cMD=0
mat summ_ctMD=0
set seed 1102012


***
*** Start the subsampling
***
drop c csum list
set more off
forvalues l=1/107{
	qui run "$out\balancehelp118.do"
 display in gr _skip(1) `l'-1+1 _skip(1) _c
}
mat drop z0t robs


***
*** Results
***

* Subsampling: differences in beta --- summ_cM summ_cMD summ_ctM summ_ctMD
mat summ_cM = summ_cM[2..108,1]
mat summ_ctM = summ_ctM[2..108,1]
mat summ_cMD = summ_cMD[2..108,1]
mat summ_ctMD = summ_ctMD[2..108,1]
svmat summ_cM 
svmat summ_ctM
svmat summ_cMD
svmat summ_ctMD


* Mean and Median beta-differences in the N=118 sample
egen diffbetas1m=mean(diffbetas1)
egen diffbetas1md=median(diffbetas1)
egen diffbetast1m=mean(diffbetast1)
egen diffbetast1md=median(diffbetast1)


*** constant only case
* Mean and Median Estimates
sum diffbetas1, de
scalar db1m=r(mean)
scalar db1md=r(p50)
gen bmcodif=log(12)*(abs((summ_cM1)-(diffbetas1m)))
sort bmcodif
cumul bmcodif, gen(cum_cM1)
gen difconfic=abs(cum_cM1-0.95)
sort difconfic
scalar db1ci=bmcodif[1]
drop difconfic bmcodif
* db1m and db1ci
gen bmcodif=log(12)*(abs((summ_cMD1)-(diffbetas1md)))
sort bmcodif
cumul bmcodif, gen(cum_cMD1)
gen difconfic=abs(cum_cMD1-0.95)
sort difconfic
scalar db1dci=bmcodif[1]
drop difconfic bmcodif
* db1md and db1dci

*** constant and trend case
* Mean and Median Estimates
sum diffbetast1, de
scalar dbt1m=r(mean)
scalar dbt1md=r(p50)
gen bmcodif=log(12)*(abs((summ_ctM1)-(diffbetast1m)))
sort bmcodif
cumul bmcodif, gen(cum_ctM1)
gen difconfic=abs(cum_ctM1-0.95)
sort difconfic
scalar dbt1ci=bmcodif[1]
drop difconfic bmcodif
* dbt1m and dbt1c1
gen bmcodif=log(12)*(abs((summ_ctMD1)-(diffbetast1m)))
sort bmcodif
cumul bmcodif, gen(cum_ctMD1)
gen difconfic=abs(cum_ctMD1-0.95)
sort difconfic
scalar dbt1dci=bmcodif[1]
drop difconfic bmcodif
* dbt1md and dbt1dci



** CI for mean (constant)
scalar lowM	=db1m-((1/log(118))*db1ci)
scalar upM	=db1m+((1/log(118))*db1ci)
** CI for mean (constant+trend)
scalar lowMt=dbt1m-((1/log(118))*dbt1ci)
scalar upMt	=dbt1m+((1/log(118))*dbt1ci)

** CI for median (constant)
scalar lowMD	=db1md-((1/log(118))*db1dci)
scalar upMD		=db1md+((1/log(118))*db1dci)
** CI for median (constant+trend)
scalar lowMDt	=dbt1md-((1/log(118))*dbt1dci)
scalar upMDt	=dbt1md+((1/log(118))*dbt1dci)


display in gr _newline _newline _newline _skip(1) "Mean Balance (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Mean" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowM)/2  _col(56)  db1m/2    _col(71)  (upM)/2 _continue ///
	in gr _newline _newline _newline _skip(1) "Mean Balance (constant & trend)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Mean" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMt)/2  _col(56)  dbt1m/2    _col(71)  (upMt)/2	_continue ///
	in gr _newline _newline _newline _skip(1) "Median Balance (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Median" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMD)/2  _col(56)  db1md/2    _col(71)  (upMD)/2 _continue ///
	in gr _newline _newline _newline _skip(1) "Median Balance (constant & trend)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Median" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMDt)/2  _col(56)  dbt1md/2   _col(71)  (upMDt)/2	_continue ///
	in gr _newline _newline _newline

drop diffbetas1-cum_ctMD1 
scalar drop _all

