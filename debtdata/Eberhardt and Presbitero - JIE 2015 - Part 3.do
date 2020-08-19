******************************************************************************
*													
*                    Eberhardt and Presbitero (2015)
*          Public debt and growth: Heterogeneity and non-linearity		
*													
******************************************************************************
*			                  Summability
*                         Appendix Table TA-7				
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
* You require the do-file summhelp118.do to do the subsampling.
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

* Squares and cubes
gen ldebts2=ldebts^2
gen ldebts3=ldebts^3

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
*** Summability
***
gen c=0
sort nwbcode year
by nwbcode: replace c=1 if !missing(ly,ldebts,lk)
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
drop c csum list

mata:
numeric matrix excludemissing(numeric matrix A)
{
    return(select(A, rowmissing(A):==0))
}
end


******************************************************************************************************************************
******************************************************************************************************************************
******************************************************************************************************************************

** EACH TIME RUN/DO THE ENTIRE DO-FILE!!!

***
*** Pick the variable for analysis in the full sample and the subsample
***

* ly lk ldebts ldebts2 ldebts3  dly dldebts dlk

mkmat dlk , matrix(z0t)
gen testvar=dlk 


***
*** Full-sample --- Country-specific estimates for beta and delta: CONSTANT ONLY CASE
***

mata: betas = J(118, 1, 0)
mata: deltas = J(118, 1, 0)


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
*** Full-sample --- Country-specific estimates for beta and delta: LINEAR TREND CASE
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

// double partial demeaning to get rid of the constant and linear trend
	L = rows(z0ti)
	zt = J(L, 1, 0)
	z1t = J(L, 1, 0)
	for (l=1; l<=L; l++) {
		z1t[l]  = z0ti[l] - (sum(z0ti[1::l])/l)
		zt[l]  = z0tm[l] - (sum(z0tm[1::l])/l) - ((2/l)*sum(z1t[1::l]))
	}
	zt=zt[3..L]
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

svmat betas, name(bt)
svmat deltas, name(dt)

* End of full sample analysis 
* From this we obtain b1 and d1 (constant), bt1 and dt1 (constant and trend) --- vectors of length N


***
*** Subsample analysis: Create matrix for countries included
***
sort nwbcode year
mat summ_cM=0
mat summ_ctM=0
mat summ_cMD=0
mat summ_ctMD=0
set seed 1102012


***
*** Start the subsampling --- here: b=12
***
set more off
forvalues l=1/107{
	qui run "$out\summhelp118.do"
 display in gr _skip(1) `l'-1+1 _skip(1) _c
}
drop testvar
mat drop z0t robs

* End of subsample analysis
* From this we obtain summ_cM1 summ_ctM1 (means) summ_cMD1 summ_ctMD1 (medians) -- betas for 95 iterations each


***
*** Results
***

mat summ_cM = summ_cM[2..108,1]
mat summ_ctM = summ_ctM[2..108,1]
mat summ_cMD = summ_cMD[2..108,1]
mat summ_ctMD = summ_ctMD[2..108,1]
svmat summ_cM 
svmat summ_ctM
svmat summ_cMD
svmat summ_ctMD


* Mean and Median beta-differences in the N=105 sample
egen betas1m=mean(b1)
egen betas1md=median(b1)
egen betast1m=mean(bt1)
egen betast1md=median(bt1)


*** constant only case
* Mean and Median Estimates
sum b1, de
scalar b1m=r(mean)
scalar b1md=r(p50)
sum d1, de
scalar d1m=r(mean)
scalar d1md=r(p50)
* adjust to log(12) for b=12 (same below)
gen bmco=log(12)*(abs((summ_cM1)-(betas1m)))
sort bmco
cumul bmco, gen(cum_cM1)
gen difconfic=abs(cum_cM1-0.95)
sort difconfic
scalar b1ci=bmco[1]
drop difconfic bmco
* b1m and b1ci
gen bmco=log(12)*(abs((summ_cMD1)-(betas1md)))
sort bmco
cumul bmco, gen(cum_cMD1)
gen difconfic=abs(cum_cMD1-0.95)
sort difconfic
scalar b1dci=bmco[1]
drop difconfic bmco
* b1md and b1dci

*** constant and trend case
* Mean and Median Estimates
sum bt1, de
scalar bt1m=r(mean)
scalar bt1md=r(p50)
sum dt1, de
scalar dt1m=r(mean)
scalar dt1md=r(p50)
gen bmco=log(12)*(abs((summ_ctM1)-(betast1m)))
sort  bmco
cumul  bmco, gen(cum_ctM1)
gen difconfic=abs(cum_ctM1-0.95)
sort difconfic
scalar bt1ci= bmco[1]
drop difconfic bmco
* Yields bt1m and bt1ci
gen bmco=log(12)*(abs((summ_ctMD1)-(betast1md)))
sort bmco
cumul bmco, gen(cum_ctMD1)
gen difconfic=abs(cum_ctMD1-0.95)
sort difconfic
scalar bt1dci=bmco[1]
drop difconfic bmco
* Yields bt1md and bt1dci

** CI for mean (constant)
scalar lowM	=b1m-(1/log(118))*b1ci
scalar upM	=b1m+(1/log(118))*b1ci
** CI for mean (constant+trend)
scalar lowMt=bt1m-(1/log(118))*bt1ci
scalar upMt	=bt1m+(1/log(118))*bt1ci

** CI for median (constant)
scalar lowMD	=b1md-(1/log(118))*b1dci
scalar upMD		=b1md+(1/log(118))*b1dci
** CI for median (constant+trend)
scalar lowMDt	=bt1md-(1/log(118))*bt1dci
scalar upMDt	=bt1md+(1/log(118))*bt1dci

display in gr _newline _newline _newline _skip(1) "Mean Summability (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Mean" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowM-1)/2  _col(56)  d1m    _col(71)  (upM-1)/2 _continue ///
	in gr _newline _newline _newline _skip(1) "Mean Summability (constant & trend)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Mean" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMt-1)/2  _col(56)  dt1m    _col(71)  (upMt-1)/2	_continue ///
	in gr _newline _newline _newline _skip(1) "Median Summability (constant only)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Median" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMD-1)/2  _col(56)  d1md    _col(71)  (upMD-1)/2 _continue ///
	in gr _newline _newline _newline _skip(1) "Median Summability (constant & trend)" _continue ///
	in gr  _col(41)  "IC-lower" _col(56)  "Median" _col(71)  "IC-upper"   _newline _continue ///
	in ye  _col(41) (lowMDt-1)/2  _col(56)  dt1md    _col(71)  (upMDt-1)/2	_continue ///
	in gr _newline _newline _newline
	
qui{
	drop summ_cM1 summ_ctM1 summ_cMD1 summ_ctMD1  cum_cM1 cum_cMD1 cum_ctM1 cum_ctMD1
	drop  b1 d1 bt1 dt1 
	scalar drop _all
}





