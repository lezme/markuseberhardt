preserve
	sort year nwbcode
	gen rando=runiform() if (year==1996) 
	centile rando, c(10)
	replace rando=. if rando>r(c_1)
	sum rando
	sort nwbcode year
	by nwbcode: egen boot=mean(rando)
	drop if missing(boot) 
	by nwbcode: gen cc=1 if _n==1
	gen ccsum=sum(cc)
	sum ccsum
	local cmax=r(max)
	drop cc ccsum boot rando 
	mata: betas = J(`cmax', 1, 0)
	mata: deltas = J(`cmax', 1, 0)
	sort nwbcode year

	
	** this guarantees the first observation for each country is never missing
	gen c=0
	sort nwbcode year
	by nwbcode: replace c=1 if !missing(dly,dldebts,dlk)
	by nwbcode: gen csum=sum(c)
	drop if csum==0
	drop c
	gen c=1 if csum==1
	gen list=sum(c)
	drop c csum
	order list
	gen c=1
	by nwbcode: egen csum=sum(c)
	by nwbcode: replace csum=. if _n!=1
	mkmat csum, matrix(robs)
	drop c csum

***
*** Pick the variable to be tested for summability
***
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

mata:
z0tm  = st_matrix("z0t")

// Loop over countries starts here
for (i=1; i<=`cmax'; i++) {

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

svmat betas, name(diffbeta)

qui sum diffbeta1, de
mat db1M=r(mean)
mat db1MD=r(p50)
mat summ_cM = summ_cM \ db1M
mat summ_cMD = summ_cMD \ db1MD 

mat drop betas

***
*** Country-specific estimates for beta and delta: LINEAR TREND CASE
***
mata: betas = J(`cmax', 1, 0)


// Data matrix (LHS, summed RHS)
mat z0t   = z0t1 , zedt

mata:
// determining the start and end values for each country i (includes missing values)
z0tm = st_matrix("z0t")

for (i=1; i<=`cmax'; i++) {
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

svmat betas, name(diffbetat)
qui sum diffbetat1, de
mat dbt1M=r(mean)
mat dbt1MD=r(p50)
mat summ_ctM = summ_ctM \ dbt1M
mat summ_ctMD = summ_ctMD \ dbt1MD 

mat drop dbt1MD db1M db1MD dbt1M betas
* End of bootstrap
restore
