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
	by nwbcode: replace c=1 if !missing(ly,ldebts,lk)
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
mkmat testvar, matrix(z0t)

***
*** Country-specific estimates for beta and delta: CONSTANT ONLY CASE
***
mata:
// determining the start and end values for each country i (includes missing values)
z0tm	= st_matrix("z0t")
for (i=1; i<=10; i++) {
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
	T = rows(zt)
	b = ceil(sqrt(T))+ 1
	ones = J(T, 1, 1)
	xdtn = log(runningsum(ones))
	ydtn = log(runningsum(zt):^2)
	ydtn = ydtn:-ydtn[1]
// getting rid of missing data
	datamat = zt , ydtn , xdtn 
	datamat = excludemissing(datamat)
	ydtn = datamat[.,2]
	xdtn = datamat[.,3]
	bmcon = invsym(xdtn'*xdtn)*(xdtn'*ydtn)
	betas[i] = bmcon
	deltas[i] = (bmcon-1)/2
}
st_matrix("betas",betas)
st_matrix("deltas",deltas)
end

svmat betas, name(b_)
svmat deltas, name(d_)

qui sum b_1, de
mat b1M=r(mean)
mat b1MD=r(p50)
mat summ_cM = summ_cM \ b1M
mat summ_cMD = summ_cMD \ b1MD 


***
*** Country-specific estimates for beta and delta: LINEAR TREND CASE
***

mata:
// determining the start and end values for each country i (includes missing values)
z0tm	= st_matrix("z0t")
for (i=1; i<=10; i++) {
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

svmat betas, name(bt_)
svmat deltas, name(dt_)

qui sum bt_1, de
mat bt1M=r(mean)
mat bt1MD=r(p50)
mat summ_ctM = summ_ctM \ bt1M
mat summ_ctMD = summ_ctMD \ bt1MD 

mat drop bt1MD b1M b1MD bt1M deltas betas
* End of subsampling procedure
restore
