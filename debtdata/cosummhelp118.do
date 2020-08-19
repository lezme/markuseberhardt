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
	* Make adjustment between deMG1,deMG2,deCMG1,deCMG2 on the one hand 
	* and deCMG1c,deCMG2c or deCMG1cc,deCMG2cc on the other (latter have shorter panel)
	by nwbcode: replace c=1 if !missing(deCMG2cc)

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
//mkmat eMG1, matrix(z0t)
// mkmat eMG2, matrix(z0t)
// mkmat eCMG1, matrix(z0t)
// mkmat eCMG2, matrix(z0t)
// mkmat eCMG1c, matrix(z0t)
// mkmat eCMG2c, matrix(z0t)
// mkmat eCMG1cc, matrix(z0t)
mkmat eCMG2cc, matrix(z0t)


mata:
z0tm	= st_matrix("z0t")

// Loop over countries st_matrix("z0t")
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


svmat betas, name(b_)
svmat deltas, name(d_)
qui sum b_1, de
mat b1M=r(mean)
mat b1MD=r(p50)
mat summ_cM = summ_cM \ b1M
mat summ_cMD = summ_cMD \ b1MD 

mat drop b1M b1MD betas deltas
drop d_1 b_1
restore
