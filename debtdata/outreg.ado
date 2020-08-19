*! version 3.1.9  2aug08  by jgallup@cal.berkeley.edu  
*! Write formatted regression output to a text file
*
* Changelog: 
* 4/22/05 - fixed bug in adjr2 option for areg command
* 8/22/08 - fixed bug for appending more than 4 estimates (with new independent variables)
*
* syntax: outreg [varlist] using filename [,
*
*                /* text options */
*                noLabel TItle(passthru) CTitle(passthru)
*                noNOTes ADDNote(passthru)  
*
*                /* coefficient options */
*                BDec(numlist int >=0 <=11) BFmt(passthru) COEfastr
*
*                /* t stat/se/pvalue/ci options */
*                SE Pvalue CI BEtaco Level(integer $S_level) 
*					  TDec(numlist int >=0 <=11 max=1) noPAren BRacket 
*					  noASter 3aster 10pct SIgsymb(passthru)
* 
*                /* statistics options */
*                noCONs noNOBs noNI noR2 ADJr2 RDec(numlist int >=0 <=11 max=1) 
*                ADDStat(passthru) ADEc(numlist int >=0 <=11)
*                EForm Margin1 Margin2(string) Onecol Xstats 
*
*                /* other options */
*                COMma APpend REPLACE Quote]
*

capture program drop outreg
program define outreg
   version 6.0
   * write formatted regression output to file

   syntax [varlist(default=none ts)] using [, noLabel TItle(passthru)        /*
    */ CTitle(passthru) noNOTes ADDNote(passthru) BDec(numlist int >=0 <=11) /*
    */ BFmt(passthru) COEfastr SE Pvalue CI BEtaco Level(integer $S_level)   /*
    */ Tdec(numlist int >=0 <=11 max=1) noPAren BRacket noASter 3aster       /*
    */ 10pct SIgsymb(passthru) noCONs noNOBs noNI noR2 ADJr2                 /*
    */ RDec(numlist int >=0 <=11 max=1) ADDStat(passthru)                    /*
    */ ADEc(numlist int >=0 <=11) EForm Margin1 Margin2(string) Onecol       /*
    */ Xstats COMma APpend REPLACE Quote] 

   if ("`se'"!="")+("`pvalue'"!="")+("`ci'"!="")+("`betaco'"!="")>1 { 
      di in red "choose only one of se, pvalue, ci, or beta"
      exit 198 
   }
   if `level'<10 | `level'>99 { 
      di in red "level() invalid"
      exit 198
   }
   if "`aster'"=="noaster" & /*
    */ ("`3aster'"!="" | "`10pct'"!="" | "`coefastr'"!="" | "`sigsymb'"!="") {
      if "`3aster'"!="" {
         di in red "cannot choose both noaster and 3aster options"}
      else if "`10pct'"!="" {
         di in red "cannot choose both noaster and 10pct options"}
      else if "`coefastr'"!="" {
         di in red "cannot choose both noaster and coefastr options"}
      else {di in red "cannot choose both noaster and sigsymb options"}
      exit 198
   }
   else if ("`3aster'"!="" & "`10pct'"!="") {
      di in red "cannot choose both 3aster and 10pct options"
      exit 198
   }
   if (`"`title'"'!="" & "`append'"!="") {
      di in red "warning: title ignored in appended regressions"
   }
   if (`"`addnote'"'!="" & "`append'"!="") {
      di in red "warning: addnote ignored in appended regressions"
   }
   if "`bdec'"=="" {local bdec = 3}

   if "`tdec'"=="" {
      if "`ci'"=="ci" {local tdec = `bdec'}
      else if "`pvalue'"=="pvalue" {local tdec = 3}
      else {local tdec = 2}
   }
   if "`rdec'"=="" {local rdec = `tdec'}
   if (`"`addstat'"'=="" & "`adec'"!="")  {
      di in red "cannot choose adec option without addstat option"
      exit 198
   }  
   if "`adec'"=="" {local adec = `tdec'}
   if "`onecol'"!="" & "`varlist'"!="" {
      di in red "cannot specify varlist with the onecol option"
      exit 198
   }
   if "`quote'"!="quote" {local quote "noquote"}
   
   tempname regN df_r rsq numi r2mat
   tempname varname coefcol b vc b_alone

   /* for svy commands with subpop(), N_sub is # of obs used for estimation */
   local cmd = e(cmd)
   local svy = substr("`cmd'",1,3)
   if "`svy'"=="svy" & e(N_sub) != . {scalar `regN' = e(N_sub)}  
   else {scalar `regN'   = e(N)}
   
   scalar `df_r'   = e(df_r)
   scalar `numi'   = e(N_g)

   local   depvar  = e(depvar)
   local   robust  = e(vcetype)
   if "`robust'"=="." {local robust "none"}
   local   ivar    = e(ivar)
   mat     `b'     = e(b)
   mat     `vc'    = e(V)
   local   univar  = 1          /* univariate (single equation) regression */
   local   bcols   = colsof(`b')  /* cols of b */
   local   bocols  = `bcols'      /* cols of b only, w/o other stats */
   capture local fracpol = (e(fp_cmd)=="fracpoly")
    

   if "`margin1'"~="" | "`margin2'"~="" {
      local margin = "margin"
      if "`margin2'"~="" {
         local margucp "margucp(_`margin2')"
         scalar `df_r' = .
         if "`margin1'"~="" {
            di in red "may not specify both margin and margin()"
            exit 198
         }
      }
      else {
         if "`cmd'"=="tobit" {
            di in red "dtobit requires margin({u|c|p}) after dtobit command"
            exit 198
         }
      }
   }

   /* parse addstat to convert possible r(), e(), and s() macros to numbers */ 
   /*   (to avoid conflicts with r-class commands used in this program) */
   if `"`addstat'"'!="" {
      gettoken part rest: addstat, parse(" (")
      gettoken rest zilch: rest, parse(" (") match(parns) /* strip off "addstat()" */
      local i = 1
      while `"`rest'"' != "" {
         gettoken name rest : rest, parse(",") quote
         if `"`name'"'=="" {
            di in red "empty strings not allowed in addstat() option"
            exit 6
         }
         gettoken acomma rest : rest, parse(",") 
         gettoken valstr rest : rest, parse(",")
         /* strip off comma between addstats */
         if `"`rest'"' != "" {gettoken comma2 rest: rest, parse(",")}
         else {local comma2 ""}
         local value = `valstr'
         capture confirm number `value'
         if _rc!=0 {
            di in red `"`valstr' found where number expected in addstat() option"'
            exit 7
         }
         local aadec : word `i' of `adec'
         if "`aadec'"=="" {local aadec `prvadec'}
         local valstr = string(`value',"%12.`aadec'f") 
         local newadd `"`newadd'`name'`acomma'`valstr'`comma2'"'
         local prvadec = `aadec'
         local i = `i'+1
      }
      local addstat `"`newadd'"'
   }

   /* logistic reports coeffients in exponentiated form (odds ratios) */
   if "`cmd'"=="logistic" {local eform "eform"}
   if "`eform'"=="eform" {local cons "nocons"}
   
   /* multi-equation models */
   if "`cmd'"=="mvreg" | "`cmd'"=="sureg" | "`cmd'"=="reg3" | "`cmd'"=="gologit2" {
      local univar  = 0      /* multivariate (multiple equation) regression */
      if "`onecol'" != "onecol" {
         local neq     = e(k_eq)
         local eqlist  "`e(eqnames)'"
         local depvar "`eqlist'"

			if "`cmd'"=="mvreg" | "`cmd'"=="sureg" | "`cmd'"=="reg3" {
				mat `r2mat'   = e(b)   /* get column labels */
				if "`cmd'"=="mvreg" {local r2list  = e(r2)}
				local eq = 1
				while `eq' <= `neq' {
					if "`cmd'"=="mvreg" { 
						local r2str: word `eq' of `r2list'
						scalar `rsq' = real("`r2str'")
					}
					else {scalar `rsq' = e(r2_`eq')}
					mat `r2mat'[1,`eq'] = `rsq'
					local eq = `eq' + 1
				}
			}
      }
      else {   /* if onecol */
         local r2 = "nor2"  
         scalar `rsq' = .
      }
   }     /* `rsq' after `r2list' to avoid type mismatch  */
   else if "`adjr2'"!="" {
      scalar `rsq' = e(r2_a)
      if `rsq' == . {
         scalar `rsq' = e(ar2)
         if `rsq' == . {
            di in red "Adjusted R-squared (e(r2_a)) not defined; cannot use adjr2 option"
            exit 198
         }
      }
   }
   else {scalar `rsq' = e(r2)}

   if ("`cmd'"=="intreg" | "`cmd'"=="svyintrg" | "`cmd'"=="xtintreg"){ 
      local depvar : word 1 of `depvar'  /* 2 depvars listed */
   }

   /* nolabels for anova and fracpoly */
   else if ("`cmd'"=="anova" | `fracpol' | "`cmd'"=="nl") {  
      local label "nolabel"
   }

   /* margin or dprobit: substitute marginal effects into b and vc */
   else if ("`cmd'"=="dprobit" | "`margin'"=="margin") {
      if "`cmd'"=="dlogit2" | "`cmd'"=="dprobit2" | "`cmd'"=="dmlogit2" {
         di in red "warning: margin option not needed"
      }
      else {
         marginal, b(`b') vc(`vc') `se' `margucp'
         local bcols  = colsof(`b')  /* cols of b */
         local bocols = `bcols'      /* cols of b only, w/o other stats */
         if "`cmd'"=="dprobit" {local cons "nocons"}
      }
   }

   /* dlogit2 and dprobit2: no e(cmd) */
   else if ("$S_E_cmd"=="dlogit2" | "$S_E_cmd"=="dprobit2") {
      local depvar $S_E_depv
   }
   
   /* multi-equation: arch, biprobit, gnbreg, heckman, heckprob, zinb, zip */
   else if ("`cmd'"=="arch" | "`cmd'"=="biprobit" | "`cmd'"=="gnbreg" | /*
      */ "`cmd'"=="heckman" | "`cmd'"=="heckprob" | "`cmd'"=="treatreg" | /*
      */  "`cmd'"=="zinb"   | "`cmd'"=="zip") {
      local univar  = 0      /* multivariate (multiple equation) regression */
      if ("`onecol'"=="") {
         local neq = 2
         if "`cmd'"=="arch" {local depvar "`depvar' ARCH"}
         else if "`cmd'"=="gnbreg" {local depvar "`depvar' lnalpha"}
         else if "`cmd'"=="heckman" {local depvar "`depvar' select"}
         else if ("`cmd'"=="zinb" | "`cmd'"=="zip") {
            local depvar "`depvar' inflate"
         }
         local eqlist "`depvar'"
      }
      else {local depvar : word 1 of `depvar'}  /* if onecol */
   }

   /* multinomial regression */
   else if ("`cmd'"=="mlogit" | "`cmd'"=="svymlog" | "$S_E_cmd"=="dmlogit2") {
      local univar = 0    /* multivariate (multinomial) regression */
      if "$S_E_cmd"=="dmlogit2" {local depvar $S_E_depv}
      if "`onecol'"!="onecol" {
         if "$S_E_cmd"=="dmlogit2" {
            qui tab `depvar'
            local neq = r(r)-1
            local cmd "dmlogit2"
         }
         else {local neq   = e(k_cat)-1}
         mat `vc' = e(V)
         local vcols = colsof(`vc')
         local k = `vcols'/`neq'
         /* workaround for failure of -local catlist : coleq(`vc')- */
         /*   when category labels have spaces in them */
         local rest : colfullnames(`vc')
         local c 1
         while `c' <= `vcols' {
            gettoken acat rest : rest, parse(:)
            capture confirm number `acat'
            if _rc==0 {
               capture confirm integer number `acat'
               if _rc!=0 {
                  di in red "outreg only accepts integer-valued dependent variables for `cmd'"
                  exit 7
               }
            }
            local catlist `"`catlist' "`acat'""'
            gettoken colon rest : rest, parse(:)
            gettoken avar rest : rest, parse(" ")
            local c = `c' + 1
         }
         local depvar ""
         local cat = 1
         while `cat' <= `neq' {
            local catwrdn = (`cat'-1)*`k' + 1 
            local catname: word `catwrdn' of `catlist'
            local depvar `"`depvar' "`catname'""'
            local cat = `cat' + 1
         }
         local eqlist `"`depvar'"'
      }
   }

   /* when e(b) includes extra statistics (which don't have variable labels) */
   capture mat `b_alone' = `b'[1,"`depvar':"]
   if _rc==0 {local bocols = colsof(`b_alone')}
   else if ("`cmd'"=="ologit" | "`cmd'"=="oprobit") {
      local bocols = e(df_m)
      mat `b_alone' = `b'[1,1..`bocols']
   }
   else if ("`cmd'"=="cnreg" | ("`cmd'"=="tobit" & "`margin'"~="margin")) {
      local bocols = `bocols'-1 /* last element of e(b) is not est coef */
      mat `b_alone' = `b'[1,1..`bocols']
   }
   else if ("`cmd'"=="intreg" | "`cmd'"=="svyintrg") {
      mat `b_alone' = `b'[1,"model:"]
      local bocols = colsof(`b_alone')
   }
   else if ("`cmd'"=="truncreg") {
      mat `b_alone' = `b'[1,"eq1:"]
      local bocols = colsof(`b_alone')
   }
   if `univar' & `bcols'>`bocols' & "`xstats'"!="xstats" {
      mat `b' = `b_alone'
      mat `vc' = `vc'[1..`bocols',1..`bocols']
   }
   if "`xstats'"=="xstats"& `bcols'==`bocols' {
      di in red "warning: no extra statistics - xstats ignored"
   }
   
   
   /* create table with coeftxt and append to existing table */   
   quietly {
      preserve
      
      if `univar' | "`onecol'"=="onecol" {  /* make univariate regression table */
			coeftxt `varlist', `se' `pvalue' `ci' `betaco' level(`level') `label' bdec(`bdec') `bfmt' tdec(`tdec') `paren' `bracket' `aster' `3aster' `10pct' `sigsymb' `coefastr' `cons' `eform' `nobs' `ni' `r2' `adjr2' rdec(`rdec') `title' `ctitle' addstat(`addstat') adec(`adec') `notes' `addnote' `apptmp' regN(`regN') df_r(scalar(`df_r')) rsq(scalar(`rsq')) numi(scalar(`numi')) ivar(`ivar') `onecol' depvar(`depvar') robust(`robust') borows(`bocols') b(`b') vc(`vc') varname(`varname') coefcol(`coefcol') univar(`univar')
        	if "`append'"=="append" {
            appfile `using', varname(`varname') coefcol(`coefcol')
            outsheet v* `coefcol' `using', nonames `quote' `comma' replace
            drop v*
         }
         else {
            noi outsheet `varname' `coefcol' `using', nonames `quote' `comma' `replace'
         }
      }
      else {   /* multiple equation regression table */
         if `"`ctitle'"'!="" { /* parse ctitle */
            partxtl `"`ctitle'"'
            local nct = r(numtxt) /* number of ctitles */
            local n = 1
            while `n'<=`nct' {
               local ctitl`n' `r(txt`n')'
               local n = `n'+1
            }
         }
         else {local nct=0}
         tempname b_eq vc_eq
         local depvlst `"`depvar'"'
         local eq = 1
         while `eq' <= `neq' {
            if `eq' <= `nct' {local ctitle `"ctitle(`ctitl`eq'')"'}
            local eqname: word `eq' of `eqlist'
            local depvar: word `eq' of `depvlst'

            /* r2mat doesn't exist for mlogit => "capture" */
            capture scalar `rsq' = `r2mat'[1,`eq'] 
            mat `b_eq'  = `b'[.,"`eqname':"]
            matrix colnames `b_eq' = _:  /* remove roweq from b_eq for explicit varlist */
            mat `vc_eq' = `vc'["`eqname':","`eqname':"]
            local bocols = colsof(`b_eq')
   
            if `eq'>1 {local addstat ""}
            if `eq' == 1 & "`append'"!="append" {local apptmp ""}
            else {local apptmp "append"}
				coeftxt `varlist', `se' `pvalue' `ci' `betaco' level(`level') `label' bdec(`bdec') `bfmt' tdec(`tdec') `paren' `bracket' `aster' `3aster' `10pct' `sigsymb' `coefastr' `cons' `eform' `nobs' `ni' `r2' `adjr2' rdec(`rdec') `title' `ctitle' addstat(`addstat') adec(`adec') `notes' `addnote' `apptmp' regN(`regN') df_r(scalar(`df_r')) rsq(scalar(`rsq')) numi(scalar(`numi')) ivar(`ivar') `onecol' depvar(`depvar') robust(`robust') borows(`bocols') b(`b_eq') vc(`vc_eq') varname(`varname') coefcol(`coefcol') univar(`univar')
            if `eq' == 1 & "`append'"!="append" {
               noi outsheet `varname' `coefcol' `using', nonames `quote' `comma' `replace'
            }
            else {
               appfile `using', varname(`varname') coefcol(`coefcol') 
               outsheet v* `coefcol' `using', nonames `quote' `comma' replace
               drop v*
            }
            local eq = `eq' + 1
            restore, preserve  /* to access var labels after first equation */
         } 
      }
      
   }  /* for quietly */
end

capture program drop coeftxt
program define coeftxt
   version 6.0
   * convert vectors b and vc into coefficient and t statistics text + variable name text

   syntax [varlist(default=none ts)] [, SE Pvalue CI BEtaco Level(integer $S_level) noLabel BDEC(numlist) BFmt(passthru) TDEC(numlist) noPAren BRacket noASter 3aster 10pct SIgsymb(passthru) COEfastr noCONs EForm noNOBs noNI noR2 ADJr2 RDec(numlist) TItle(passthru) CTitle(string) ADDStat(passthru) ADEC(numlist) noNOTes ADDNote(passthru) APpend regN(string) df_r(string) rsq(string) numi(string) ivar(string) Onecol depvar(string) robust(string) BOROWS(string) b(string) vc(string) varname(string) coefcol(string) univar(string)] 

   tempvar beta st_err tstat astrix mrgrow
   if "`betaco'"=="betaco" {tempname betcoef}
   tempfile bcopy 
   tempname b_alone vc_alon b_xtra vc_xtra tcrit1 tcrit5 tcrit10 t_alpha

   mat `b' = `b''
   mat `vc' = vecdiag(`vc')
   mat `vc' = `vc''
   local brows = rowsof(`b')

   if `"`ctitle'"' != "" {local coltitl `"`ctitle'"'}
   else {
      if "`label'"=="" {capture local coltitl : var label `depvar'}
      if `"`coltitl'"'=="" {local coltitl "`depvar'"}
   }

   if (`numi'!=. & "`ni'"!="noni") {
      if ("`label'"=="" & "`ivar'"!=".") {local iname : var label `ivar'}
      if `"`iname'"'==""  {local iname "`ivar'"}
      if `"`iname'"'=="." {local iname "groups"}
   }

   /* fill in "beta" "st. err." & "tstat" variables from regression output */
   /* in case varlist is specified:  */
   mat `b_alone' = `b'[1..`borows',1] /* use to avoid _cons in xtra stats */
   if `brows'>`borows' {  /* allow for xtra stats */
      local borows1 = `borows'+1
      mat `b_xtra' = `b'[`borows1'...,1...]
      mat `vc_xtra' = `vc'[`borows1'...,1...]
   }
   if "`varlist'"!="" {
      tempname arow testnew newb newvc
        /* add the constant unless "nocons" is chosen */
      if "`cons'"!="nocons" {local varlist "`varlist' _cons"}
      local vname : word 1 of `varlist'
      local i=1
      while "`vname'"!="" {
         local j = rownumb(`b_alone',"`vname'")
         if `j'!=. {
            matrix `arow'  = `b'[`j',1...]  /* "..." needed to get rownames */
            matrix `newb'  = nullmat(`newb')\ `arow'
            matrix `arow'  = `vc'[`j',1...]
            matrix `newvc' = nullmat(`newvc')\ `arow'
         }
         else if (`univar' & "`vname'"!="_cons") {
            di in red "`vname' not found in regression coefficients"
            exit 111
         }
         local i = `i'+1
         local vname : word `i' of `varlist'
      }
      mat `b_alone' = `newb'
      if `brows'>`borows' {  /* allow for xtra stats */
         mat `newb'  = `newb'\ `b_xtra'
         mat `newvc' = `newvc'\ `vc_xtra'
      }
      mat `b' = `newb'
      mat `vc' = `newvc'
   }
   else if "`cons'"=="nocons" {   /* delete the constant if "nocons" is chosen */
      local j_1 = rownumb(`b_alone',"_cons")-1
      if `j_1'!=. {  /* in case there is no constant in b */
         mat `b_alone' = `b_alone'[1..`j_1',1...]
         mat `vc_alon' = `vc'[1..`j_1',1...]
         if `brows'==`borows' {  
            mat `b' = `b_alone'
            mat `vc' = `vc_alon'
         }
         else {  /* allow for xtra stats */
            mat `b' = `b_alone' \ `b_xtra'
            mat `vc' = `vc_alon' \ `vc_xtra'
         }
      }
   }
   local borows = rowsof(`b_alone')
   local brows = rowsof(`b')  /* reset brows */

   gen double `beta' = matrix(`b'[_n, 1]) in 1/`brows'
   gen double `st_err' = matrix(`vc'[_n, 1]) in 1/`brows'
   replace `st_err' = sqrt(`st_err')
   gen double `tstat' = abs(`beta'/`st_err')
   if "`pvalue'"=="pvalue" {
      if `df_r'==. {replace `tstat' = 2*(1-normprob(`tstat'))}
      else{replace `tstat' = tprob(`df_r', `tstat')}
   }
   if "`eform'"=="eform" {  /* exponentiate beta and st_err */
      replace `beta' = exp(`beta')
      replace `st_err' = `beta'*`st_err'
   }
   if "`betaco'"=="betaco" { /* create "beta" coefficients */
      sum `depvar' if e(sample)
      gen `betcoef' =  `beta'/r(sd) in 1/`brows'
   }

   /* fill in variables names column */
   gen str31 `varname' = ""
   local bnames : rownames(`b_alone')
   local bname : word 1 of `bnames'
   if "`onecol'" != "" {
      local beqs : roweq(`b_alone')
      local beq : word 1 of `beqs'
   }
   local i 1
   while "`bname'"!="" {
      if "`bname'"!="_cons" {
         if "`label'"!="nolabel" {
            local dotloc = index("`bname'", ".") 
            if `dotloc'!=0 {    /* deal w/time series prefixes */
               local tspref = substr("`bname'",1,`dotloc')
               local bname = substr("`bname'",`dotloc'+1,.)
               capture local blabel : var label `bname'
               if (_rc!=0) {local blabel "`bname'"}
               local blabel = `"`tspref' `blabel'"'
            }
            else {
            	capture local blabel : var label `bname'
            	if (_rc!=0) {local blabel "`bname'"}
            }
         }
         else {local blabel ""}
         
         if "`betaco'"=="betaco" { /* create "beta" coefficients */
            sum `bname' if e(sample)
            replace `betcoef' =  r(sd)*`betcoef' if `i'==_n
         }
      }
      else { 
         local blabel "Constant" 
      }
      if `"`blabel'"'=="" { 
         local blabel "`bname'" 
      }
      if "`onecol'" != "" {local blabel `"`beq':`blabel'"'}
      replace `varname' = trim(`"`blabel'"') if `i'==_n
      local i = `i'+1
      local bname : word `i' of `bnames'
      local beq : word `i' of `beqs'
   }
   if `brows'>`borows' {  /* allow for xtra stats */
      local borows1 = `borows'+1
      mat `b_xtra' = `b'[`borows1'...,1...]
      local bnames : rownames(`b_xtra')
      local beqs : roweq(`b_xtra')
      local bname : word 1 of `bnames'
      local beq : word 1 of `beqs'
      local i 1
      while "`bname'"!="" {
         if "`bname'"!="_cons" {
            if "`beq'"=="_" {local blabel "`bname'"}
            else {local blabel "`beq':`bname'"}
         }
         else {local blabel "`beq':Constant"}
         replace `varname' = `"`blabel'"' if `i'+`borows'==_n
         local i = `i'+1
         local bname : word `i' of `bnames'
         local beq : word `i' of `beqs'
      }
   }

   /* keep data until after 'betaco' */
   keep if `beta'!=.
   /* get rid of original data since labels already accessed */
   keep `varname' `beta' `tstat' `st_err' `betcoef' 
   
   /* calculate asterisks for t stats (or standard errors) */
   if `"`title'"'!="" {  /* parse title */
      partxtl `"`title'"'
      local titlrow = r(numtxt)
      local t = 1
      while `t'<=`titlrow' {
         local titl`t' `r(txt`t')'
         local t = `t'+1
      }
   }
   else{local titlrow=0}
   
   if "`append'"=="append" {local appottl = 1}
   else {local appottl = `titlrow'}
     /* either an appended column (not the first regression) or has a title */
     /* i.e. need an extra line above the coefficients */
   gen `mrgrow' = 2*_n + 1 + `appottl'

   if "`aster'"!="noaster" {
      if `"`sigsymb'"'!="" {  /* load custom significance symbols into astr1-ast10 */
         gettoken part rest: sigsymb, parse(" (") 
         gettoken part rest: rest, parse(" (")     /* strip off "sigsymb(" */
         gettoken astr1 rest : rest, parse(",)")
         gettoken acomma rest : rest, parse(",")     /* strip off comma */
         gettoken astr5 rest : rest, parse(",)")
         gettoken acomma rest : rest, parse(",")     /* strip off comma */
         gettoken astr10 rest : rest, parse(",)")
      }
      else if "`3aster'"=="3aster" {
         local astr1 "***"
         local astr5 "**"
         local astr10 "*"
      }
      else {
         local astr1 "**"
         local astr5 "*"
         if "`10pct'"=="10pct" {local astr10 "+"}
      }
      if "`astr1'"!="" {local astrtxt "`astr1' significant at 1%"}
      if "`astr5'"!="" {local astrtxt "`astr5' significant at 5%; `astrtxt'"}
      if ("`10pct'"=="10pct" | "`3aster'"=="3aster") {
         local astrtxt "`astr10' significant at 10%; `astrtxt'"
      }

      if "`pvalue'"=="pvalue" {
         gen str3 `astrix' = "`astr10'" if (`tstat'<=0.10 & `tstat'!=.)
         replace  `astrix' = "`astr5'"  if (`tstat'<=0.05 & `tstat'!=.)
         replace  `astrix' = "`astr1'"  if (`tstat'<=0.01 & `tstat'!=.)
      }
      else {
         if `df_r'==. {
            scalar `tcrit1'  = invnorm(0.995)
            scalar `tcrit5'  = invnorm(0.975)
            scalar `tcrit10' = invnorm(0.95)
         }
         else {
            scalar `tcrit1'  = invt(`df_r',0.99)
            scalar `tcrit5'  = invt(`df_r',0.95)
            scalar `tcrit10' = invt(`df_r',0.90)
         }
         gen str3 `astrix' = "`astr10'" if (`tstat'>=`tcrit10' & `tstat'!=.)
         replace  `astrix' = "`astr5'"  if (`tstat'>=`tcrit5'  & `tstat'!=.)
         replace  `astrix' = "`astr1'"  if (`tstat'>=`tcrit1'  & `tstat'!=.)
      }
   }
   else {
      gen str2 `astrix' = ""
   }      
   if "`paren'"!="noparen" {
      if "`bracket'"=="bracket" {
         local lparen "["
         local rparen "]"
      }
      else {
         local lparen "("
         local rparen ")"
      }
   }
   else {
      local lparen ""
      local rparen ""
   }
   if "`ci'"=="ci" {
      if `df_r'==. {scalar `t_alpha'  = invnorm( 1-(1-`level' /100)/2 )}
      else { 
      scalar `t_alpha'  = invt(`df_r', `level' /100)}
   }

   if `"`bfmt'"'=="" {
      local bfmt1 "fc"
      local bfmtcnt = 1
   }
   else{  /* parse bfmt */
      local fmttxt "fcgce"
      partxtl `"`bfmt'"'
      local bfmtcnt = r(numtxt)
      local b = 1
      while `b'<=`bfmtcnt' {
         local bfmt`b' `r(txt`b')'
         if index("`fmttxt'","`bfmt`b''")==0 {
            di in red `"bfmt element "`bfmt`b''" is not a valid number format (f,fc,e,g, or gc)"'
            exit 198
         }
         local b = `b'+1
      }
   }

   local bdeccnt : word count `bdec'   
   if `bdeccnt'>1 | `bfmtcnt'>1 {  /* fill in bdec# & bfmt# */
      local b = 1
      while `b'<=_N {
         local bdec`b' : word `b' of `bdec'
         if "`bdec`b''"=="" {local bdec`b' = `prvbdec'}
         local prvbdec "`bdec`b''"
         local b = `b'+1
      }
      local b = `bfmtcnt'+1
      while `b'<=_N {
         local b_1 = `b'-1
         local bfmt`b' "`bfmt`b_1''"
         local b = `b'+1
      }
   }


   /* first put the t statistics (or se | ci| pvalue | beta) in `coefcol' */
   if "`se'"!="se" & "`ci'"!="ci" & "`betaco'"!="betaco" {
      gen str12 `coefcol' = string(`tstat',"%12.`tdec'f")
   }  
   else if `bdeccnt'==1 & `bfmtcnt'==1 {
      if "`ci'"=="ci" {
         if "`eform'"!="eform" {
            gen str27 `coefcol' =   /*
              */ string(`beta'-`st_err'*`t_alpha',"%12.`bdec'`bfmt1'") + /*
              */ " - " + string(`beta'+`st_err'*`t_alpha',"%12.`bdec'`bfmt1'")
         }
         else {
            gen str27 `coefcol' =   /*
              */ string(exp(ln(`beta')-`t_alpha'*`st_err'/ `beta'),"%12.`bdec'`bfmt1'") + /*
              */ " - " + string(exp(ln(`beta')+`t_alpha'*`st_err'/ `beta'),"%12.`bdec'`bfmt1'") 
         }
      }
      else if "`betaco'"=="betaco" {
         gen str12 `coefcol' = string(`betcoef',"%12.`bdec'`bfmt1'")
      }
      else {gen str12 `coefcol' = string(`st_err',"%12.`bdec'`bfmt1'")}
   }
   else {
      gen str12 `coefcol' = ""
      local i 1
      while `i'<=_N {
         if "`ci'"=="ci" {
            if "`eform'"!="eform" {
               gen str27 `coefcol' =   /*
                 */ string(`beta'-`t_alpha'*`st_err',"%12.`bdec`i''`bfmt`i''") + /*
                 */ " - " + string(`beta'+`t_alpha'*`st_err',"%12.`bdec`i''`bfmt`i''") /*
                 */ if `i'==_n
            }
            else {
               gen str27 `coefcol' =   /*
                 */ string(exp(ln(`beta')-`t_alpha'*`st_err'/ `beta'),"%12.`bdec`i''`bfmt`i''") + /*
                 */ " - " + string(exp(ln(`beta')+`t_alpha'*`st_err'/ `beta'),"%12.`bdec`i''`bfmt`i''") /*
                 */ if `i'==_n
            }
         }
         else if "`betaco'"=="betaco" {
            replace str12 `coefcol' = string(`betcoef',"%12.`bdec`i''`bfmt`i''") if `i'==_n
         }
         else {
            replace `coefcol' = string(`st_err',"%12.`bdec`i''`bfmt`i''") if `i'==_n
         }
         local i = `i'+1
      }
   }
   if "`coefastr'" != "coefastr" {
      replace `coefcol' = "`lparen'" + `coefcol' + "`rparen'" + `astrix'
   }
   else {replace `coefcol' = "`lparen'" + `coefcol' + "`rparen'"}
   sort `mrgrow'
   save "`bcopy'", replace         /* double quotes needed for Macintosh   */
                                   /*   problems with spaces in file names */

   /* then offset the row number and put the coefficients in coefcol */
   replace `mrgrow' = `mrgrow'-1   
   if `bdeccnt'==1 & `bfmtcnt'==1 { 
      replace `coefcol' = string(`beta',"%12.`bdec'`bfmt1'")
   }
   else {
      local i 1
      while `i'<=_N {
         replace `coefcol' = string(`beta',"%12.`bdec`i''`bfmt`i''") if `i'==_n
         local i = `i'+1
      }
   }
   if "`coefastr'" == "coefastr" {replace `coefcol' = `coefcol' + `astrix'}
   sort `mrgrow'

   /* interleave the coefficients with the t statistics */
   merge `mrgrow' using "`bcopy'"
   replace `varname' = " " if _merge==2  /* no variable names next to tstats */
   drop `beta' `st_err' `tstat' `astrix' `betcoef' _merge
   sort `mrgrow'
   save "`bcopy'", replace

   drop `varname' `coefcol'


   /* offset row numbers again to add header and other statistics */

   /* first find number of new rows for addstat() */
   if `"`addstat'"'!="" {
      partxtl `"`addstat'"'
      local naddst = int((real(r(numtxt))+1)/2)
      local n = 1
      while `n'<=`naddst' {
         local t = (`n'-1)*2+1
         local astnam`n' `r(txt`t')'
         local t = `t'+1
         local astval`n' `r(txt`t')' /* pair: stat name & value */
         local n = `n'+1
      }
   }
   else{local naddst=0}
   
   /* find number of new rows for addnote() */
   if (`"`addnote'"'!="" & "`append'"!="append") {
      partxtl `"`addnote'"'
      local naddnt = r(numtxt)
      local n = 1
      while `n'<=`naddnt' {
         local anote`n' `r(txt`n')'
         local n = `n'+1
      }
   }
   else{local naddnt=0}
 
   /* calculate total number of rows in table */
   local coefrow = 2*`brows'+1+`appottl'
   local totrows = `coefrow' + ("`nobs'"!="nonobs") + ("`ni'"!="noni"&`numi'!=.) + /*
   */ ("`r2'"!="nor2"&`rsq'!=.&`df_r'!=.) + `naddst' +              /*
   */ ("`notes'"!="nonotes"&"`append'"!="append")*(1+("`aster'"!="noaster")) /*
   */ + `naddnt'

   local minobs = max(5000,`brows')
   set obs `minobs'  /* obs must be set to >= 100 (depending on matsize?) */
   keep if _n <= `totrows'
   replace `mrgrow' = _n
   sort `mrgrow'
   merge `mrgrow' using "`bcopy'"
   sort `mrgrow'
   drop _merge
   
   if `appottl' {
      if (`"`title'"'!="" & "`append'"!="append") {
         local row = 1
         while `row'<=`titlrow' {
            replace `varname' = `"`titl`row''"' if _n==`row'
            local row = `row'+1
         }
         replace `coefcol' = `"`coltitl'"' if _n==`row'
      }
      else {replace `coefcol' = `"`coltitl'"' if _n==2}
   }
   else {
      replace `coefcol' = `"`coltitl'"' if _n==1
   }
   if "`nobs'"!="nonobs" {
      local coefrow = `coefrow'+1
      replace `varname' = "Observations" if _n==`coefrow'
      replace `coefcol' = string(`regN') if _n==`coefrow'
   }
   if (`numi'!=. & "`ni'"!="noni") {
      local coefrow = `coefrow'+1
      replace `varname' = "Number of " + rtrim(`"`iname'"') if _n==`coefrow'
      replace `coefcol' = string(`numi') if _n==`coefrow'
   }
   if "`r2'"!="nor2" & `rsq'!=. & `df_r'!=. { /* if df_r=., not true r2 */
      local coefrow = `coefrow'+1
      replace `coefcol' = string(`rsq',"%12.`rdec'f") if _n==`coefrow'
      replace `varname' = "R-squared" if _n==`coefrow'
      if "`adjr2'"=="adjr2" {
         replace `varname' = "Adjusted " + `varname' if _n==`coefrow'
      }
   }
   if `"`addstat'"'!="" {
      local i 1
      local adeccnt : word count `adec'
      while `i'<=`naddst' {
         local coefrow = `coefrow'+1
         local aadec : word `i' of `adec'
         if "`aadec'"=="" {local aadec `prvadec'}
         if `"`astval`i''"'!="" {
            replace `coefcol' = "`astval`i''" if _n==`coefrow'
         }
         replace `varname' = trim(`"`astnam`i''"') if _n==`coefrow'
         local i = `i'+1
         local prvadec `aadec'
      }
   }
   if ("`notes'"!="nonotes" & "`append'"!="append") {
      local coefrow = `coefrow'+1
      if "`bracket'"=="bracket" {local par_bra "brackets"}
      else {local par_bra "parentheses"}
      if "`pvalue'"=="pvalue" {local statxt "p values"}
      else if "`se'"=="se" {local statxt "Standard errors"}
      else if "`ci'"=="ci" {local statxt "`level'% confidence intervals"}
      else if "`betaco'"=="betaco" {local statxt "Normalized beta coefficients"}
      else {
         if `df_r'!=. {local t_or_z "t"}
         else {local t_or_z "z"}
         local statxt "`t_or_z' statistics"
         if "`robust'"=="none" {local statxt "Absolute value of `statxt'"}
      }
      if "`robust'"=="Robust" {local statxt = "Robust " + lower("`statxt'")}
      replace `varname' = "`statxt' in `par_bra'" if _n==`coefrow'
      if "`aster'"!="noaster" {
         local coefrow = `coefrow'+1
         replace `varname' = "`astrtxt'" if _n==`coefrow'
      }
   }
   if (`"`addnote'"'!="" & "`append'"!="append") {
      local i 1
      while `i'<=`naddnt' {
         local coefrow = `coefrow'+1
         replace `varname' = `"`anote`i''"' if _n==`coefrow'
         local i = `i'+1
      }
   }
end  /* coeftxt */

capture program drop appfile
program define appfile
   version 6.0
   * append regression results to pre-existing file
   
   syntax using/, varname(string) coefcol(string)

   tempname vord1 vord2 vartmp varsml v2plus vorder merge2 
   tempfile tmpf1
   gen str80 `vartmp' = substr(`varname',1,79) /* room for "!" at end */
   replace `vartmp' = "0" if _n==1
   gen str10 `varsml' = trim(`vartmp')
   replace `vartmp' = `vartmp'[_n-1]+"!" if `varsml'==""
      /* add "!" to variable name to make it sort after previous variable name */
      /* will cause bug if both "varname" and "varname!" already exist */
   count if (`varsml'=="" | (`varsml'[_n+1]=="" & _n!=_N))
   local ncoeff2 = r(N) /* number of estimated coefficients in file 2 */
   local N2 = _N        /* number of lines in file 2 */
   gen `vord2' = _n     /* ordering variable for file 2 */
   sort `vartmp'
   keep `varname' `coefcol' `vartmp' `vord2'
   save `"`tmpf1'"', replace
   
   if index(`"`using'"', ".")==0 {local using `"`using'.out"'}
   insheet using `"`using'"', nonames clear
   describe, short
   local numcol = r(k) /* number of columns already in file 1 */
   gen str80 `vartmp' = substr(v1,1,79)  /* room for "!" at end */
   local titlrow = (v1[1]!="")
   if `titlrow' {   /* count title rows */
      while v1[`titlrow'+1]!="" {
         local titlrow = `titlrow'+1
      }
   }
   local frstrow = 1 + `titlrow' /* first non-title row */
   replace `vartmp' = "0" if _n==`frstrow' & v2=="(1)"
   replace `vartmp' = "0!" if _n==`frstrow' & v2!="(1)"
   replace `vartmp' = `vartmp'[_n-1]+"!" if `vartmp'==""
   gen long `vord1' = _n
   gen str80 `v2plus' = trim(v2)
   local col = 3
   while `col'<=`numcol' {
      replace `v2plus' = `v2plus' + trim(v`col')
      local col = `col'+1
   }
   count if ((v1==""&`v2plus'!="")  /* i.e. a t stat or column heading
         */ | (v1[_n+1]=="" & (`v2plus'[_n+1]!=""|_n==1) & _n!=_N)) 
         /* i.e. a coefficient (next row is a t stat) */
   local ncoeff1 = r(N)
   gen str10 `varsml' = `vartmp'
   summ `vord1' if `vord1'>`ncoeff1' & `v2plus'!="" /* v2plus for addstat */
   local endsta1 = r(max)  /* calc last row of statistics before notes */
   if `endsta1'==. {local endsta1 = `ncoeff1'}
   drop `varsml'

   sort `vartmp'
   merge `vartmp' using "`tmpf1'"
   gen str10 `varsml' = `vartmp'
   gen `vorder' = 1 if (`vord1'<=`ncoeff1' | `vord2'<=`ncoeff2') /* coefficients */
   replace `vorder' = 0 if ((`vord1'<=`titlrow') | (`vartmp'=="0" & _merge==2))  /* has title or has no column numbers */
   replace `vorder' = 2 if (`varsml'=="Constant" | `varsml'=="Constant!") /* constant */
   replace `vorder' = 3 if `vorder'==. & (`vord1'<=`endsta1' | `vord2'<=`N2') /* statistics */
   replace `vorder' = 4 if `vorder'==.  /* notes below statistics */
   gen byte `merge2' = _merge==2
   sort `vorder' `vord1' `vord2' `merge2'

   replace v1 = `varname' if v1=="" & `varname'!=""
   drop `varname' `vartmp' `varsml' `vorder' `vord1' `vord2' `merge2' _merge
   if (`numcol'==2) {
      replace v2 = "(1)" if _n==`frstrow'
      replace `coefcol' = "(2)" if _n==`frstrow'
   }
   else {
      replace `coefcol' = "(" + string(`numcol') + ")" if _n==`frstrow'
   }
end  /* appfile */

capture program drop marginal
program define marginal
   version 6.0
   * put marginal effects (dfdx) into b and vc matrices 
   syntax , b(string) vc(string) [se margucp(string)]

   tempname dfdx se_dfdx new_vc dfdx_b2      
   capture mat `dfdx' = e(dfdx`margucp')
   if _rc==0 {
      local cnam_b : colnames `dfdx'
      local cnam_1 : word 1 of `cnam_b'
   }
   if _rc!=0 {
      if "`cnam_1'"=="c1" {
         di in red `"Update dprobit ado file: type "help update" in Stata"'
      }
      else {di in red "margin option invalid: no marginal effects matrix e(dfdx`margucp') exists"}
      exit
   }
   /* create matrix of diagonals for vc */
   if "`se'"=="se" {
      if e(cmd)=="dprobit" | e(cmd)=="tobit" {
         if e(cmd)=="dprobit" {local margucp "_dfdx"}
         mat `se_dfdx' = e(se`margucp')
         mat `vc' = diag(`se_dfdx')
         mat `vc' = `vc' * `vc'
      }
      else {mat `vc' = e(V_dfdx)}
   mat colnames `vc' = `cnam_b'
   }
   else {
   /* if t or p stats reported then trick `cv' into giving the right t stat */
      local coldfdx = colsof(`dfdx')
      mat `new_vc' = J(`coldfdx',`coldfdx',0)
      local i = 1
      while `i' <= `coldfdx' {
         scalar `dfdx_b2' = (el(`dfdx',1,`i')/el(`b',1,`i'))^2
         mat `new_vc'[`i',`i'] = `dfdx_b2'*`vc'[`i',`i']
         local i = `i'+1
      }
      mat colnames `new_vc' = `cnam_b'
      mat `vc' = `new_vc'
   }  
   mat `b' = `dfdx'
end

capture program drop partxtl
program define partxtl, rclass   
/* parse text list to find number of text elements and return them */
   local ntxt = 0
   gettoken part rest: 1, parse(" (") 
   gettoken rest zilch: rest, parse(" (") match(parns) /* strip off "option()" */
   while `"`rest'"' != "" {
      local ntxt = `ntxt'+1
      gettoken part rest: rest, parse(",")
      return local txt`ntxt' `"`part'"'
      gettoken part rest: rest, parse(",") /* strip off "," */
   }
   return local numtxt `ntxt'
end  /* end partxt1 */

