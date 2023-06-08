// Trevor Williams, 28 August 2019. Revised by Giuseppe Moscarini, June 2021

set more off
set matsize 11000
clear mata
capture log close main
capture log close results
clear

global workingdir "D:\Kevin\Firm_Growth_Project\_Final\Sales"


global moments     "$workingdir\Output\moments"
global densities   "$workingdir\Output\densities"
global figures     "$workingdir\Output\figures"
global singledir   "$figures\single\ihs"
global paperdir    "$figures\paper\ihs"
capture mkdir      "$moments"
capture mkdir      "$densities"
capture mkdir      "$figures"
capture mkdir      "$singledir"
capture mkdir      "$paperdir"
* quietly log using "$moments/compute_moments.log", replace name(main)

tempfile temp_ADM
tempfile autocorr_cell
tempfile autocorr_age_class
tempfile autocorr_size_class

use "$workingdir\Data\clean\new_vat_data_cleaned.dta", clear


// use ihs or halt changes for the growth rate

drop stzd_gro_log
drop stzd_gro_simpl
drop stzd_gro_halt
rename stzd_gro_ihs stzd_gro


drop if entry == 1

**************************************************************
*  COMPUTATION OF GROWTH RATE MOMENTS BY SIZE AND AGE CLASS  *
**************************************************************

tsset firm_id age

* lag and forward of growth rate
* NOTE: original standardized growth rate stzd_gro is forward, between t and t+1
gen stzd_f_gro = f.stzd_gro /* growth rate between t+1 and t+2 */
gen stzd_l_gro = l.stzd_gro /* growth rate between t-1 and t   */

* tag sets indicator firm=1 in the firm's first year in the sample
egen firm=tag(firm_id)

save `temp_ADM', replace

* compute changes in growth rate in adjacent periods
gen double f_dgro=stzd_f_gro-stzd_gro  /* growth rate from t+1 to t+2 minus growth rate from t to t+1 */
gen double dgro=stzd_gro-stzd_l_gro   /* growth ratre from t to t+1 minus growth rate from t-1 to t */

* total just summarizes variables 
total firm size_sh
local tot_firm=_b[firm] // total number of firms
local tot_size=_b[size_sh] // total size share (across all years)

****************************
*  AUTOCORRELATION BY AGE  *
****************************
**** current persistence
preserve
	bys age: egen sd_stzd_l_gro = sd(stzd_l_gro) /* stdev of lagged growth rate by age (NOTE: Not equal to 1, because standardization is by calendar quarter, not by age)*/
	bys age: egen sd_stzd_gro = sd(stzd_gro)     /* stdev of current (forward) growth rate by age */
	gen stzd_l_gro_unitsd = stzd_l_gro/sd_stzd_l_gro /* normalize lagged growth rate by its stdev */
	gen stzd_gro_unitsd = stzd_gro/sd_stzd_gro       /* normalize current growth rate by its stdev */
 	statsby rho_c=_b[stzd_l_gro_unitsd] se = _se[stzd_l_gro_unitsd], by(age) clear: regress stzd_gro_unitsd stzd_l_gro_unitsd
	label var rho_c "Current persistence"
	label var se "Standard error"
	save "$moments/current_persistence.dta", replace
restore



******* Computation of first 4 moments by: *************************************
		* age class
		* size class
		* age and size class (cell)

save `temp_ADM', replace
		
foreach class of varlist age_class size_class cell {
	use `temp_ADM', clear
	preserve
	
        * counts of firms 
		keep cell size_sh stzd_gro f_dgro age_class size_class exit entry firm firm_id
		bys `class': gegen double pjk=total(firm)
		replace pjk=pjk/`tot_firm' // cell's share of all firms, by number
		bys `class': gegen double sharejk=total(size_sh)
		replace sharejk=sharejk/`tot_size' // cell's share of all sales
		
		bys `class': gegen double tot_size=total(size_sh)
		gen double share=size_sh/tot_size // firm's share of cell's sales
		bys `class': gen wexit = exit*share
		bys `class': gegen double xjk=sum(wexit) // exit rate in cell
		
		* mean growth, weighted by sale shares, all
		* weighs each individual growth rate by sale share, then adds up
		gen double wgro=share*stzd_gro
		bys `class': gegen double gjk=total(wgro)
		
        * mean growth, weighted by sale shares, survivors
		* first, recompute normalized sales share of survivors only, within the class
		
		bys `class': gegen double tot_size_surv=total(size_sh) if exit==0
		gen double share_surv=size_sh/tot_size_surv if exit==0
	   
		******************************

		* If exit==1, double tot_size_surv will be missing,
		* hence so will the renormalized slaes share of survivors, share_surv
		* therefore, the growth rate multiplied by sale share will also be missing
		gen double wgro_surv=share_surv*stzd_gro
		* and not contribute to the total, which is the mean growth rate of survivors
		bys `class': gegen double gjk_surv=total(wgro_surv)
		
		******************************

		* variance of growth, weighted by sale shares, all and survivors
		gen double vijk=((stzd_gro-gjk)^2)*share
		bys `class': gegen double vjk=total(vijk)
		bys `class': gegen double vdjk=sd(f_dgro) [iweight=size_sh] // stdev of forward growth rates
		replace vdjk=vdjk^2
		gen double vijk_surv=((stzd_gro-gjk_surv)^2)*share_surv
		bys `class': gegen double vjk_surv=total(vijk_surv) // * does total include exiting firms?
		
		******************************

		* skewness, weighted by sale shares, all and survivors
		gen skwijk=((stzd_gro-gjk)^3)*share
		bys `class': gegen double skwjk=total(skwijk)
		replace skwjk=skwjk/(vjk^1.5)
		gen skwijk_surv=((stzd_gro-gjk_surv)^3)*share_surv
		bys `class': gegen double skwjk_surv=total(skwijk_surv)
		replace skwjk_surv=skwjk_surv/(vjk_surv^1.5)
		
		* kurtosis, weighted by sale shares, all and survivors
		gen kurijk=((stzd_gro-gjk)^4)*share
		bys `class': gegen double kurjk=total(kurijk)
		replace kurjk=kurjk/(vjk^2)
		gen kurijk_surv=((stzd_gro-gjk_surv)^4)*share_surv
		bys `class': gegen double kurjk_surv=total(kurijk_surv)
		replace kurjk_surv=kurjk_surv/(vjk_surv^2)
		
******* Variance of moments ****************************************************
		
		bys `class': gen double var_xjk=xjk*(1-xjk)/(pjk*`tot_firm')
		gen double share_sq=share^2
		gen double share_surv_sq=share_surv^2
		bys `class': egen double sum_share_sq=total(share_sq)
		bys `class': egen double sum_share_surv_sq=total(share_sq)
		bys `class': gen double var_gjk=sum_share_sq*vjk
		bys `class': gen double var_gjk_surv=sum_share_surv_sq*vjk_surv
		
		bys `class': gen nfirms = _N
		
    // compute standard errors of exit rate, mean growth, and variance of growth
		
		/* exit rate (stderr of a Bernoulli) */
		bys `class': gen se_xjk = xjk*(1-xjk) / sqrt(nfirms)
		
		/* mean growth (stdev divided by sqrt of sample size) */
		bys `class': gen se_gjk = sqrt(vjk) / sqrt(nfirms)
		
		/* variance of growth  (formula for chi^2( */
		bys `class': gen se_vjk = vjk * sqrt(kurjk - (nfirms-3)/(nfirms-1)) / sqrt(nfirms)
		
		keep cell size_class age_class pjk sharejk xjk gj* vj* vdj* skwj* kurj* var_* nfirms se*
		order cell size_class age_class pjk sharejk xjk gj* vj* vdj* skwj* kurj* var_*, first
		bys `class': keep if _n==1
	// 	replace gjk=gjk + $mean_gro
	// 	replace gjk_surv=gjk_surv + $mean_gro_surv
		
		save "$moments/moments_`class'.dta", replace
	restore
	
	
******* Autocorrelations of growth rates* *************************************	

	keep cell size_class age_class stzd_gro stzd_f_gro stzd_l_gro f_dgro dgro
	gen double rho_c=.
	gen double rho_f=.
	gen double rho_dg=.
	
	levelsof(`class'), local(levels)
	foreach jk of local levels {
		preserve
				keep if `class' == `jk'
				capture {
					quietly regress stzd_gro stzd_l_gro
					estadd summ, sd
					estadd ysumm, sd
					matrix sd=e(sd)
					quietly replace rho_c=_b[stzd_l_gro]*sd[1,1]/e(ysd)
				}
				capture {
					quietly regress stzd_f_gro stzd_gro
					estadd summ, sd
					estadd ysumm, sd
					matrix sd=e(sd)
					quietly replace rho_f=_b[stzd_gro]*sd[1,1]/e(ysd)
				}
				capture {
					quietly regress f_dgro dgro
					estadd summ, sd
					estadd ysumm, sd
					matrix sd=e(sd)
					quietly replace rho_dg=_b[dgro]*sd[1,1]/e(ysd)
				}
				keep if _n==1
				keep size_class age_class rho_* cell
				capture append using `autocorr_`class''
				save `autocorr_`class'', replace
		restore
	}
	use `autocorr_`class'', clear
	merge 1:1 `class' using "$moments/moments_`class'.dta", nogen
	order pjk sharejk xjk gj* vj* vdj* skwj* kurj*, first
	outsheet using "$moments/moments_`class'.txt", comma nolabel replace
	save "$moments/moments_`class'.dta", replace // version(13)
}



*******************
* PREPARE VARIABLES FOR PLOTTING
*******************

rename gjk mean_growth
rename gjk_surv mean_growth_surv
rename xjk exit_rate
rename pjk firm_count_distr
rename sharejk size_distr
rename vjk variance_growth
rename vjk_surv variance_growth_surv
rename vdjk variance_dgrowth
rename skwjk skew_growth
rename skwjk_surv skew_growth_surv
rename kurjk kurt_growth
rename kurjk_surv kurt_growth_surv
rename rho_c  persist_c_growth_f
rename rho_f  persist_f_growth_f
rename rho_dg persist_dgrowth_f
destring *, replace



order size_class age_class firm* size* mean* variance* skew* kurt* exit* persist*
	  	  		  
* label empirical variables to use their proper names in graph titles		  
lab var size_class "size class"
lab var age_class "age class"
lab var mean_growth "mean growth rate"
lab var mean_growth_surv "mean growth rate, continuing"
lab var exit_rate "exit rate,"
lab var firm_count_distr "firm count distribution,"
lab var size_distr "size distribution,"
lab var variance_growth "variance of growth rates"
lab var variance_growth_surv "variance of g-rates, continuing"
lab var variance_dgrowth "variance of changes in growth rates"
lab var skew_growth "skewness of growth rates"
lab var skew_growth_surv "skewness of g-rates, continuing"
lab var kurt_growth "kurtosis of growth rates"
lab var kurt_growth_surv "kurtosis of g-rates, continuing"
lab var persist_c_growth_f "current growth persistence"
lab var persist_f_growth_f "forward growth persistence"
lab var persist_dgrowth_f "forward dgrowth persistence"
	  
* choose moments to plot and collect their names in a macro `empirical_moments'
local empirical_moments "mean_growth mean_growth_surv exit_rate firm_count_distr size_distr variance_growth variance_growth_surv variance_dgrowth skew_growth skew_growth_surv kurt_growth kurt_growth_surv persist_c_growth_f persist_f_growth_f persist_dgrowth_f"
keep `empirical_moments' *class



* save variable labels in local macros before reshape, to use them in graph titles

foreach var of local empirical_moments {
	local lab`var': variable label `var'
} 

* PLOT MOMENTS BY FIRM SIZE CLASS FOR EACH AGE CLASS
****************************************************

preserve
drop if size_class==.
drop if age_class==.

reshape wide `empirical_moments', i(size_class) j(age_class)


* merge with unconditional moments

sort size_class

merge m:1 size_class using "${moments}/moments_size_class.dta", nogen

**

*rename
rename gjk mean_growth
rename gjk_surv mean_growth_surv
rename xjk exit_rate
rename pjk firm_count_distr
rename sharejk size_distr
rename vjk variance_growth
rename vjk_surv variance_growth_surv
rename vdjk variance_dgrowth
rename skwjk skew_growth
rename skwjk_surv skew_growth_surv
rename kurjk kurt_growth
rename kurjk_surv kurt_growth_surv
rename rho_c  persist_c_growth_f
rename rho_f  persist_f_growth_f
rename rho_dg  persist_dgrowth_f

destring *, replace


* plot moments one by one, by size, for all age classes in one graph per moment
foreach var of local empirical_moments {
	global var `var'
    local all_ages

* collect moment for all ages in a local macro `all_ages'
	foreach i of numlist 1/9 {
	local all_ages `all_ages' `var'`i'
		}
	
* create graph of current variable `var' by firm size, for all ages in same graph
	line `all_ages' size_class, graphregion(fcolor(white)) ///
	name(`var'_by_size, replace) title("`lab`var'' by size", size(medsmall)) ///
	xtitle("size decile") leg(off)
    graph export "${singledir}/_${var}_by_size_${name}.pdf", replace	
	
* add unconditional moments by size (average across ages) and include it in plot
    local all_ages `var' `all_ages'
	line `all_ages' size_class, graphregion(fcolor(white)) lpattern(dash) lwidth(vthick) ///
	name(`var'_by_size_all, replace) title("`lab`var'' by size", size(medsmall)) ///
	xtitle("size decile") leg(off)	
    graph export "${singledir}/_${var}_by_size_all_${name}.pdf", replace	
}

restore

* PLOT MOMENTS BY AGE CLASS FOR EACH SIZE CLASS
***********************************************

preserve
drop if size_class==.
drop if age_class==.

reshape wide "`empirical_moments'", i(age_class) j(size_class)

sort age_class

merge m:1 age_class using "${moments}/moments_age_class.dta"
drop _merge

*rename
rename gjk mean_growth
rename gjk_surv mean_growth_surv
rename xjk exit_rate
rename pjk firm_count_distr
rename sharejk size_distr
rename vjk variance_growth
rename vjk_surv variance_growth_surv
rename vdjk variance_dgrowth
rename skwjk skew_growth
rename skwjk_surv skew_growth_surv
rename kurjk kurt_growth
rename kurjk_surv kurt_growth_surv
rename rho_c  persist_c_growth_f
rename rho_f  persist_f_growth_f
rename rho_dg  persist_dgrowth_f

destring *, replace



foreach var of local empirical_moments {
	global var `var'
 
* collect moment for all sizes in a local macro	`all_sizes'
    local all_sizes 
	foreach i of numlist 1/10 {
	local all_sizes `all_sizes' `var'`i'
		}
				
* create graph of current variable `var' by firm age, for all sizes in same graph
	line `all_sizes' age_class, graphregion(fcolor(white)) ///
	name(`var'_by_age, replace) title("`lab`var'' by age", size(medsmall)) ///
	xtitle("age class") leg(off)
	graph export "${singledir}/_${var}_by_age_${name}.pdf", replace
	
* add unconditional moments by age (average across sizes) and include it in plot
    local all_sizes `var' `all_sizes'
	line `all_sizes' age_class, graphregion(fcolor(white)) lpattern(dash) lwidth(vthick) ///
	name(`var'_by_age_all, replace) title("`lab`var'' by age", size(medsmall)) ///
	xtitle("age class") leg(off)
	graph export "${singledir}/_${var}_by_age_all_${name}.pdf", replace

}

restore
	
* CREATE AND SAVE COMBINED GRAPHS BY SIZE AND AGE FOR THE PAPER
***************************************************************	

* exit rate by size and age	  

	graph combine exit_rate_by_size exit_rate_by_age, ///
	graphregion(fcolor(white)) name(exit_rate, replace)  ycommon ///
	cols(2) iscale(1.5) ysize(3) graphregion(margin(zero)) 
	graph export "${paperdir}/_exit_rate_${name}.pdf", replace	

	graph combine exit_rate_by_size_all exit_rate_by_age_all, ///
	graphregion(fcolor(white)) name(exit_rate_all, replace)  ycommon ///
	cols(2) iscale(1.5) ysize(3) graphregion(margin(zero)) 
	graph export "${paperdir}/_exit_rate_all_${name}.pdf", replace	

	
*   mean, variance, skew, kurtosis of growth rates by size and age


* EDIT: remove ycommon after name(`var'..., replace)
foreach var in mean variance skew kurt {	
	global var `var' 
    graph combine `var'_growth_surv_by_size `var'_growth_surv_by_age  ///
	`var'_growth_by_size `var'_growth_by_age , ///
	cols(2) iscale(.7273) ysize(6) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(`var', replace)
	
	graph export "${paperdir}/_${var}_growth_${name}.pdf", replace

    graph combine `var'_growth_surv_by_size_all `var'_growth_surv_by_age_all ///
	`var'_growth_by_size_all `var'_growth_by_age_all, ///
	cols(2) iscale(.7273) ysize(6) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(`var'_all, replace)
	
	graph export "${paperdir}/_${var}_growth_all_${name}.pdf", replace
	
}

    graph combine variance_dgrowth_by_size variance_dgrowth_by_age, ///
	cols(2) iscale(.7273) ysize(6) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(variance_dgrowth, replace) ycommon
	
	graph export "${paperdir}/_variance_dgrowth_${name}.pdf", replace

    graph combine variance_dgrowth_by_size_all variance_dgrowth_by_age_all, ///
	cols(2) iscale(.7273) ysize(6) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(variance_dgrowth_all, replace) ycommon
	
	graph export "${paperdir}/_variance_dgrowth_all_${name}.pdf", replace


* current persistence of growth rates, by size and age
	
	graph combine persist_c_growth_f_by_size persist_c_growth_f_by_age, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_growth_persist, replace)  ycommon
	graph export "${paperdir}/_current_persist_growth_${name}.pdf", replace

	graph combine persist_c_growth_f_by_size_all persist_c_growth_f_by_age_all, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_growth_persist_all, replace)  ycommon
	graph export "${paperdir}/_current_persist_growth_all_${name}.pdf", replace		

* forward persistence of growth rates, by size and age
	
	graph combine persist_f_growth_f_by_size persist_f_growth_f_by_age, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_growth_persist, replace)  ycommon
	graph export "${paperdir}/_forward_persist_growth_${name}.pdf", replace

	graph combine persist_f_growth_f_by_size_all persist_f_growth_f_by_age_all, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_growth_persist_all, replace)  ycommon
	graph export "${paperdir}/_forward_persist_growth_all_${name}.pdf", replace		

* forward persistence of growth rates, by size and age
	
	graph combine persist_dgrowth_f_by_size persist_dgrowth_f_by_age, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_dgrowth_persist, replace)  ycommon
	graph export "${paperdir}/_forward_persist_dgrowth_${name}.pdf", replace

	graph combine persist_dgrowth_f_by_size_all persist_dgrowth_f_by_age_all, ///
	cols(2) iscale(.7273) ysize(3) graphregion(margin(zero)) ///
	graphregion(fcolor(white)) name(f_dgrowth_persist_all, replace)  ycommon
	graph export "${paperdir}/_forward_persist_dgrowth_all_${name}.pdf", replace		
		
	graph close

	
