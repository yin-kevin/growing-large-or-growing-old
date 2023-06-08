// compute the moments

// Francesco Beraldi and Giuseppe Moscarini, Spring 2021


set more off
set matsize 11000
clear mata
capture log close main
capture log close results
clear

global workingdir "D:\Kevin\Firm_Growth_Project\_Final\Employment"


global moments     "$workingdir\Output\Moments"
global densities   "$workingdir\Output\Distributions_and_Densities"
global figures     "$workingdir\Output\Figures"
global singledir   "$figures\Single\halt"
global paperdir    "$figures\Paper\halt"
* quietly log using "$moments/compute_moments.log", replace name(main)

tempfile temp_ADM
tempfile autocorr_cell
tempfile autocorr_age_class
tempfile autocorr_size_class

timer clear

use "$workingdir\Data\clean\empl_new_vat_data_cleaned.dta", clear

// use ihs changes for the growth rate

drop stzd_gro_log
drop stzd_gro_ihs

*******************construct dummy for presence in balanced database ***********

* here "balanced" means, for comparison with PSS, firms that survive to age 20
* we use also left-censored records, firms that enter the dataset where they are already n<20 years old
* and use their 20-n observations available in the dataset, at ages n, n+1, ....20


* define min and max age we are looking at (increase by 1 because it must be used as matrix row index)
local min_t=1
local max_t=78
* define min and max age/time horizon over which we look at autocorrelation
local min_a=1

* CHANGE: 20 -> 14 (20 as opposed to 21 was to replicate PSS which we can't do with employment data because max is 14)
local max_a=14

* choose max age that the firm must reach to be in the balanced panel
local maxage=`max_a'

sort firm_id
bys firm_id: egen maxagefirm=max(age)
bys firm_id: gen is_balanced = (maxagefirm>=`maxage') 
drop maxagefirm


/* method #2: Mata

gen is_balanced=1 

mata
st_view(Z=.,.,("firm_id","age","is_balanced")) 
i_begin=1
	for(i=1; i<=(rows(Z)-1); i++) { 
if (Z[i,1]!=Z[i+1,1]) { /*if firm_id in row i is different from firm_id in row i+1 (so the age of the firm here is as old as it gets) */
if (Z[i,2]<`maxage') { /* if the maximum age reached by the firm is less than 'maxage' */
Z[|i_begin,3\i,3|]=J(i-i_begin+1,1,0)  /* then for all the rows corresponding to that firm, substitute a 0 in the is_balanced column */
}
i_begin=i+1
}
}
if (Z[rows(Z),2]<`maxage'){
Z[|i_begin,3\rows(Z),3|]=J(rows(Z)-i_begin+1,1,0)
}
"1"  /* printing this prevents "unexpected end of line" error */
end

*/


****************************************************************************
******************* REPLICATING PSS (LOG SIZE) *****************************
****************************************************************************

gen log_empl=log(empl)

** compute stardard deviation of log empl by age (unbalanced panel)
preserve
	bys age: egen  sd_log_empl = sd(log_empl) /* assigns stdev of log empl by age to all observations with the same age, duplicated */
    keep age sd_log_empl
	bysort age: keep if _n==1 /* keeps only one observation of stdev per age (the first); the others are identical, duplicates */
*	label var se "Standard error"
	save "$moments/sd_log_empl_unbal.dta", replace
restore


** compute stardard deviation of log empl by age (balanced panel)
preserve
keep if is_balanced==1
	bys age: egen  sd_log_empl = sd(log_empl)
    keep age sd_log_empl
	bysort age: keep if _n==1
*	label var se "Standard error"
	save "$moments/sd_log_empl_bal.dta", replace
restore


** compute autocorrelation of log empl at ages t, a>t (unbalanced panel)
mata
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end

* function for assigning to correct position
mata
void impute_val(t,a,auto_corr){
st_view(Z=.,.,"rho")
auto_corr[t,a+1]=Z[1,1]
}
end

* drop age=0 duplicate observations, built to have fempl as empl before entry
drop if empl==0

tsset firm_id age


forval a = 1/14 {
gen  F`a'_log_empl=f`a'.log_empl
}




forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr log_empl F`a'_log_empl if age==`t'-1
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho

}
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocorr_log_empl_unbal.dta", replace
restore


** compute autocorrelation of log empl at ages t, a>t (balanced panel)

mata
* CHANGE: 22 -> 15
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end

tsset firm_id age

forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr log_empl F`a'_log_empl if (age==`t'-1 & is_balanced==1)
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocorr_log_empl_bal.dta", replace
restore

** autocorrelations: save with PSS shape

preserve
drop _all
use "${moments}/autocorr_log_empl_unbal.dta"

* CHANGE: 23 -> 16, 22 -> 15, 21 -> 14
mata
Z=st_data(Z=.,.,.)
X=J(15,14,.)

	for(i=1; i<=15; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,16-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_log_empl_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocorr_log_empl_bal.dta"

* CHANGE: 23 -> 16, 22 -> 15, 21 -> 14
mata
Z=st_data(Z=.,.,.)
X=J(15,14,.)

	for(i=1; i<=14; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,16-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_log_empl_bal_PSSshape.dta", replace
restore



/*
* AUTOCOVARIANCES

mata
auto_corr=J(100,22,.)
end


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr log_empl F`a'_log_empl if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)

drop rho
}
}


forval t = `min_t'/`max_t' {


corr log_empl log_empl if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)

drop rho
}



preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_log_empl_unbal.dta", replace
restore

mata
auto_corr=J(100,22,.)
end




forval t = `min_t'/`max_t' {
forval a = `min_a'/`max_a' {

corr log_empl F`a'_log_empl if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


forval t = `min_t'/`max_t' {

corr log_empl log_empl if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)


drop rho
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_log_empl_bal.dta", replace
restore



forval a = `min_a'/`max_a' {
drop F`a'_log_empl
}



** autocovariances: save with PSS shape
preserve
drop _all
use "${moments}/autocov_log_empl_unbal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_log_empl_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocov_log_empl_bal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_log_empl_bal_PSSshape.dta", replace
restore


** REPEATING WITH LEVELS OF empl


** compute autocorrelation of empl at ages t, a>t (unbalanced panel)
mata
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end


* function for assigning to correct position
mata
void impute_val(t,a,auto_corr){
st_view(Z=.,.,"rho")
auto_corr[t,a+1]=Z[1,1]
}
end

tsset firm_id age

forval a = `min_a'/`max_a' {
gen  F`a'_empl=f`a'.empl
}


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr empl F`a'_empl if age==`t'-1
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho

}
}


preserve
drop _all
getmata a*=auto_corr /* a* assigns a1 a2 a3 etc to columns 1 2 3; t is on the rows */
	save "$moments/autocorr_empl_unbal.dta", replace
restore


** compute autocorrelation at ages t, a>t (balanced panel)

mata
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr empl F`a'_empl if (age==`t'-1 & is_balanced==1)
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocorr_empl_bal.dta", replace
restore


* AUTOCOVARIANCES
mata
auto_corr=J(100,22,.)
end


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr empl F`a'_empl if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)

drop rho
}
}


forval t = `min_t'/`max_t' {


corr empl empl if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)

drop rho
}





preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_empl_unbal.dta", replace
restore

mata
auto_corr=J(100,22,.)
end




forval t = `min_t'/`max_t' {
forval a = `min_a'/`max_a' {

corr empl F`a'_empl if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


forval t = `min_t'/`max_t' {

corr empl empl if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)


drop rho
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_empl_bal.dta", replace
restore



forval a = `min_a'/`max_a' {
drop F`a'_empl
}

** autocorrelations: save with PSS shape 
* transpose to bring a on rows, and truncate because PSS stop at t+a=20
preserve
drop _all
use "${moments}/autocorr_empl_unbal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_empl_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocorr_empl_bal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_empl_bal_PSSshape.dta", replace
restore



** autocovariances: save with PSS shape
preserve
drop _all
use "${moments}/autocov_empl_unbal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_empl_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocov_empl_bal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_empl_bal_PSSshape.dta", replace
restore



*/




***********************************
** PLOTTING PSS FIGURES LOG empl **
***********************************


* ST DEV OF LOG empl BY AGE
***************************
preserve
drop _all
use "${moments}/sd_log_empl_bal.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/sd_log_empl_unbal.dta"

getmata b*=X_balanced

drop b1
rename b2 sd_log_empl_b

save "${moments}/sd_log_empl_PSSshape.dta", replace

keep if age<41

line sd_log_empl_b sd_log_empl age, graphregion(fcolor(white)) title("standard deviation of log empl by firm age") xtitle("firm age") lp(solid dash) lwidth(medium thick) note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) 
label var sd_log_empl  "Unbalanced"
label var sd_log_empl_b  "Balanced"
label var age  "Age"
graph export "$singledir/sd_log_empl_PSSshape.pdf", replace

restore




* AUTOCORRELATION OF LOG empl BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocorr_log_empl_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocorr_log_empl_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocorr_log_empl_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocorr of log empl") 	xtitle("Age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off)    lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)

* CHANGE: 21 -> 14
keep if x_axis<14

gen age=x_axis-1

* CHANGE: a20 -> a14, b20 -> b14
* CHANGE: deleted 6 of the 'dash' arguments
line a1-a14 b1-b14 age, graphregion(fcolor(white)) title("autocorrelation of log empl between firm age a and a+h") xtitle("firm age a+h") note("solid=balanced, dashed=unbalanced; each line corresponds to an initial age a", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocorr_log_empl_PSSshape.pdf", replace

restore


/*

* AUTOCOVARIANCE OF LOG empl BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocov_log_empl_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocov_log_empl_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocov_log_empl_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocov of log empl") xtitle("Age") note("solid=unbalanced, dashed=balanced", pos(6)) legend(off)  lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocov of log empl") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocov_log_empl_PSSshape.pdf", replace

restore

* AUTOCORRELATION OF empl BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocorr_empl_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocorr_empl_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocorr_empl_empl_PSSshape.dta", replace


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocorr of empl") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocorr_empl_empl_PSSshape.pdf", replace

restore



* AUTOCOVARIANCE OF empl BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocov_empl_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocov_empl_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocov_empl_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocov of empl") xtitle("Age") note("solid=balanced, dashed=unbalanced", pos(6)) legend(off) lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocov of empl") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocov_empl_PSSshape.pdf", replace

restore


*/







****************************************************************************
***** REPLICATING PSS BUT THIS TIME WITH GROWTH RATES **********************
****************************************************************************

* rename growth rate variable for ease of use

rename stzd_gro_halt growth

** compute stardard deviation of empl growth by age (unbalanced panel)
preserve
	bys age: egen sd_growth = sd(growth) /* assigns stdev of log empl by age to all observations with the same age, duplicated */
    keep age sd_growth
	bysort age: keep if _n==1 /* keeps only one observation of stdev per age (the first); the others are identical, duplicates */
*	label var se "Standard error"
	save "$moments/sd_growth_unbal.dta", replace
restore


** compute stardard deviation of log empl by age (balanced panel)
preserve
keep if is_balanced==1
	bys age: egen  sd_growth = sd(growth) /* assigns stdev of log empl by age to all observations with the same age, duplicated */
    keep age sd_growth
	bysort age: keep if _n==1 /* keeps only one observation of stdev per age (the first); the others are identical, duplicates */
*	label var se "Standard error"
	save "$moments/sd_growth_bal.dta", replace
restore

** compute autocorrelation of growth rates at ages t, a>t (unbalanced panel)
mata
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end


* function for assigning to correct position
mata
void impute_val(t,a,auto_corr){
st_view(Z=.,.,"rho")
auto_corr[t,a+1]=Z[1,1]
}
end

tsset firm_id age

forval a = `min_a'/`max_a' {
gen  F`a'_growth=f`a'.growth
}


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr growth F`a'_growth if age==`t'-1
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho

}
}


preserve
drop _all
getmata a*=auto_corr /* a* assigns a1 a2 a3 etc to columns 1 2 3; t is on the rows */
	save "$moments/autocorr_growth_unbal.dta", replace
restore


** compute autocorrelation of growth rates at ages t, a>t (balanced panel)

mata
auto_corr=J(100,22,.)
auto_corr[|1,1\.,1|]=J(100,1,1)
end

tsset firm_id age

forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr growth F`a'_growth if (age==`t'-1 & is_balanced==1)
gen rho = r(rho)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocorr_growth_bal.dta", replace
restore


** autocorrelations: save with PSS shape 
* transpose to bring a on rows, and truncate because PSS stop at t+a=20
preserve
drop _all
use "${moments}/autocorr_growth_unbal.dta"

* CHANGE: 23 -> 16, 22 -> 15, 21 -> 14
mata
Z=st_data(Z=.,.,.)
X=J(15,14,.)

	for(i=1; i<=14; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,16-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_growth_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocorr_growth_bal.dta"

* CHANGE: 23 -> 16, 22 -> 15, 21 -> 14
mata
Z=st_data(Z=.,.,.)
X=J(15,14,.)

	for(i=1; i<=14; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,16-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocorr_growth_bal_PSSshape.dta", replace
restore




/*

* AUTOCOVARIANCES
mata
auto_corr=J(100,22,.)
end


forval t = `min_t'/`max_t' {

forval a = `min_a'/`max_a' {

corr growth F`a'_growth if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)

drop rho
}
}


forval t = `min_t'/`max_t' {


corr growth growth if age==`t'-1, covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)

drop rho
}





preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_growth_unbal.dta", replace
restore

mata
auto_corr=J(100,22,.)
end




forval t = `min_t'/`max_t' {
forval a = `min_a'/`max_a' {

corr growth F`a'_growth if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',`a',auto_corr)


drop rho
}
}


forval t = `min_t'/`max_t' {

corr growth growth if (age==`t'-1 & is_balanced==1), covariance
gen rho = r(cov_12)   

mata impute_val(`t',0,auto_corr)


drop rho
}


preserve
drop _all
getmata a*=auto_corr
	save "$moments/autocov_growth_bal.dta", replace
restore



forval a = `min_a'/`max_a' {
drop F`a'_growth
}



** autocovariances: save with PSS shape
preserve
drop _all
use "${moments}/autocov_growth_unbal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_growth_unbal_PSSshape.dta", replace
restore


preserve
drop _all
use "${moments}/autocov_growth_bal.dta"

mata
Z=st_data(Z=.,.,.)
X=J(22,21,.)

	for(i=1; i<=21; i++) {
	
X[|i,i\.,i|]=Z[|i,1\i,23-i|]'
}
end

drop _all

getmata a*=X

save "$moments/autocov_growth_bal_PSSshape.dta", replace
restore


*/

*******************************************
** PLOTTING PSS FIGURES FOR GROWTH RATES **
*******************************************


* ST DEV OF GROWTH BY AGE
***************************
preserve
drop _all
use "${moments}/sd_growth_bal.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/sd_growth_unbal.dta"

getmata b*=X_balanced

drop b1
rename b2 sd_growth_b

save "${moments}/sd_growth_PSSshape.dta", replace

keep if age<41

line sd_growth_b sd_growth age, graphregion(fcolor(white)) title("standard deviation of growth rates by firm age") xtitle("firm age") lp(solid dash) lwidth(medium thick) note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) 
label var sd_growth  "Unbalanced"
label var sd_growth_b  "Balanced"
label var age  "Age"
graph export "$singledir/sd_growth_PSSshape.pdf", replace

restore




* AUTOCORRELATION OF GROWTH RATES BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocorr_growth_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocorr_growth_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocorr_growth_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocorr of growth rates") 	xtitle("Age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off)    lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)

* CHANGE: 21 -> 14
keep if x_axis<14

gen age=x_axis-1

* CHANGE: a20 -> a14, b20 -> b14
* CHANGE: deleted 6 of the 'dash' arguments
line a1-a14 b1-b14 age, graphregion(fcolor(white)) title("autocorrelation of growth rates between firm age a and a+h") xtitle("firm age a+h") note("solid=balanced, dashed=unbalanced; each line corresponds to an initial age a", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocorr_growth_PSSshape_all.pdf", replace

* excluding first point, =1

* CHANGE: 
forval i=1/14 {
replace a`i'=. if a`i'==1
replace b`i'=. if b`i'==1
}

* CHANGE: a20 -> a14, b20 -> b14
* CHANGE: deleted 6 of the 'dash' arguments
line a1-a14 b1-b14 age, graphregion(fcolor(white)) title("autocorrelation of growth rates between firm age a and a+h") xtitle("firm age a+h") note("solid=balanced, dashed=unbalanced; each line corresponds to an initial age a", pos(6)) legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocorr_growth_PSSshape.pdf", replace


restore

/*

* AUTOCOVARIANCE OF GROWTH RATES BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocov_growth_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocov_growth_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocov_growth_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocov of growth rates") xtitle("Age") note("solid=unbalanced, dashed=balanced", pos(6)) legend(off)  lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocov of growth rates") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocov_growth_PSSshape.pdf", replace

restore

* AUTOCORRELATION OF GROWTH RATES BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocorr_growth_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocorr_growth_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocorr_growth_PSSshape.dta", replace


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocorr of growth rates") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocorr_growth_PSSshape.pdf", replace

restore



* AUTOCOVARIANCE OF GROWTH RATES BETWEEN AGES a AND a+h, FOR EACH a by h
*********************************************************************

preserve
drop _all
use "${moments}/autocov_growth_bal_PSSshape.dta"

mata
X_balanced=st_data(.,.)
end

drop _all

use "${moments}/autocov_growth_unbal_PSSshape.dta"

getmata b*=X_balanced

gen x_axis=0

mata
st_view(Z=.,.,"x_axis")
	for(i=1; i<=rows(Z); i++) {
Z[i]=i	
	}
end

save "${moments}/autocov_growth_PSSshape.dta", replace

* line `longlist' x_axis, graphregion(fcolor(white)) title("Autocov of empl") xtitle("Age") note("solid=balanced, dashed=unbalanced", pos(6)) legend(off) lp(solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid solid dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash ) lwidth(medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium medium thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick thick)


keep if x_axis<21

gen age=x_axis-1

line a1-a20 b1-b20 age, graphregion(fcolor(white)) title("autocov of growth rates") xtitle("firm age") note("solid=balanced, dashed=unbalanced", pos(6))  legend(off) lp(dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash dash)

graph export "$singledir/autocov_growth_PSSshape.pdf", replace

restore

*/

