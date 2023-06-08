/*
 
Trevor Williams, 7 September 2019. Revised by Giuseppe Moscarini, June 2021
Recent edits by Kevin Yin, February 2022

This do-file computes a panel of firms with age, standardized growth rate, and a standardized size decile. 
Annual values are quarterly averages. It produces a cleaned dataset with the following variables:

1. firm_id. 
The original dataset contains "vat" which is a VAT tax number. 
Occasionally, a vat reports no sales for several consecutive quarters and/or years,
then re-enters. In that case, we assume it is a new firm and assign a new firm_id to the old vat.

2. age. Ê
The original dataset contains the variable age_trunc, which is truncated at 25.
This variable, when <25, increases even if vat is mothballed (not reporting sales), i.e. it is the age of the VAT number.
The original dataset also contains birth_y, which is the birth_year of the firm,
and start, startq, which are the year and quarter in which the VAT id enters the dataset 
(first report of positive sales). 
Sometimes, a firm born in, say, 1960, enters the dataset in 2007, at age 47, 
even if the dataset begins in 1996.

We create age of the firm, defined as current year_quarter minus "birth year joined with start quarter." 

EXAMPLE: a firm born in 1960 that enters the dataset in 2008q2 has, 
in 2010q1, age (2010q1-1960q2)/4 = 199/4 =  49 (rounded).

ENTRY EXAMPLE: a firm born in 2003q2 that enters the dataset then has,
through 2004q1, age 0. Age in the first incomplete year of operations is 0.
If the firm is born in, e.g., 2005q1, then age 0 corresponds to 2005q1:2005q4, 
age 1 corresponds to 2006q1:2006q4, etc. if the firm is born in, e.g., 
2005q3, then age 0 corresponds to 2005q3:2006q2.

3. size_class: decile of the firm's standardized size in that year 
(decile by quarter, averaged across quarters and rounded)

4. stzd_gro_ihs: annual forward growth rate (at t: between t and t+1) 
using the inverse hyperbolic sine transformation: 
g_t = IHS(s_t+1) - IHS(s_t) for each quarterly s
and then standardized by the mean and sdev of growth of all firms within the quarter.
Then the annual rate is the average quarterly rate

5. stdz_gro_size_log: annual forward growth rate using the log transformation: 
g_t = log(s_t+1) - log(s_t) for each quarterly s; then the annual rate is the average quarterly rate
and then standardized by the mean and sdev of growth of all firms within the quarter.
Then the annual rate is the average quarterly rate

6. exit: equal to 1 in the last year of positive sales


notes: 
- We used the registered age of the firm. 
This will overstate the age of the firm in the case of (potentially unobserved) re-entry

*/

/******************* LOAD, ORGANIZE, AND RESHAPE THE DATA *********************/
clear
capture log close
set more off
global workingdir "D:\Kevin\Firm_Growth_Project\_Final\Sales\Data"

use "$workingdir\raw\new_vat_data_07_21_ADM.dta"


******** create dataset of forward lagged deflators and total size
bys year: egen total_size = total(size)
duplicates drop year deflator, force
drop if year == 1995

* the data has some erroneous deflators, use the deflator most common to each year
bysort year: egen flag=min(deflator)
duplicates tag year, gen(dup)
drop if dup>0 & deflator != flag
rename deflator cpi

* generate lags
tsset year
gen f_cpi = cpi[_n+1]
gen f_total_size = total_size[_n+1]

* only one incorrect observation for 1995, no future size for final year
drop if year == 1995

* save dataset of lags
keep year cpi f_cpi total_size f_total_size
save "$workingdir\raw\lagged_vars.dta", replace



******** eliminate entries with age=0 for pre-entry info

clear
use "$workingdir\raw\new_vat_data_07_21_ADM.dta"

sort firm_id age year	
drop if firm_id==firm_id[_n+1]&age==0&age[_n+1]==0&size==0
tsset firm_id age


/**************** DETECT ENTRY (AFTER 1996q1) AND EXIT (BEFORE 2019) ********/
sort firm_id year
gen entry = (firm_id~=firm_id[_n-1]) & (age==0) // dummy equal to 1 in first sample quarter the firm operates
gen exit = (firm_id~=firm_id[_n+1]) & (year<2018) // dummy equal to 1 in last year firm is in the data, if before 2018 - this will always be in the fourth calendar quarter
replace exit=1 if year==2018&f_size==0 /* firms that exit in 2018-2019 */


/***** recode missing to zeros for the computation of the growth rates? *******
recode size f_size (. = 0)
recode size f_size (-1 = 0)

****** NO! when size is missing, so should the growth rate
we leave f_size=0 when the firm is exiting

*/

recode size f_size l_size (-1 = .)
recode size l_size (0 = .)
replace f_size=. if f_size==0&exit==0


/**************** COMPUTE SIZE DECILES ****************************************/
* omit zeros for computing deciles - only count active firms
replace size = . if size <= 0

* size decile by calendar year
forvalues i = 10(10)90 {
	bys year: egen size_p`i' = pctile(size), p(`i')
}
gen decile=1 if size<=size_p10 & !missing(size)
forvalues i = 2 (1) 9 {
	local j = `i'-1
	replace decile=`i' if size>size_p`j'0 & size<=size_p`i'0 & !missing(size)
}
replace decile=10 if size>size_p90 & !missing(size)
drop size_p*



/*** CREATE SIZE CLASSES, AGE CLASSES, AGExSIZE 'CELLS' ***********************/
gen size_class = round(decile)
gen age_class=.
replace age_class=1 if age==0
replace age_class=2 if age==1
replace age_class=3 if age==2
replace age_class=4 if age>=3 & age<=4
replace age_class=5 if age>=5 & age<=6
replace age_class=6 if age>=7 & age<=9
replace age_class=7 if age>=10 & age<=14
replace age_class=8 if age>=15 & age<=19
replace age_class=9 if age>=20 & age!=.
gen cell=0
replace cell=10*size_class+age_class /* encoding ageXsize class */



/********** COMPUTE STANDARDIZED SIZE (SALE SHARES FOR WEIGHTING) **************/
merge m:1 year using "$workingdir\raw\lagged_vars.dta"

* the standardized size cannot be negative, since it will be used as a weight later
* standardization method: each firm's standardized size is the average share of total revenue among all firms in the economy in the firm's year
* note that this doesn't standardize the variance of firm sizes within each year, but it's probably the best method 
gen size_sh = (1000 * size) / total_size



/**************** COMPUTE GROWTH RATES FORWARD ******************************/

* compute simple total growth rate of all sales
gen total_size_d = total_size / (cpi/(103.5083))
gen f_total_size_d = f_total_size / (f_cpi/(103.5083))

gen total_gro = (f_total_size_d - total_size_d) / total_size_d


****** Compute firm growth rates 4 different ways ****
* first deflate sales
gen size_d = size / (cpi/(103.5083))
gen f_size_d = f_size / (f_cpi/(103.5083))

* (1) simple growth rates
gen gro_simpl = (f_size_d - size_d) / size_d

* (2) log differences
gen gro_log = log(f_size_d) - log(size_d)

* (3) inverse hyperbolic sine
gen gro_ihs = log(f_size_d + sqrt(1+f_size_d^2)) - log(size_d + sqrt(1+size_d^2))
replace gro_ihs = . if size==0 & f_size==0

* (4) Haltiwanger, divide by average of current and future sales
gen gro_halt = (f_size_d - size_d) / ((f_size_d + size_d)/2)



/*** Add growth rates for entrants for IHS and Haltiwanger growth rates **/ 

expand 2 if age==0, gen(dup) // duplicate observation of new firms
sort firm_id year
replace age=-1 if dup == 1  // age=-1 means the firm had not entered yet
replace year=year-1 if dup == 1
replace gro_halt = 2 if age==-1 // here 'growth' is haltiwanger
replace gro_ihs = log(size_d + sqrt(1+size_d^2)) if age==-1 // here 'growth' is ihs, size_d is actually future size

** NOTE: the variable 'entry' is now doubled, make it so only age == -1 firms are entrants
replace entry = 0 if age == 0

* entry firms are going to all be in age_class 1 otherwise
replace age_class = 0 if entry == 1
*replace cell = . if entry == 1

** Now can drop 2018 because we needed it to construct 2017 entrants
* drop final year since we can't compute a total_gro for it
drop if year == 2018
drop if year == 1995



/****** STANDARDIZE GROWTH RATES OF SURVIVING FIRMS BY CALENDAR QUARTER *******/

* NOTE: we calculate mean and stdev of growth rates in that each quarter across all 
* surviving firms only (before they exit), because the forward growth rate in the last quarter just
* before exit either is equal to minus log size (ihs case) or is not defined (log(0) in the log case).
* We use those mean and stdev to standardize growth rates of all firms.

* Entrants are not an issue because growth rates are computed forward, so the first 
* growth rate of an entrant is between the first and second quarter of operations;
* if a firm enters for only one quarter, the growth rate is missing.

* After computing standardized growth rates for each calendar quarter, we average them 
* over the 4 quarters of a firm's age; so the 4 quarters of a firm of age 3 in 1998
* are associated to the 4 quarters of a firm of age 3 in 2007, because standardization
* kills time effects and leaves only age effects

* we don't remove exiting firms (ie. exit == 0)
foreach transf in "simpl" "log" "ihs" "halt" {
	gen std_gro_`transf' = gro_`transf' - total_gro
	bys firm_id age: egen stzd_gro_`transf' = mean(std_gro_`transf')
	drop std_gro_`transf'
}





/*** SAVE PANEL OF STANDARDIZED GROWTH RATES; DROP YEAR AS NO LONGER RELEVANT */
keep year firm_id age size stzd_gro* total_gro gro_simpl gro_log gro_ihs gro_halt size_class age_class cell size_sh entry exit 
tsset firm_id age
label var year "Year"
label var firm_id "Firm identifier (from anonymized VAT id)"

label var stzd_gro_ihs "Growth in IHS differences (standardized)"
label var stzd_gro_log "Growth in log differences (standardized)"
label var stzd_gro_simpl "Simple growth (standardized)"
label var stzd_gro_halt "Haltiwanger growth (standardized)"

label var total_gro "Total simple sales growth"
label var gro_simpl "Simple growth rate"
label var gro_log "Growth in log differences"
label var gro_ihs "Growth in ihs differences"
label var gro_halt "Haltiwanger growth rate"

label var exit "=1 in the year after the firm exits"
label var size_class "standardized decile of size distribution of active firms"
label var age "age of the firm (entry is at t=0)"
label var size_sh "average share of economywide firm revenue (x1000)"
label var age_class "age class"
label var cell "sizeXage cell (first digit size class, second digit age class)"
label var size "average quarterly sales"
* label var nq "number of quarters in this firm-age observation"


save "$workingdir\clean\new_vat_data_cleaned.dta", replace
* log close



