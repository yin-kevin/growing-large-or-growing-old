/*
 
Trevor Williams, 7 September 2019. Revised by Giuseppe Moscarini, June 2021
Adjusted to use employment data by Kevin Yin, January 2022.

This do-file computes a panel of firms with age, standardized growth rate, and a standardized empl decile. 
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

3. empl_class: decile of the firm's standardized empl in that year 
(decile by quarter, averaged across quarters and rounded)

4. stzd_gro_ihs: annual forward growth rate (at t: between t and t+1) 
using the inverse hyperbolic sine transformation: 
g_t = IHS(s_t+1) - IHS(s_t) for each quarterly s
and then standardized by the mean and sdev of growth of all firms within the quarter.
Then the annual rate is the average quarterly rate

5. stdz_gro_empl_log: annual forward growth rate using the log transformation: 
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
global workingdir "D:\Kevin\Firm_Growth_Project\_Final\Employment\Data"

use "$workingdir\raw\data_empl.dta"

******** eliminate entries with age=0 for pre-entry info 

sort firm_id age year	
drop if firm_id==firm_id[_n+1]&age==0&age[_n+1]==0&size==0
tsset firm_id age



/**************** DETECT ENTRY AND EXIT ********/
sort firm_id year
gen entry = (firm_id~=firm_id[_n-1]) & (age==0) // dummy equal to 1 in first sample quarter the firm operates
gen exit = (firm_id~=firm_id[_n+1]) & (year<2018) // dummy equal to 1 in last year firm is in the data, if before 2018 - this will always be in the fourth calendar quarter

/**************** GENERATE EMPLOYMENT LEVEL VARIABLES ****************/
*** empl[i] is the annual average over quarters [i,i+3]
*** we want to set our employment level to empl[i] such that i is the starting quarter


* present employment level
gen empl = empl1 if startq == 1
replace empl = empl2 if startq == 2
replace empl = empl3 if startq == 3
replace empl = empl4 if startq == 4

* future employment level
gen f_empl = f_empl1 if startq == 1
replace f_empl = f_empl2 if startq == 2
replace f_empl = f_empl3 if startq == 3
replace f_empl = f_empl4 if startq == 4
replace f_empl = 0 if f_empl == .
replace exit = 1 if year == 2018 & f_empl == 0

* past employment level
gen l_empl = l_empl1 if startq == 1
replace l_empl = l_empl2 if startq == 2
replace l_empl = l_empl3 if startq == 3
replace l_empl = l_empl4 if startq == 4

label variable empl "average employment level for current year"
label variable f_empl "average employment level for future (next) year"
label variable l_empl "average employment level for last year"


****** When empl_lvl is missing the growth rate should be too. We leave f_empl=0 when the firm is exiting
* none of these make any changes
recode empl f_empl l_empl (-1 = .)
recode empl l_empl (0 = .)
replace f_empl=. if f_empl==0&exit==0



/**************** COMPUTE EMPLOYMENT BINS ****************************************/
* omit zeros for computing deciles - only count active firms
replace empl = . if empl <= 0

* instead of deciles, use employment bins defined as follows:
gen empl_class = 1
replace empl_class = 2 if empl > 4
replace empl_class = 3 if empl > 9
replace empl_class = 4 if empl > 19
replace empl_class = 5 if empl > 49
replace empl_class = 6 if empl > 99
replace empl_class = 7 if empl > 249
replace empl_class = 8 if empl > 999
replace empl_class = 9 if empl > 2499
replace empl_class = 10 if empl > 9999

label define empl_class_lab 1 "<= 4" 2 "(4,9]" 3 "(9,19]" 4 "(19,49]" 5 "(49,99]" 6 "(99,249]" 7 "(249,999]" 8 "(999,2499]" 9 "(2499,9999]" 10 "> 9999"
label values empl_class empl_class_lab

/***
* employment decile by calendar year
forvalues i = 10(10)90 {
	bys year: egen empl_p`i' = pctile(empl), p(`i')
}



gen decile=1 if empl<=empl_p10 & !missing(empl)
forvalues i = 2 (1) 9 {
	local j = `i'-1
	replace decile=`i' if empl>empl_p`j'0 & empl<=empl_p`i'0 & !missing(empl)
}



replace decile=10 if empl>empl_p90 & !missing(empl)
drop empl_p*

* 1 employee puts you in all of the first 3 deciles so we randomly split them up
gen random = runiform() if decile == 1
replace decile = 2 if rando > 1/3 & random != .
replace decile = 3 if rando > 2/3 & random != .

* the 5th decile is disproportionately small because it requires you have more than 2 but less than 2.25 so we balance it with the 4th
replace random = .
replace rando = runiform() if decile == 4
replace decile = 5 if rando < 0.35981 & random != .
***/


/*** CREATE EMPLOYMENT CLASSES, AGE CLASSES, AGExEMPLOYMENT 'CELLS' ***********************/
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
replace cell=10*empl_class+age_class /* encoding ageXempl class */



/********** COMPUTE STANDARDIZED EMPLOYMENT LEVEL (SALE SHARES FOR WEIGHTING) **************/

* the standardized empl cannot be negative, since it will be used as a weight later
* standardization method: each firm's standardized empl is the average share of total revenue among all firms in the economy in the firm's year
* note that this doesn't standardize the variance of firm empls within each year, but it's probably the best method 
bys year: egen total_empl = total(empl)
gen empl_sh = (1000 * empl) / total_empl



/**************** COMPUTE GROWTH RATES _FORWARD_ ******************************/

* quarterly growth rates (using unstandardized empls)


* compute growth rates using the inverse hyperbolic sine and the log

* inverse hyperbolic sine
gen gro_ihs = log(f_empl + sqrt(1+f_empl^2)) - log(empl + sqrt(1+empl^2))
replace gro_ihs = . if empl==0 & f_empl==0 // these are extra observations

* logarithm
gen gro_log = log(f_empl) - log(empl)

* haltiwanger
gen gro_halt = (f_empl - empl) / ((f_empl + empl)/2)


/****** STANDARDIZE GROWTH RATES OF SURVIVING FIRMS BY CALENDAR QUARTER *******/

* NOTE: we calculate mean and stdev of growth rates in that each quarter across all 
* surviving firms only (before they exit), because the forward growth rate in the last quarter just
* before exit either is equal to minus log empl (ihs case) or is not defined (log(0) in the log case).
* We use those mean and stdev to standardize growth rates of all firms.

* Entrants are not an issue because growth rates are computed forward, so the first 
* growth rate of an entrant is between the first and second quarter of operations;
* if a firm enters for only one quarter, the growth rate is missing.

* After computing standardized growth rates for each calendar quarter, we average them 
* over the 4 quarters of a firm's age; so the 4 quarters of a firm of age 3 in 1998
* are associated to the 4 quarters of a firm of age 3 in 2007, because standardization
* kills time effects and leaves only age effects

foreach transf in "ihs" "log" "halt" {
	bys year: egen mean_gro_`transf' = mean(gro_`transf')
	gen std_gro_`transf' = (gro_`transf' - mean_gro_`transf')
	bys firm_id age: egen stzd_gro_`transf' = mean(std_gro_`transf')
	drop mean_gro_`transf' std_gro_`transf'
}

/************ COMPUTE MEAN AND VARIANCE OF LOG SALES OF ENTRANTS **************/



/*** SAVE PANEL OF STANDARDIZED GROWTH RATES; DROP YEAR AS NO LONGER RELEVANT */
keep year firm_id age empl stzd_gro* empl_class age_class cell empl_sh entry exit 
tsset firm_id age
label var year "Year"
label var firm_id "Firm identifier (from anonymized VAT id)"
label var stzd_gro_ihs "Growth in IHS differences of employment (standardized)"
label var stzd_gro_log "Growth in log differences of employment (standardized)"
label var stzd_gro_halt "Growth in Haltiwanger rates of employment (standardized)"
label var exit "=1 in the year after the firm exits"
label var empl_class "classification bin of employment level"
label var age "age of the firm (entry is at t=0)"
label var empl_sh "average share of economywide firm employment (x1000)"
label var age_class "age class"
label var cell "emplXage cell (first digit empl class, second digit age class)"
label var empl "average annual employment level for current year"
* label var nq "number of quarters in this firm-age observation"


save "$workingdir\clean\empl_new_vat_data_cleaned.dta", replace
* log close



