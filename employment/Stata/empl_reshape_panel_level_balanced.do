clear
set more off
use D:\Kevin\Firm_Growth_Project\_Final\Employment\Data\clean\empl_new_vat_data_cleaned.dta, clear

gen log_size=log(empl)

// drop if exit==1
keep firm_id age log_size



replace age = age + 1

// Keep only firms that have reached age 20

local maxage=20
bys firm_id: egen maxagefirm=max(age)
bys firm_id: gen is_balanced = (maxagefirm>=`maxage')
drop maxagefirm
drop if is_balanced==0
drop is_balanced
 
local maxT = 100
drop if age > `maxT'
sort firm_id age 

reshape wide log_size, i(firm_id) j(age)

// for balanced data, log_sizes 1-5 don't exist, generate them
/*
forvalues i = 1/5 {
	gen log_size`i' = .
}
*/

// survivors
export delimited "D:\Kevin\Firm_Growth_Project\_Final\Employment\Data\clean\reshaped\empl_new_vat_data_reshaped_balanced.csv", replace

//keep firm_id std_size_gro1-std_size_gro22

//export delimited "../../Data/clean/reshaped/new_vat_data_reshaped_surv_small_level.csv", replace

