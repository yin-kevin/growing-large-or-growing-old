clear
set more off
use D:\Kevin\Firm_Growth_Project\_Final\Sales\Data\clean\new_vat_data_cleaned.dta, clear

drop if entry == 1
drop if exit == 1

gen log_size=log(size)

// drop if exit==1 // no need to do this step for levels

keep firm_id age log_size
local maxT = 100
replace age = age + 1
drop if age > `maxT'
sort firm_id age 

reshape wide log_size, i(firm_id) j(age)

// survivors
export delimited "D:\Kevin\Firm_Growth_Project\_Final\Sales\Data\clean\reshaped\new_vat_data_reshaped_level_unbalanced.csv", replace

//keep firm_id std_size_gro1-std_size_gro22

//export delimited "../../Data/clean/reshaped/new_vat_data_reshaped_surv_small_level.csv", replace

