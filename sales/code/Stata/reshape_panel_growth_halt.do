clear
set more off
use D:\Kevin\Firm_Growth_Project\_Final\Sales\Data\clean\new_vat_data_cleaned.dta, clear

// we drop entry and exit since they have outlier values for growth
drop if entry == 1
drop if exit == 1

drop stzd_gro_log
rename stzd_gro_halt stzd_gro

keep firm_id age stzd_gro

local maxT = 100
replace age = age + 1
drop if age > `maxT'
sort firm_id age

reshape wide stzd_gro, i(firm_id) j(age)

// survivors
export delimited "D:\Kevin\Firm_Growth_Project\_Final\Sales\Data\clean\reshaped\new_vat_data_reshaped_growth_halt.csv", replace

// keep firm_id std_size_gro1-std_size_gro22

// export delimited "../../Data/clean/reshaped/new_vat_data_reshaped_surv_small.csv", replace

