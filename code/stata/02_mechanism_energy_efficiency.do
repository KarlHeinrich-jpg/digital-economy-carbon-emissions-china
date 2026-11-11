* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================*
* 能源效率(lnee)中介效应三步法，M取对数
*   (1) 总效应 c ：lnCEI ~ lnDEI + 6控制
*   (2) 直接效应 c' + b：lnCEI ~ lnDEI lnee + 6控制
*   (3) a 路径：lnee ~ lnDEI + 6控制
*==============================================================*
clear all
set more off
capture log close

cd "`c(pwd)'"

import delimited using "panel_AI5.csv", encoding(UTF-8) case(preserve) clear

unab allvars : *
local yearvar : word 2 of `allvars'
rename `yearvar' year

encode sh, gen(id)
gen lnCEI = log(CEI)
gen lnDEI = log(DEI)
gen lnee = log(ee)
xtset id year

global controls lnpgdp fisc dens fin urb mkt

cap which esttab
if _rc ssc install estout, replace

eststo clear

eststo c:  reg lnCEI lnDEI $controls       i.id i.year, vce(cluster id)
eststo cp: reg lnCEI lnDEI lnee $controls i.id i.year, vce(cluster id)
eststo a:  reg lnee lnDEI $controls     i.id i.year, vce(cluster id)

esttab c cp a, b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI lnee $controls) order(lnDEI lnee) ///
    stats(N r2 r2_a, labels("N" "R2" "adj_R2") fmt(%9.0f %9.3f %9.3f)) ///
    mtitles("(1) 总效应c" "(2) 直接效应c'+b" "(3) a路径") nonotes

esttab c cp a using "三步法_lnee.rtf", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI lnee $controls) order(lnDEI lnee) ///
    stats(N r2 r2_a, labels("N" "R2" "adj_R2") fmt(%9.0f %9.3f %9.3f))

eststo clear