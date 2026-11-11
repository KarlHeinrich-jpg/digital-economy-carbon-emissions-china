* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================*
* 表1  DEI 对 CEI 影响的基准回归 —— 已在 Stata18 实跑验证
* 数据：panel_AI5.csv（30 省 x 2011-2022，N=360）
*
* 【口径说明】双向固定效应采用 LSDV 写法：reg y x i.id i.year
*   - 其系数与 xtreg,fe 完全相同；
*   - 但文档报告的聚类标准误(0.114/0.092/0.094/0.086)、整体R2(0.98x)、
*     调整R2(0.97x) 正是 reg+i.id i.year 这一口径，已逐位对上文档；
*   - 若改用 xtreg ...,fe，系数不变，但 SE 与 R2 报告口径会与文档略异。
*==============================================================*
clear all
set more off
capture log close

cd "`c(pwd)'"

* 导入；case(preserve) 保留 CEI/DEI 大写，否则会被转成小写
import delimited using "panel_AI5.csv", encoding(UTF-8) case(preserve) clear

* 第二列是中文"年份"，Stata 变量名不允许中文，按位置重命名为 year
unab allvars : *
local yearvar : word 2 of `allvars'
rename `yearvar' year

encode sh, gen(id)
gen lnCEI = log(CEI)
gen lnDEI = log(DEI)

xtset id year

* 全局控制变量（已锁定）
global controls lnpgdp fisc dens fin urb mkt

cap which esttab
if _rc ssc install estout, replace

eststo clear

* Panel A：被解释变量 lnCEI（对数）
eststo m1: reg lnCEI lnDEI,                       vce(cluster id)
eststo m2: reg lnCEI lnDEI           i.id i.year,  vce(cluster id)
eststo m3: reg lnCEI lnDEI $controls i.id i.year,  vce(cluster id)

* Panel B：被解释变量 CEI（水平值）
eststo m4: reg CEI DEI,                           vce(cluster id)
eststo m5: reg CEI DEI            i.id i.year,     vce(cluster id)
eststo m6: reg CEI DEI $controls  i.id i.year,     vce(cluster id)

* 屏幕显示
esttab m1 m2 m3 m4 m5 m6, b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI DEI $controls) ///
    stats(N r2 r2_a, labels("样本量 N" "R2 整体" "调整 R2") fmt(%9.0f %9.3f %9.3f)) ///
    mtitles("(1)" "(2)" "(3)" "(4)" "(5)" "(6)") nonotes

* 导出 RTF（Word 可打开）
esttab m1 m2 m3 m4 m5 m6 using "表1_基准回归_复现.rtf", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI DEI $controls) ///
    stats(N r2 r2_a, labels("N" "R2" "adj_R2") fmt(%9.0f %9.3f %9.3f))

eststo clear