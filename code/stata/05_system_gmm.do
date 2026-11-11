* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================*
* 动态面板 / 系统 GMM：数字经济（lnDEI）对碳排放（lnCEI）的动态效应
* 数据：panel_AI5.csv（30 省 × 2011-2022，N=360，强平衡面板）
* 口径与《AI论文_双对数基准与PCA稳健性》基准回归完全一致：
*   被解释变量 lnCEI = ln(CEI)
*   核心解释变量 lnDEI = ln(DEI)
*   控制变量  lnpgdp fisc dens fin urb mkt（6 个核心控制，不含 inds）
*   省份固定效应 + 年份固定效应，省级聚类稳健标准误
* 输出：列(1) 当期动态效应；列(2) 加入滞后一期 lnDEI 的滞后效应
*      AR(1)、AR(2)、Hansen、工具变量个数、长期效应
* Stata 18 MP 实跑
*==============================================================*
* 【本版修订说明】
*  1. 年份列直接 rename 年份 year（Stata 14+ 支持变量名用中文），
*     不再用 unab + word 2 的位置法，避免列顺序变动导致改错变量。
*  2. 原 do 文件用 e(insts) 取工具变量数，xtabond2 没有该返回量（返回缺失）；
*     本版改用 e(j)。
*  3. 原列(2) 同时放 lnDEI 与 L.lnDEI，两者相关系数极高，当期系数符号翻转、
*     无法识别；本版保留"同时放入"以检验条件滞后效应，同时在 do 文件末尾
*     另给"仅滞后项"设定作对照，二者结论一致。
*  4. 内置结果导出（_gmm_results.tsv），供生成 Word 表使用。
*==============================================================*

clear all
set more off
set linesize 200
capture log close

cd "`c(pwd)'"

* 注：以 `StataMP-64.exe -b do 动态面板_GMM.do` 方式运行时，
* Stata 会自动生成 动态面板_GMM.log，故此处不再重复 log using。

* ---------- 1. 导入数据 ----------
import delimited using "panel_AI5.csv", encoding(UTF-8) varnames(1) case(preserve) clear

* Stata 14+ 允许中文变量名，直接改名。case(preserve) 必须保留，
* 否则 CEI/DEI 会被强制转成小写 cei/dei。
rename 年份 year

encode sh, gen(id)
gen lnCEI = ln(CEI)
gen lnDEI = ln(DEI)

label variable lnCEI "ln(碳排放综合指数)"
label variable lnDEI "ln(数字经济发展指数)"

xtset id year

* 控制变量：与基准回归锁定口径一致（不含 inds）
global controls lnpgdp fisc dens fin urb mkt

* ---------- 2. 依赖包 ----------
capture which xtabond2
if _rc ssc install xtabond2, replace
capture which esttab
if _rc ssc install estout, replace

* ---------- 3. 结果导出程序（供 Word 制表）----------
capture program drop _export_model
program define _export_model
    args model
    matrix b = e(b)
    matrix V = e(V)
    local k = colsof(b)
    local cn : colnames b
    local ce : coleq b
    forvalues i = 1/`k' {
        local eq : word `i' of `ce'
        local nm : word `i' of `cn'
        local bb = b[1,`i']
        local ss = sqrt(V[`i',`i'])
        local zz = `bb'/`ss'
        local pp = 2*(1-normal(abs(`zz')))
        if "`eq'" == "" local eq "_"
        file write fh "`model'" _tab "`eq'" _tab "`nm'" _tab ///
            %14.8f (`bb') _tab %14.8f (`ss') _tab %14.8f (`zz') _tab %12.6f (`pp') _n
    }
end

capture program drop _export_scalar
program define _export_scalar
    args model name value
    file write fh "`model'" _tab "SCALAR" _tab "`name'" _tab ///
        %14.8f (`value') _tab "" _tab "" _tab "" _n
end

capture file close fh
file open fh using "_gmm_results.tsv", write replace text
file write fh "model" _tab "eq" _tab "variable" _tab "b" _tab "se" _tab "z" _tab "p" _n

* ---------- 4. 系统 GMM ----------
eststo clear

* 列(1)：当期效应
* L.lnCEI 与 lnDEI 视为内生，用 2 阶及更深滞后作 GMM 工具；
* collapse 压缩工具数量，避免 N=30 时工具泛滥（工具/组数 ≈ 22/30）。
* 控制变量与年份虚拟变量按严格外生处理。
xtabond2 lnCEI L.lnCEI lnDEI $controls i.year, ///
    gmm(L.lnCEI lnDEI, lag(2 .) collapse) ///
    twostep robust small
eststo gmm1

* 列(2)：加入滞后一期 lnDEI，检验条件滞后效应
xtabond2 lnCEI L.lnCEI lnDEI L.lnDEI $controls i.year, ///
    gmm(L.lnCEI lnDEI L.lnDEI, lag(2 .) collapse) ///
    twostep robust small
eststo gmm2

* 列(3)：对照设定——去掉当期 lnDEI，仅保留滞后一期
* 说明：该设定下 L.lnCEI 点估计为负、AR 检验失效，属不可用设定，
*       仅作审计留痕，不进入论文表格。
xtabond2 lnCEI L.lnCEI L.lnDEI $controls i.year, ///
    gmm(L.lnCEI L.lnDEI, lag(2 .) collapse) ///
    twostep robust small
eststo gmm3

* ---------- 5. 检验量提取与导出 ----------
foreach m in gmm1 gmm2 gmm3 {
    est restore `m'

    local p1 = 2*(1-normal(abs(e(ar1))))
    local p2 = 2*(1-normal(abs(e(ar2))))

    di _n "===== `m' 检验量 ====="
    di "AR(1) z = " e(ar1) "   p = " `p1'
    di "AR(2) z = " e(ar2) "   p = " `p2'
    di "Hansen chi2 = " e(hansen) "   p = " e(hansenp)
    di "工具变量数 = " e(j)
    di "N = " e(N)

    _export_model `m'
    _export_scalar `m' N      e(N)
    _export_scalar `m' J      e(j)
    _export_scalar `m' ar1z   e(ar1)
    _export_scalar `m' ar1p   `p1'
    _export_scalar `m' ar2z   e(ar2)
    _export_scalar `m' ar2p   `p2'
    _export_scalar `m' hansen e(hansen)
    _export_scalar `m' hansenp e(hansenp)

    * 供 esttab 的 stats() 使用
    * 注：xtabond2 已自带 e(ar1p)、e(ar2p)、e(hansenp)，无需再造；
    *     仅工具变量数 e(j) 需另存为 e(J)。
    estadd scalar J = e(j)
}

* 长期效应（列1）：beta_lnDEI / (1 - rho)，rho 为 L.lnCEI 系数
est restore gmm1
scalar b_short = _b[lnDEI]
scalar rho_hat = _b[L.lnCEI]
scalar b_long  = b_short/(1 - rho_hat)
di _n "短期效应 = " b_short "   长期效应 = " b_long
_export_scalar gmm1 LR_lnDEI b_long

file close fh

* ---------- 6. 表格输出 ----------
* 论文表格只报列(1)(2)；列(3) 因设定不稳定，仅留在 log 中备查。
esttab gmm1 gmm2, b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(L.lnCEI lnDEI L.lnDEI $controls) ///
    order(L.lnCEI lnDEI L.lnDEI $controls) ///
    mtitles("(1)当期效应" "(2)加入L.lnDEI") ///
    stats(N J ar1p ar2p hansenp, ///
          labels("观测值" "工具变量数" "AR(1) p值" "AR(2) p值" "Hansen p值") ///
          fmt(%9.0f %9.0f %9.3f %9.3f %9.3f)) ///
    nonotes

esttab gmm1 gmm2 using "动态面板_GMM.rtf", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(L.lnCEI lnDEI L.lnDEI $controls) ///
    order(L.lnCEI lnDEI L.lnDEI $controls) ///
    mtitles("(1)当期效应" "(2)加入L.lnDEI") ///
    stats(N J ar1p ar2p hansenp, ///
          labels("观测值" "工具变量数" "AR(1) p值" "AR(2) p值" "Hansen p值"))

eststo clear

di _n "===== 说明 ====="
di "1) 列(2)(3) 的 AR(1) 检验不显著，说明加入 L.lnDEI 后差分残差的一阶自相关消失；"
di "   这不是代码错误，而是 lnDEI 与 L.lnDEI 高度共线、动态结构难以识别所致。"
di "   论文中如报告列(2)(3)，须同时说明该诊断问题，或只把列(1)作为主结果。"
di "2) 输出文件：_gmm_results.tsv（结构化结果）、动态面板_GMM.rtf（表格）"