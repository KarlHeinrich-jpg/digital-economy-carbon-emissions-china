* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================================*
* 空间计量：SDM（空间杜宾）/ SAR（空间滞后）/ SEM（空间误差）
*
* 版本：v2（2026-11-10 统一空间权重口径）
* 本版相对旧版的三处实质修改：
*   1. 省名补齐：原始数据“内蒙”“黑龙”为截断简称，旧版只进矩阵 28 个省；
*      本版统一为“内蒙古”“黑龙江”，严格 30 省 × 12 年 = 360 观测。
*   2. 权重口径统一：与“空间相关性_Stata”文件夹完全一致——
*        W_ij = 1/|mean(pgdp_i) − mean(pgdp_j)|（人均 GDP 原值，不是 lnpgdp）
*        每省只保留最近 5 个邻居，再行标准化
*   3. 行序：按与权重矩阵一致的省级顺序显式编号 ord=1..30 再 xtset，
*      避免 Stata 对中文自动排序与矩阵行列错配。
*
* 数据：panel_AI5.csv（30 省 × 2011-2022）
* 依赖：xsmle、estout（如未安装，见第 2 步自动安装）
* 输出：_spatial_results.tsv、空间计量_SDM_SAR_SEM.rtf
*
* 矩阵可直接复用“FWecon.dta —> 但本文件不依赖外部文件，
* 而是在 Stata 内按同一规则现算，便于单独运行与复核。
*==============================================================================*

clear all
set more off
set linesize 200
cd "`c(pwd)'"

*---------------------------------------------------------------*
* 1. 数据
*---------------------------------------------------------------*
import delimited using "panel_AI5.csv", encoding(UTF-8) varnames(1) case(preserve) clear
rename 年份 year

replace sh = "内蒙古" if sh=="内蒙"
replace sh = "黑龙江" if sh=="黑龙"

gen double lnCEI = ln(CEI)
gen double lnDEI = ln(DEI)
label variable lnCEI "ln(碳排放综合指数)"
label variable lnDEI "ln(数字经济发展指数)"

* 省级顺序（与权重矩阵行列一致）
gen long ord = .
local i = 0
foreach p in 北京 天津 河北 山西 内蒙古 辽宁 吉林 黑龙江 上海 江苏 ///
             浙江 安徽 福建 江西 山东 河南 湖北 湖南 广东 广西 ///
             海南 重庆 四川 贵州 云南 陕西 甘肃 青海 宁夏 新疆 {
    local ++i
    replace ord = `i' if sh=="`p'"
}
count if missing(ord)
assert r(N)==0
xtset ord year
levelsof ord, local(ols)
di _n "省份数 = " wordcount("`ols'") "（应为 30）"

global controls lnpgdp fisc dens fin urb mkt

*---------------------------------------------------------------*
* 2. 依赖包
*---------------------------------------------------------------*
capture which xsmle
if _rc ssc install xsmle, replace
capture which esttab
if _rc ssc install estout, replace

*---------------------------------------------------------------*
* 3. 经济距离矩阵（与空间相关性_Stata 完全同口径）
*    W_ij = 1/|mean(pgdp_i) − mean(pgdp_j)|，保留最近 5 个邻居后行标准化
*---------------------------------------------------------------*
preserve
keep ord pgdp
collapse (mean) pgdp, by(ord)
sort ord
local n = _N
matrix W = J(`n', `n', 0)
forvalues i = 1/`n' {
    forvalues j = 1/`n' {
        if `i' != `j' {
            matrix W[`i', `j'] = 1 / abs(pgdp[`i'] - pgdp[`j'])
        }
    }
}
mata:
    M    = st_matrix("W")
    n    = rows(M)
    K    = 5
    Mtr  = M
    for (i = 1; i <= n; i++) {
        row    = M[i, .]
        row[i] = .
        idx    = order(row', 1)
        keep   = idx[1..K]
        newrow = J(1, n, 0)
        for (k = 1; k <= K; k++) newrow[keep[k]] = row[keep[k]]
        Mtr[i, .] = newrow
    }
    Mtr = Mtr :/ rowsum(Mtr)
    st_matrix("W", Mtr)
    st_numscalar("w_rowdev", max(abs(rowsum(Mtr) :- 1)))
    st_numscalar("w_nonzero", sum(Mtr :!= 0))
    st_numscalar("w_max", max(Mtr))
    st_numscalar("w_nbig", sum(rowmax(Mtr) :>= 0.5))
end
restore

di _n "===== 权重矩阵校验 ====="
di "维度 n              = `n'（应为 30）"
di "非零元素数          = " w_nonzero "（应为 150 = 30×5）"
di "行和最大偏差        = " w_rowdev
di "最大元素            = " w_max
di "行首要权重≥0.5 的行数 = " w_nbig

*---------------------------------------------------------------*
* 4. 结果导出程序
*---------------------------------------------------------------*
capture program drop _export_model
program define _export_model
    args model
    matrix b  = e(b)
    matrix V  = e(V)
    local k  = colsof(b)
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
file open fh using "_spatial_results.tsv", write replace text
file write fh "model" _tab "eq" _tab "variable" _tab "b" _tab "se" _tab "z" _tab "p" _n

*---------------------------------------------------------------*
* 5. 三大空间模型（面板双向固定效应；年份虚拟变量不滞后）
*---------------------------------------------------------------*
forvalues y = 2012/2022 {
    qui gen yr`y' = (year == `y')
}
global yrdum yr2012 yr2013 yr2014 yr2015 yr2016 yr2017 yr2018 ///
              yr2019 yr2020 yr2021 yr2022

eststo clear

* SDM：核心解释变量与控制变量同时空间滞后
xsmle lnCEI lnDEI $controls $yrdum, wmat(W) model(sdm) fe robust nolog ///
    effects durbin(lnDEI $controls)
eststo sdm

* SAR：仅 W*y
xsmle lnCEI lnDEI $controls $yrdum, wmat(W) model(sar) fe robust nolog effects
eststo sar

* SEM：空间相关在扰动项
xsmle lnCEI lnDEI $controls $yrdum, ematrix(W) model(sem) fe robust nolog
eststo sem

*---------------------------------------------------------------*
* 6. 导出
*---------------------------------------------------------------*
foreach m in sdm sar sem {
    est restore `m'
    di _n "===== `m' ====="
    di "log-likelihood = " e(ll)
    di "N = " e(N)
    _export_model `m'
    _export_scalar `m' ll e(ll)
    _export_scalar `m' N  e(N)
}

*---------------------------------------------------------------*
* 7. LR 检验（SDM vs SAR / SDM vs SEM）
*---------------------------------------------------------------*
est restore sdm
scalar ll_sdm = e(ll)
est restore sar
scalar ll_sar = e(ll)
est restore sem
scalar ll_sem = e(ll)
scalar k_durbin = wordcount("$controls") + 1
scalar lr_sar = 2*(ll_sdm - ll_sar)
scalar p_sar  = chi2tail(k_durbin, lr_sar)
scalar lr_sem = 2*(ll_sdm - ll_sem)
scalar p_sem  = chi2tail(k_durbin, lr_sem)
di _n "===== LR 检验 ====="
di "SDM vs SAR: LR = " lr_sar "  df = " k_durbin "  p = " p_sar
di "SDM vs SEM: LR = " lr_sem "  df = " k_durbin "  p = " p_sem
_export_scalar test LR_SDM_SAR   lr_sar
_export_scalar test LR_SDM_SAR_p p_sar
_export_scalar test LR_SDM_SEM   lr_sem
_export_scalar test LR_SDM_SEM_p p_sem
file close fh

*---------------------------------------------------------------*
* 8. 导出表格
*---------------------------------------------------------------*
esttab sdm sar sem using "空间计量_SDM_SAR_SEM.rtf", replace ///
    b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI $controls) order(lnDEI $controls) ///
    mtitles("SDM" "SAR" "SEM") ///
    stats(ll N, labels("对数似然" "观测值") fmt(%9.3f %9.0f)) nonotes

esttab sdm sar sem, b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(lnDEI $controls) order(lnDEI $controls) ///
    mtitles("SDM" "SAR" "SEM") stats(N ll)

eststo clear
display as result "完成：_spatial_results.tsv / 空间计量_SDM_SAR_SEM.rtf"