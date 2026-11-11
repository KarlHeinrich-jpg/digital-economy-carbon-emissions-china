* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================================*
* 复合（地理邻接 × 经济距离）嵌套权重下的空间面板模型：SDM / SAR / SEM
*
* 权重构造（在已有三套权重基础上补充，不改动原始权重文件）：
*   A  : 重建后的 Queen 地理邻接二值矩阵（Wadj_bin.txt，
*        30×30 对称，132 条有向权重 / 66 对相邻省份，海南与广东、广西隔海补边）
*   E  : 经济距离逆矩阵 E_ij = 1/|mean(pgdp_i) − mean(pgdp_j)|（i≠j，不做 k 近邻截断）
*   复合权重 C_ij = A_ij × E_ij，对角线为 0，再按行标准化得到 W*
*        —— 邻居范围由地理邻接界定，邻居间权重按经济距离倒数赋权。
*
* 模型：xsmle，面板双向固定效应（省份 + 年份虚拟），稳健标准误，年份虚拟不滞后。
* 样本：30 省 × 2011—2022 年，N=360；控制变量 lnpgdp fisc dens fin urb mkt。
* 输出：本 do 运行日志（复合权重_SDM_SAR_SEM.log）；
*       表_复合权重空间面板结果.docx（putdocx 三线表）。
* 依赖：xsmle、estout（如未安装，见第 2 步自动安装）。
*==============================================================================*

clear all
set more off
set linesize 200

* 切到本 do 文件所在目录，保证相对路径可复现
local full "`c(filename)'"
local full : subinstr local full "\" "/" , all
local slash = strrpos("`full'", "/")
if `slash' > 0 {
    local dotdir = substr("`full'", 1, `slash'-1)
    capture cd "`dotdir'"
}

*---------------------------------------------------------------*
* 1. 定位面板数据（实验目录与归档目录两种命名均兼容）
*---------------------------------------------------------------*
global csv1 "../空间杜宾，滞后，误差/panel_AI5.csv"
global csv2 "../空间杜宾，滞后/panel_AI5.csv"
capture confirm file "$csv1"
if !_rc {
    global csv "$csv1"
}
else {
    capture confirm file "$csv2"
    if !_rc global csv "$csv2"
    else display as error "找不到panel_AI5.csv"
}
global adj "Wadj_bin.txt"

if "$csv"=="" {
    display as error "找不到panel_AI5.csv，终止"
    exit
}

import delimited using "$csv", encoding(UTF-8) varnames(1) case(preserve) clear
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

global controls lnpgdp fisc dens fin urb mkt

* 年份虚拟变量（不进入空间滞后项）
capture drop yr*
forvalues y = 2012/2022 {
    qui gen yr`y' = (year == `y')
}
global yrdum yr2012 yr2013 yr2014 yr2015 yr2016 yr2017 yr2018 ///
              yr2019 yr2020 yr2021 yr2022

*---------------------------------------------------------------*
* 2. 依赖包
*---------------------------------------------------------------*
capture which xsmle
if _rc ssc install xsmle, replace
capture which esttab
if _rc ssc install estout, replace

*---------------------------------------------------------------*
* 3. 构造复合嵌套权重 W*
*    A：二值邻接（读入已有 Wadj_bin.txt，行序 = ord 1..30）
*    E：人均 GDP 均值经济距离倒数
*    C = A ⊙ E，行标准化
*---------------------------------------------------------------*
preserve
keep ord pgdp
collapse (mean) pgdp, by(ord)
sort ord
mata:
    n = st_nobs()
    y = st_data(., "pgdp")
    E = J(n, n, 0)
    for (i=1; i<=n; i++) {
        for (j=1; j<=n; j++) {
            if (i!=j) E[i,j] = 1/abs(y[i]-y[j])
        }
    }
    adjpath = st_global("adj")
    fh = fopen(adjpath, "r")
    A  = J(n, n, 0)
    for (i=1; i<=n; i++) {
        toks = tokens(fget(fh))
        v = strtoreal(toks)
        for (j=1; j<=cols(toks); j++) A[i,j] = v[j]
    }
    fclose(fh)
    C = A :* E
    rs = rowsum(C)
    W  = C :/ rs
    st_matrix("W", W)
    printf("复合权重：非零元素 = %g；行和最大偏差 = %g；对称 = %g\n", ///
        sum(C:!=0), max(abs(rowsum(W):-1)), allof((W:-W'):==0, 1))
end
restore

* 矩阵指纹校验
matrix Wcheck = W
mata:
    M = st_matrix("Wcheck")
    d = rowsum(M:!=0)
    printf("每行邻居数：最少 %g 个，最多 %g 个\n", min(d), max(d))
end

*---------------------------------------------------------------*
* 4. 三大空间模型（面板双向固定效应；年份虚拟变量不滞后）
*---------------------------------------------------------------*
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
* 5. LR 检验（SDM vs SAR / SDM vs SEM）
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
di _n "===== LR 检验（复合权重） ====="
di "SDM vs SAR: LR = " lr_sar "  df = " k_durbin "  p = " p_sar
di "SDM vs SEM: LR = " lr_sem "  df = " k_durbin "  p = " p_sem

*---------------------------------------------------------------*
* 6. 提取结果（供 putdocx 三线表使用）
*---------------------------------------------------------------*
program define starof, rclass
    args z
    local s ""
    if abs(`z') > 2.576 local s "***"
    else if abs(`z') > 1.96  local s "**"
    else if abs(`z') > 1.645 local s "*"
    return local s "`s'"
end

* 取某模型某方程某变量的 系数(标准误) 字符串；缺失返回 —
program define cellof, rclass
    args model eq var
    est restore `model'
    matrix b  = e(b)
    matrix V  = e(V)
    local cn : colnames b
    local ce : coleq b
    local k = colsof(b)
    local found = 0
    forvalues i = 1/`k' {
        local eqi : word `i' of `ce'
        local nmi : word `i' of `cn'
        if ("`eqi'"=="`eq'" & "`nmi'"=="`var'") {
            local bb = b[1,`i']
            local ss = sqrt(V[`i',`i'])
            local zz = `bb'/`ss'
            starof `zz'
            local str = string(`bb',"%9.3f") + "`r(s)'" + " (" + string(`ss',"%9.3f") + ")"
            local found = 1
            continue, break
        }
    }
    if `found'==0 local str "—"
    return local str "`str'"
end

*---------------------------------------------------------------*
* 7. Word 三线表：表_复合权重空间面板结果.docx
*---------------------------------------------------------------*
local vars "lnDEI lnpgdp fisc dens fin urb mkt"
* 数据行：Panel A 1标签+7；Panel B 1标签+7；Panel C 1标签+3；空间参数 1；拟合指标 3
local nrow = 8 + 8 + 4 + 1 + 3
putdocx begin
putdocx paragraph, halign(center)
putdocx text ("表  复合权重（地理邻接 × 经济距离）下空间面板模型估计结果"), bold

putdocx table t = (`=`nrow'+1', 4), border(all, nil) border(top, single) border(bottom, single) ///
    width(92%) layout(autofitcontents)
putdocx table t(1,1) = ("变量"), halign(center) bold
putdocx table t(1,2) = ("SDM"), halign(center) bold
putdocx table t(1,3) = ("SAR"), halign(center) bold
putdocx table t(1,4) = ("SEM"), halign(center) bold
putdocx table t(1,.), border(bottom, single)

local r = 1

* Panel A：主回归系数
putdocx table t(`++r',1) = ("Panel A：主回归系数 β"), bold
foreach v of local vars {
    local ++r
    cellof sdm Main `v'
    local a "`r(str)'"
    cellof sar Main `v'
    local b "`r(str)'"
    cellof sem Main `v'
    local c "`r(str)'"
    putdocx table t(`r',1) = ("`v'"), halign(center)
    putdocx table t(`r',2) = ("`a'"), halign(center)
    putdocx table t(`r',3) = ("`b'"), halign(center)
    putdocx table t(`r',4) = ("`c'"), halign(center)
}

* Panel B：解释变量空间滞后项 W·X（仅 SDM）
putdocx table t(`++r',1) = ("Panel B：空间滞后项 W·X（仅 SDM）"), bold
foreach v of local vars {
    local ++r
    cellof sdm Wx `v'
    local a "`r(str)'"
    putdocx table t(`r',1) = ("W·`v'"), halign(center)
    putdocx table t(`r',2) = ("`a'"), halign(center)
    putdocx table t(`r',3) = ("—"), halign(center)
    putdocx table t(`r',4) = ("—"), halign(center)
}

* Panel C：lnDEI 长期效应分解
putdocx table t(`++r',1) = ("Panel C：lnDEI 长期效应分解"), bold
foreach eq in LR_Direct LR_Indirect LR_Total {
    local ++r
    cellof sdm `eq' lnDEI
    local a "`r(str)'"
    cellof sar `eq' lnDEI
    local b "`r(str)'"
    local lab : di cond("`eq'"=="LR_Direct","直接效应", ///
                 cond("`eq'"=="LR_Indirect","间接效应","总效应"))
    putdocx table t(`r',1) = ("`lab'"), halign(center)
    putdocx table t(`r',2) = ("`a'"), halign(center)
    putdocx table t(`r',3) = ("`b'"), halign(center)
    putdocx table t(`r',4) = ("—"), halign(center)
}

* 空间参数 ρ / λ
local ++r
cellof sdm Spatial rho
local ar "`r(str)'"
cellof sar Spatial rho
local br "`r(str)'"
cellof sem Spatial lambda
local cl "`r(str)'"
putdocx table t(`r',1) = ("空间参数 ρ / λ"), halign(center)
putdocx table t(`r',2) = ("`ar'"), halign(center)
putdocx table t(`r',3) = ("`br'"), halign(center)
putdocx table t(`r',4) = ("`cl'"), halign(center)

* 拟合指标
est restore sdm
local ll_sdm = e(ll)
local r2_sdm = e(r2_w)
local N      = e(N)
est restore sar
local ll_sar = e(ll)
local r2_sar = e(r2_w)
est restore sem
local ll_sem = e(ll)
local r2_sem = e(r2_w)

local ++r
putdocx table t(`r',1) = ("对数似然 LL"), halign(center)
putdocx table t(`r',2) = ("`=string(`ll_sdm',"%9.3f")'"), halign(center)
putdocx table t(`r',3) = ("`=string(`ll_sar',"%9.3f")'"), halign(center)
putdocx table t(`r',4) = ("`=string(`ll_sem',"%9.3f")'"), halign(center)
local ++r
putdocx table t(`r',1) = ("组内 R²"), halign(center)
putdocx table t(`r',2) = ("`=string(`r2_sdm',"%9.3f")'"), halign(center)
putdocx table t(`r',3) = ("`=string(`r2_sar',"%9.3f")'"), halign(center)
putdocx table t(`r',4) = ("`=string(`r2_sem',"%9.3f")'"), halign(center)
local ++r
putdocx table t(`r',1) = ("观测值 N"), halign(center)
putdocx table t(`r',2) = ("`=string(`N',"%6.0f")'"), halign(center)
putdocx table t(`r',3) = ("`=string(`N',"%6.0f")'"), halign(center)
putdocx table t(`r',4) = ("`=string(`N',"%6.0f")'"), halign(center)

putdocx paragraph
putdocx text ("注：被解释变量为 lnCEI；括号内为稳健标准误；***、**、* 分别表示 1%、5%、10% 水平显著。"), italic
putdocx text ("样本为 30 省 2011—2022 年面板（N=360），控制省份与年份固定效应及 lnpgdp、fisc、dens、fin、urb、mkt。"), italic
putdocx text ("复合权重 W* = A（Queen 邻接）⊙ E（1/|人均 GDP 均值差|），行标准化；SAR、SEM 不含解释变量空间滞后项，"), italic
putdocx text ("SEM 不提供直接/间接/总效应分解；Panel C 为长期（long-run）效应分解。"), italic
putdocx save "表_复合权重空间面板结果.docx", replace

eststo clear
display as result "完成：复合权重_SDM_SAR_SEM.do 已输出 表_复合权重空间面板结果.docx"