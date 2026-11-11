* Repository note: this do-file was archived from the working project and
* adapted for this repository. Run it with the working directory set to
* the repository's data/ folder, or edit the cd line below.
* The main panel is data/panel_AI5.csv (originally 面板_AI5(1).csv).
* Stata 17 or later is required; spatial models use xsmle/spatwmat.

*==============================================================================*
* 1_wmat_moran.do
* 空间权重矩阵（spatwmat）+ 全局空间自相关检验（spatgsa：Moran's I / Geary's c）
*
* 数据：panel.txt（30 省 × 2011-2022，空格分隔，无表头）
* 矩阵：FWadj.dta / FWecon.dta / FWdist.dta（30 × 30，列名 w01..w30，行标准化）
* 环境：Stata 18 MP + sg162（Pisati, 2001）
*       net install sg162, from("http://www.stata.com/stb/stb60") replace
*
* 用法：把本 do 文件与上述数据放在同一文件夹，Stata 中 cd 到该文件夹后
*       do 1_wmat_moran.do
*
* 三套权重矩阵（均已行标准化）：
*   Wadj  地理邻接：Queen 邻接，海南按与广东、广西隔海相邻处理
*   Wdist 地理距离：W_ij = 1/d_ij，省会大圆距离，保留最近 5 个邻居
*   Wecon 经济距离：W_ij = 1/|mean(pgdp_i) − mean(pgdp_j)|，保留最近 5 个邻居
*
* 输出：
*   moran_long.txt   逐年 Moran's I 长表：year W Wname variable I EI sd z p
*   geary_long.txt   逐年 Geary's c 长表
*   moran_width.txt  论文用宽表（DEI、CEI × 三个矩阵的 I 与 z）
*==============================================================================*

clear all
set more off
set linesize 200

*---------------------------------------------------------------*
* 0. 面板数据
*---------------------------------------------------------------*
qui import delimited using "panel.txt", delimiter(" ") clear ///
    varnames(nonames) stringcols(1 3 4)
qui rename (v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 ///
            v16 v17 v18 v19 v20 v21 v22) ///
           (province year region inc_grp trend DEI CEI CE lnCE ///
            pgdp lnpgdp taxb fisc conslv dens fin urb inds ee lnpfa mkt gov)
* 原始数据“内蒙”“黑龙”为截断简称，统一补齐
replace province = "内蒙古" if province=="内蒙"
replace province = "黑龙江" if province=="黑龙"
gen double lnDEI = ln(DEI)
gen double lnCEI = ln(CEI)
label var lnDEI "ln(DEI)"
label var lnCEI "ln(CEI)"

*---------------------------------------------------------------*
* 1. 省序编号 ord = 1..30（关键）
*    spatgsa 按数据集“行顺序”对应权重矩阵的行，因此必须先 sort ord。
*    这里不依赖 encode/排序的默认行为（Stata 对中文的显示与排序可能异常），
*    直接用显式映射写出省序；缺失即为省名不匹配，assert 会报错。
*---------------------------------------------------------------*
gen long ord = .
local i = 0
foreach p in 北京 天津 河北 山西 内蒙古 辽宁 吉林 黑龙江 上海 江苏 ///
             浙江 安徽 福建 江西 山东 河南 湖北 湖南 广东 广西 ///
             海南 重庆 四川 贵州 云南 陕西 甘肃 青海 宁夏 新疆 {
    local ++i
    replace ord = `i' if province=="`p'"
}
count if missing(ord)
assert r(N)==0
label var ord "省序（与权重矩阵行列一致）"
xtset ord year
compress
tempfile pan
save `pan'

*---------------------------------------------------------------*
* 2. spatwmat：导入三套矩阵并行标准化
*---------------------------------------------------------------*
spatwmat using FWadj.dta,  name(Wadj)  standardize
spatwmat using FWdist.dta, name(Wdist) standardize
spatwmat using FWecon.dta, name(Wecon) standardize

*---------------------------------------------------------------*
* 3. 逐年 Moran's I / Geary's c
*    spatgsa 返回 r(Moran) / r(Geary)，行=变量，列=[stat mean sd z p]
*    直接导出 r-class 结果，不手抄数
*---------------------------------------------------------------*
local varlist "DEI CEI lnDEI lnCEI"

capture file close fM
file open fM using "moran_long.txt", write replace text
file write fM "year W Wname variable I EI sd z p" _n

capture file close fG
file open fG using "geary_long.txt", write replace text
file write fG "year W Wname variable c Ec sd z p" _n

foreach Wname in Wadj Wdist Wecon {
    if "`Wname'"=="Wadj"  local wcn "GeographicContiguity"
    if "`Wname'"=="Wdist" local wcn "GeographicDistance"
    if "`Wname'"=="Wecon" local wcn "EconomicDistance"

    forvalues y = 2011/2022 {
        quietly {
            use `pan', clear
            keep if year==`y'
            sort ord
        }
        noisily spatgsa `varlist', weights(`Wname') moran geary twotail
        matrix MOR = r(Moran)
        matrix GEA = r(Geary)
        local i = 0
        foreach v of local varlist {
            local ++i
            file write fM "`y' `Wname' `wcn' `v' " ///
                %10.6f (MOR[`i',1]) %10.6f (MOR[`i',2]) %9.6f (MOR[`i',3]) ///
                %8.4f  (MOR[`i',4]) %9.6f  (MOR[`i',5]) _n
            file write fG "`y' `Wname' `wcn' `v' " ///
                %10.6f (GEA[`i',1]) %10.6f (GEA[`i',2]) %9.6f (GEA[`i',3]) ///
                %8.4f  (GEA[`i',4]) %9.6f  (GEA[`i',5]) _n
        }
        matrix drop MOR
        matrix drop GEA
    }
}
file close fM
file close fG

*---------------------------------------------------------------*
* 4. 论文宽表：DEI 与 CEI × 三个矩阵
*---------------------------------------------------------------*
capture file close fW
file open fW using "moran_width.txt", write replace text
file write fW "var" _char(9) "year" _char(9) "Wecon_I" _char(9) "Wecon_z" ///
    _char(9) "Wdist_I" _char(9) "Wdist_z" _char(9) "Wadj_I" _char(9) "Wadj_z" _n
foreach v in DEI CEI {
    forvalues y = 2011/2022 {
        quietly {
            use `pan', clear
            keep if year==`y'
            sort ord
            spatgsa `v', weights(Wecon) moran twotail
            matrix e = r(Moran)
            spatgsa `v', weights(Wdist) moran twotail
            matrix d = r(Moran)
            spatgsa `v', weights(Wadj) moran twotail
            matrix a = r(Moran)
        }
        * 固定宽度：保证任何行都是 8 个字段，便于 import delimited 读入
        file write fW "`v'" _char(9) "`y'" _char(9) ///
            %8.4f (e[1,1]) _char(9) %6.3f (e[1,4]) _char(9) ///
            %8.4f (d[1,1]) _char(9) %6.3f (d[1,4]) _char(9) ///
            %8.4f (a[1,1]) _char(9) %6.3f (a[1,4]) _n
        matrix drop e
        matrix drop d
        matrix drop a
    }
}
file close fW

type moran_width.txt
display as result "完成：moran_long.txt / geary_long.txt / moran_width.txt"