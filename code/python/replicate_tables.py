"""Replicate the baseline and robustness estimates reported in the paper.

Requires: pandas, numpy, linearmodels.
Run from anywhere: python code/python/replicate_tables.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from linearmodels.panel import PanelOLS

DATA = Path(__file__).resolve().parents[2] / "data" / "panel_AI5.csv"

df = pd.read_csv(DATA, encoding="utf-8-sig")
df["lnDEI"] = np.log(df["DEI"])
df["lnCEI"] = np.log(df["CEI"])
df = df.set_index(["sh", "年份"])

controls = ["lnpgdp", "fisc", "dens", "fin", "urb", "mkt"]


def twoway_demean(s: pd.Series) -> pd.Series:
    """Remove province and year means (two-way within transformation)."""
    return (s
            - s.groupby(level=0).transform("mean")
            - s.groupby(level=1).transform("mean")
            + s.mean())


def within_r2(d: pd.DataFrame, y: str, resids: pd.Series) -> float:
    """Within R-squared on the two-way demeaned outcome, as in Table 4."""
    yt = twoway_demean(d[y])
    return float(1 - (resids ** 2).sum() / (yt ** 2).sum())


def fit(y, xs, data):
    d = data.dropna(subset=[y] + xs)
    mod = PanelOLS(d[y], d[xs], entity_effects=True, time_effects=True,
                   drop_absorbed=True)
    res = mod.fit(cov_type="clustered", cluster_entity=True, cluster_time=False)
    res._r2_within = within_r2(d, y, res.resids)
    return res


def report(label, res, xs):
    print(f"\n== {label}  N = {res.nobs}  within R2 = {res._r2_within:.3f}")
    for v in xs:
        print(f"{v:10s} {res.params[v]: .4f}  se = {res.std_errors[v]:.4f}  "
              f"p = {res.pvalues[v]:.3f}")


# --- baseline and the two channel specifications -------------------------
for y, xs, label in [
    ("lnCEI", ["lnDEI"] + controls, "baseline"),
    ("lnCE", ["lnDEI"] + controls, "outcome: ln carbon intensity"),
    ("lnCEI", ["lnDEI", "inds"] + controls, "lnDEI + inds"),
    ("lnCEI", ["lnDEI"] + controls + ["inds"], "inds added to controls"),
]:
    report(label, fit(y, xs, df), xs)

# --- lagged DEI ----------------------------------------------------------
df2 = df.copy()
df2["L_lnDEI"] = df2.groupby(level=0)["lnDEI"].shift(1)
r = fit("lnCEI", ["L_lnDEI"] + controls, df2)
print(f"\n== lagged DEI  N = {r.nobs}  within R2 = {r._r2_within:.3f}")
print(f"L_lnDEI    {r.params['L_lnDEI']:.4f}  se = {r.std_errors['L_lnDEI']:.4f}  "
      f"p = {r.pvalues['L_lnDEI']:.3f}")

# --- exclude the four municipalities ------------------------------------
mun = ["北京", "上海", "天津", "重庆"]
d3 = df[~df.index.get_level_values(0).isin(mun)]
r = fit("lnCEI", ["lnDEI"] + controls, d3)
print(f"\n== exclude municipalities  N = {r.nobs}  within R2 = {r._r2_within:.3f}")
print(f"lnDEI      {r.params['lnDEI']:.4f}  se = {r.std_errors['lnDEI']:.4f}  "
      f"p = {r.pvalues['lnDEI']:.3f}")

# --- 1% winsorization ----------------------------------------------------
d4 = df.copy()
for c in ["DEI", "CEI"]:
    lo, hi = d4[c].quantile([0.01, 0.99])
    d4[c] = d4[c].clip(lo, hi)
d4["lnDEI"] = np.log(d4["DEI"])
d4["lnCEI"] = np.log(d4["CEI"])
r = fit("lnCEI", ["lnDEI"] + controls, d4)
print(f"\n== 1% winsorization  N = {r.nobs}  within R2 = {r._r2_within:.3f}")
print(f"lnDEI      {r.params['lnDEI']:.4f}  se = {r.std_errors['lnDEI']:.4f}  "
      f"p = {r.pvalues['lnDEI']:.3f}")

# --- descriptive statistics ---------------------------------------------
print("\n== descriptive statistics")
print(df[["DEI", "CEI", "lnDEI", "lnCEI", "lnCE"]]
      .describe().T[["mean", "std", "min", "max"]].round(4))
print(f"\nprovinces = {df.index.get_level_values(0).nunique()}, "
      f"corr(DEI, CEI) = {df['DEI'].corr(df['CEI']):.3f}")
