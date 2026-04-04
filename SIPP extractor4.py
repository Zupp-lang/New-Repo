# 01_extract.py

#!/usr/bin/env python3
"""Stage 1: Extract selected columns from yearly SIPP .dta files to CSV."""
 
from pathlib import Path
import os
import pandas as pd
import pyreadstat
 
 
# ── Config loader ───────────────────────────────────────────────────
 
def read_env_config(path: Path) -> dict:
    cfg = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or "=" not in s:
            continue
        k, v = s.split("=", 1)
        cfg[k.strip()] = v.strip()
    return cfg
 
 
SCRIPT_DIR = Path(__file__).resolve().parent
CFG_PATH = SCRIPT_DIR / "pipeline_config.env"
if not CFG_PATH.exists():
    raise FileNotFoundError(f"Missing shared config: {CFG_PATH}")
CFG = read_env_config(CFG_PATH)
 
base_folder = CFG["BASE_DIR"]
years = [int(x) for x in CFG["YEARS"].split(",")]
clean_tmpl = CFG["CLEAN_EXTRACT_TEMPLATE"]
 
print("[config] BASE_DIR:", base_folder)
print("[config] YEARS:", years)
print("[config] CLEAN_EXTRACT_TEMPLATE:", clean_tmpl)
if not os.path.isdir(base_folder):
    raise FileNotFoundError(f"Configured BASE_DIR does not exist: {base_folder}")
 
 
# ── Variables to extract ────────────────────────────────────────────
# Identifiers and demographics
vars_core = [
    "SSUID", "SHHNUM", "SPUID", "TEHC_ST", "PNUM", "RHCALM", "RHCALYR",
    "TAGE", "RRELATE", "ERELATE", "ERELRPE", "EEDUC", "TMWKHRS", "TPEARN",
    "SHHADID", "MONTHCODE", "TNETWORTH", "EORIGIN", "TRACE", "SPANEL",
]
 
# TANF variables — both receipt indicator AND exit reason / lifetime limit
vars_tanf = [
    "RTANF_MNYN",       # <-- KEY: ever-received-TANF flag (month-level)
    "TTANF_AMT",        # TANF amount
    "TTANF_ERSN",       # TANF exit reason (3 = lifetime limit)
    "RTANF_LCYR",       # lifetime limit year indicator
]
 
# Employment / disability / job search
vars_employment = [
    "EJOBCANT", "ENJ_NOWRK2", "ENJ_LKWRK",
    "EERRANDS", "ESELFCARE", "EAMBULAT", "ECOGNIT", "ESEEING", "EHEARING",
]
 
# Relationship variables (for parent identification)
vars_rel = [f"RREL{i}" for i in range(1, 31)] + [f"RREL_PNUM{i}" for i in range(1, 15)]
 
# Childcare variables
vars_childcare = [
    "TPAYWK", "EPAY", "EPAYHELP",
    "EWHOPAID1", "EWHOPAID2", "EWHOPAID3", "EWHOPAID4", "EWHOPAID5",
    "EDAYCARE", "EDAYHS", "EFAM", "EGRAN", "EHEADST", "EKIDSELF",
    "ENREL", "ENUR", "ENURHS", "EOTHR", "EPAR", "EPROG", "ESELF",
    "ESIB15", "ELIST", "EWORKMORE", "ETIMELOST", "ETIMELOST_TP",
]
 
vars_to_keep = list(dict.fromkeys(
    vars_core + vars_tanf + vars_employment + vars_rel + vars_childcare
))
 
# Key diagnostic variables — we'll report their availability per year
DIAG_VARS = ["RTANF_MNYN", "TTANF_ERSN", "RTANF_LCYR", "EPAY", "TPAYWK",
             "TEHC_ST", "TAGE", "TPEARN"]
 
 
# ── Helpers ─────────────────────────────────────────────────────────
 
def coerce_ssuid(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    out = series.astype("string")
    ok = numeric.notna()
    if ok.any():
        out.loc[ok] = numeric[ok].round().astype("Int64").astype("string")
    return out
 
 
# ── Year loop ───────────────────────────────────────────────────────
 
for yr in years:
    input_path = os.path.join(base_folder, f"pu{yr}_dta", f"pu{yr}.dta")
    output_path = os.path.join(base_folder, clean_tmpl.format(year=yr))
 
    print(f"\n{'='*60}")
    print(f"  YEAR {yr}")
    print(f"{'='*60}")
    print("  input :", input_path)
    print("  output:", output_path)
 
    if not os.path.exists(input_path):
        print(f"  !! MISSING file, skipping: {input_path}")
        continue
 
    # Try reading with all requested columns; fall back to available subset
    try:
        df, _ = pyreadstat.read_dta(input_path, usecols=vars_to_keep)
    except Exception:
        df_probe, _ = pyreadstat.read_dta(input_path, row_limit=1)
        available = [v for v in vars_to_keep if v in df_probe.columns]
        missing = sorted(set(vars_to_keep) - set(available))
        if missing:
            print(f"  note: {len(missing)} requested vars absent in {yr}")
            if len(missing) <= 15:
                print(f"         {missing}")
            else:
                print(f"         first 15: {missing[:15]}")
        df, _ = pyreadstat.read_dta(input_path, usecols=available)
 
    # Type coercion
    if "TEHC_ST" in df.columns:
        df["TEHC_ST"] = pd.to_numeric(df["TEHC_ST"], errors="coerce")
    if "SSUID" in df.columns:
        df["SSUID"] = coerce_ssuid(df["SSUID"])
 
    # Per-year diagnostics for key variables
    print(f"\n  --- Diagnostics for {yr} ---")
    print(f"  Total rows: {len(df):,}")
    for dv in DIAG_VARS:
        if dv in df.columns:
            non_na = df[dv].notna().sum()
            pct = 100 * non_na / len(df) if len(df) > 0 else 0
            if dv == "RTANF_MNYN":
                n_tanf = (pd.to_numeric(df[dv], errors="coerce") == 1).sum()
                print(f"  {dv:20s}: {non_na:>8,} non-NA ({pct:5.1f}%), "
                      f"==1 (TANF receipt): {n_tanf:,}")
            elif dv == "EPAY":
                n_pay = (pd.to_numeric(df[dv], errors="coerce") == 1).sum()
                print(f"  {dv:20s}: {non_na:>8,} non-NA ({pct:5.1f}%), "
                      f"==1 (pays for care): {n_pay:,}")
            elif dv == "TPAYWK":
                vals = pd.to_numeric(df[dv], errors="coerce")
                n_pos = (vals > 0).sum()
                print(f"  {dv:20s}: {non_na:>8,} non-NA ({pct:5.1f}%), "
                      f">0: {n_pos:,}")
            else:
                print(f"  {dv:20s}: {non_na:>8,} non-NA ({pct:5.1f}%)")
        else:
            print(f"  {dv:20s}: ** ABSENT in {yr} **")
 
    df.to_csv(output_path, index=False)
    print(f"\n  Saved: {output_path}  shape: {df.shape}")
 
print(f"\n{'='*60}")
print("Done stage 1.")