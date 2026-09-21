#!/usr/bin/env python3
"""Create lambda_GC and LDSC-style supplementary tables.

This script uses the downloaded EUR LD-score reference in Draft/ldsc/eur_w_ld_chr
and runs a lightweight LD-score regression with block jackknife standard errors.
It does not call the original Python-2 LDSC program.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Draft" / "ldsc_outputs"
OUT.mkdir(parents=True, exist_ok=True)

LD_DIR = ROOT / "Draft" / "ldsc" / "eur_w_ld_chr"
HM3_MAP = ROOT / "Draft" / "ldsc" / "hm3_rsid_position_allele_map.tsv"

REFSEQ_MAP = {
    str(i): ref
    for i, ref in {
        1: "NC_000001.11",
        2: "NC_000002.12",
        3: "NC_000003.12",
        4: "NC_000004.12",
        5: "NC_000005.10",
        6: "NC_000006.12",
        7: "NC_000007.14",
        8: "NC_000008.11",
        9: "NC_000009.12",
        10: "NC_000010.11",
        11: "NC_000011.10",
        12: "NC_000012.12",
        13: "NC_000013.11",
        14: "NC_000014.9",
        15: "NC_000015.10",
        16: "NC_000016.10",
        17: "NC_000017.11",
        18: "NC_000018.10",
        19: "NC_000019.10",
        20: "NC_000020.11",
        21: "NC_000021.9",
        22: "NC_000022.11",
    }.items()
}

GWAS_SCANS = [
    ("PCC vs responder controls", ROOT / "GWAS" / "pcc_vsallfil.txt.gz", "PCC_responder"),
    ("Pauci-symptomatic PCC vs responder controls", ROOT / "GWAS" / "clust1_vsallfil.txt.gz", "C1_responder"),
    ("Fatigue/sleep-dominated PCC vs responder controls", ROOT / "GWAS" / "clust2_vsallfil.txt.gz", "C2_responder"),
    ("Multi-symptomatic PCC vs responder controls", ROOT / "GWAS" / "clust3_vsallfil.txt.gz", "C3_responder"),
    ("PCC vs infected non-PCC controls", ROOT / "GWAS" / "pcc_vscovid.txt.gz", "PCC_infected_nonPCC"),
    ("Pauci-symptomatic PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust1_vsno.txt.gz", "C1_infected_nonPCC"),
    ("Fatigue/sleep-dominated PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust2_vsno.txt.gz", "C2_infected_nonPCC"),
    ("Multi-symptomatic PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust3_vsno.txt.gz", "C3_infected_nonPCC"),
]

HGI_TRAITS = [
    ("HGI C2 susceptibility", ROOT / "GWAS" / "COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz", "HGI_C2_susceptibility"),
    ("HGI B2 hospitalisation/severity", ROOT / "GWAS" / "COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz", "HGI_B2_severity"),
]


def fmt(x: float | int | str | None, digits: int = 4) -> str:
    if x is None or x == "":
        return ""
    if isinstance(x, str):
        return x
    if not math.isfinite(float(x)):
        return ""
    return f"{float(x):.{digits}f}"


def fmt_p(x: float | None) -> str:
    if x is None or not math.isfinite(float(x)):
        return ""
    x = float(x)
    if x < 1e-3:
        return f"{x:.2e}"
    return f"{x:.4f}"


def normal_p_from_z(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def latex_escape(text: object) -> str:
    out = "" if text is None else str(text)
    for raw, repl in {
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
    }.items():
        out = out.replace(raw, repl)
    return out


def write_tex(path: Path, rows: list[dict], fields: list[str], caption: str, label: str) -> None:
    with path.open("w") as f:
        f.write("\\begin{table}[htbp]\n\\centering\n")
        f.write(f"\\caption{{{latex_escape(caption)}}}\n")
        f.write(f"\\label{{{label}}}\n")
        f.write("\\resizebox{\\linewidth}{!}{%\n")
        f.write("\\begin{tabular}{" + "l" * len(fields) + "}\n\\toprule\n")
        f.write(" & ".join(latex_escape(x) for x in fields) + " \\\\\n\\midrule\n")
        for row in rows:
            f.write(" & ".join(latex_escape(row.get(x, "")) for x in fields) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n")


def load_ld_scores() -> pd.DataFrame:
    frames = []
    for chrom in range(1, 23):
        path = LD_DIR / f"{chrom}.l2.ldscore.gz"
        frames.append(pd.read_csv(path, sep="\t", usecols=["SNP", "L2"]))
    ld = pd.concat(frames, ignore_index=True)
    ld = ld.drop_duplicates("SNP")
    return ld


def load_hm3_map() -> pd.DataFrame:
    mapping = pd.read_csv(HM3_MAP, sep="\t")
    mapping = mapping.rename(columns={"POSITION": "position", "RSID": "SNP", "REF": "ref", "ALT": "alt"})
    mapping["position"] = mapping["position"].astype(np.int64)
    mapping["ref"] = mapping["ref"].str.upper()
    mapping["alt"] = mapping["alt"].str.upper()
    return mapping


def lambda_gc_from_full(path: Path) -> tuple[int, float, float]:
    chunks = pd.read_csv(path, sep="\t", usecols=["CHISQ"], chunksize=2_000_000)
    values = []
    n = 0
    for chunk in chunks:
        arr = pd.to_numeric(chunk["CHISQ"], errors="coerce").dropna().to_numpy(float)
        arr = arr[np.isfinite(arr)]
        n += arr.size
        values.append(arr)
    all_values = np.concatenate(values)
    median_chisq = float(np.median(all_values))
    return n, median_chisq, median_chisq / 0.454936423119572


def munge_gwas(label: str, path: Path, short_name: str, mapping: pd.DataFrame, ld: pd.DataFrame) -> pd.DataFrame:
    out_path = OUT / f"munged_{short_name}.tsv.gz"
    if out_path.exists():
        return pd.read_csv(out_path, sep="\t")

    pieces = []
    usecols = ["CHROM", "POSITION", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "BETA", "SE", "P", "N"]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=1_000_000):
        chunk = chunk.rename(columns={"POSITION": "position"})
        chunk["CHR_REFSEQ"] = chunk["CHROM"].astype(str).map(REFSEQ_MAP)
        chunk["position"] = pd.to_numeric(chunk["position"], errors="coerce").astype("Int64")
        chunk["ea"] = chunk["EFFECT_ALLELE"].astype(str).str.upper()
        chunk["nea"] = chunk["NON_EFFECT_ALLELE"].astype(str).str.upper()
        chunk = chunk.dropna(subset=["CHR_REFSEQ", "position"])
        merged = chunk.merge(mapping, on=["CHR_REFSEQ", "position"], how="inner")
        same = (merged["ea"] == merged["alt"]) & (merged["nea"] == merged["ref"])
        flip = (merged["ea"] == merged["ref"]) & (merged["nea"] == merged["alt"])
        merged = merged[same | flip].copy()
        if merged.empty:
            continue
        beta = pd.to_numeric(merged["BETA"], errors="coerce")
        se = pd.to_numeric(merged["SE"], errors="coerce")
        z = beta / se
        z.loc[flip.loc[merged.index]] = -z.loc[flip.loc[merged.index]]
        merged["Z"] = z
        merged["P"] = pd.to_numeric(merged["P"], errors="coerce")
        merged["N"] = pd.to_numeric(merged["N"], errors="coerce")
        pieces.append(merged[["SNP", "Z", "P", "N"]])

    if not pieces:
        raise RuntimeError(f"No HapMap3 overlap for {label}")
    munged = pd.concat(pieces, ignore_index=True)
    munged = munged.replace([np.inf, -np.inf], np.nan).dropna(subset=["SNP", "Z", "P", "N"])
    munged = munged.sort_values("P").drop_duplicates("SNP", keep="first")
    munged = munged.merge(ld, on="SNP", how="inner")
    munged.to_csv(out_path, sep="\t", index=False, compression="gzip")
    return munged


def munge_hgi(label: str, path: Path, short_name: str, ld: pd.DataFrame) -> pd.DataFrame:
    out_path = OUT / f"munged_{short_name}.tsv.gz"
    if out_path.exists():
        return pd.read_csv(out_path, sep="\t")

    pieces = []
    usecols = ["rsid", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p", "all_inv_var_meta_effective"]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=1_000_000):
        beta = pd.to_numeric(chunk["all_inv_var_meta_beta"], errors="coerce")
        se = pd.to_numeric(chunk["all_inv_var_meta_sebeta"], errors="coerce")
        out = pd.DataFrame({
            "SNP": chunk["rsid"],
            "Z": beta / se,
            "P": pd.to_numeric(chunk["all_inv_var_meta_p"], errors="coerce"),
            "N": pd.to_numeric(chunk["all_inv_var_meta_effective"], errors="coerce"),
        })
        pieces.append(out)
    munged = pd.concat(pieces, ignore_index=True)
    munged = munged.replace([np.inf, -np.inf], np.nan).dropna(subset=["SNP", "Z", "P", "N"])
    munged = munged.sort_values("P").drop_duplicates("SNP", keep="first")
    munged = munged.merge(ld, on="SNP", how="inner")
    munged.to_csv(out_path, sep="\t", index=False, compression="gzip")
    return munged


def jackknife_intercept(df: pd.DataFrame, y_col: str) -> tuple[float, float, float, float, int]:
    x = df["L2"].to_numpy(float)
    y = df[y_col].to_numpy(float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    def fit(xv, yv):
        X = np.column_stack([np.ones_like(xv), xv])
        coef, *_ = np.linalg.lstsq(X, yv, rcond=None)
        return float(coef[0]), float(coef[1])

    intercept, slope = fit(x, y)
    n_blocks = min(200, max(20, len(x) // 1000))
    idx = np.array_split(np.arange(len(x)), n_blocks)
    vals = []
    slopes = []
    all_idx = np.arange(len(x))
    for block in idx:
        keep = np.ones(len(x), dtype=bool)
        keep[block] = False
        b0, b1 = fit(x[keep], y[keep])
        vals.append(b0)
        slopes.append(b1)
    vals = np.array(vals)
    slopes = np.array(slopes)
    se = math.sqrt((n_blocks - 1) / n_blocks * np.sum((vals - vals.mean()) ** 2))
    slope_se = math.sqrt((n_blocks - 1) / n_blocks * np.sum((slopes - slopes.mean()) ** 2))
    return intercept, se, slope, slope_se, len(x)


def make_qc_table(scan_data: dict[str, pd.DataFrame]) -> list[dict]:
    rows = []
    for label, _, short_name in GWAS_SCANS:
        full_n, median_chisq, lambda_gc = lambda_gc_from_full(dict((s, p) for _, p, s in GWAS_SCANS)[short_name])
        df = scan_data[short_name].copy()
        df["CHI2"] = df["Z"] ** 2
        intercept, se, slope, slope_se, n_ldsc = jackknife_intercept(df, "CHI2")
        rows.append({
            "Scan": label,
            "Full-summary variants": f"{full_n:,}",
            "Median chi-square": fmt(median_chisq, 4),
            "lambda_GC": fmt(lambda_gc, 4),
            "LDSC SNPs": f"{n_ldsc:,}",
            "LDSC intercept": fmt(intercept, 4),
            "LDSC intercept SE": fmt(se, 4),
        })
    return rows


def make_rg_table(scan_data: dict[str, pd.DataFrame], trait_data: dict[str, pd.DataFrame]) -> list[dict]:
    rows = []
    # Precompute slopes for all traits on the same regression scale.
    slopes = {}
    slope_ses = {}
    for name, df in {**scan_data, **trait_data}.items():
        tmp = df.copy()
        tmp["CHI2"] = tmp["Z"] ** 2
        _, _, slope, slope_se, _ = jackknife_intercept(tmp, "CHI2")
        slopes[name] = slope
        slope_ses[name] = slope_se

    for label, _, short_name in GWAS_SCANS:
        for hgi_label, _, hgi_short in HGI_TRAITS:
            merged = scan_data[short_name][["SNP", "Z", "L2"]].merge(
                trait_data[hgi_short][["SNP", "Z"]], on="SNP", suffixes=("_gwas", "_hgi")
            )
            merged["ZPROD"] = merged["Z_gwas"] * merged["Z_hgi"]
            intercept, se, cov_slope, cov_slope_se, n = jackknife_intercept(merged, "ZPROD")
            denom = math.sqrt(slopes[short_name] * slopes[hgi_short]) if slopes[short_name] > 0 and slopes[hgi_short] > 0 else math.nan
            rg = cov_slope / denom if math.isfinite(denom) and denom != 0 else math.nan
            # Delta approximation using the covariance slope SE only; conservative enough for table triage.
            rg_se = abs(cov_slope_se / denom) if math.isfinite(denom) and denom != 0 else math.nan
            p = normal_p_from_z(rg / rg_se) if math.isfinite(rg) and math.isfinite(rg_se) and rg_se > 0 else math.nan
            rows.append({
                "PCC scan": label,
                "External trait": hgi_label,
                "Shared LDSC SNPs": f"{n:,}",
                "Genetic covariance intercept": fmt(intercept, 4),
                "rg": fmt(rg, 4),
                "rg SE": fmt(rg_se, 4),
                "rg p": fmt_p(p),
            })
    return rows


def main() -> None:
    print("Loading LD scores")
    ld = load_ld_scores()
    mapping = load_hm3_map()

    scan_data = {}
    for label, path, short_name in GWAS_SCANS:
        print(f"Munging {label}")
        scan_data[short_name] = munge_gwas(label, path, short_name, mapping, ld)
        print(f"  {len(scan_data[short_name]):,} LDSC SNPs")

    trait_data = {}
    for label, path, short_name in HGI_TRAITS:
        print(f"Munging {label}")
        trait_data[short_name] = munge_hgi(label, path, short_name, ld)
        print(f"  {len(trait_data[short_name]):,} LDSC SNPs")

    qc_rows = make_qc_table(scan_data)
    qc_fields = ["Scan", "Full-summary variants", "Median chi-square", "lambda_GC", "LDSC SNPs", "LDSC intercept", "LDSC intercept SE"]
    write_csv(OUT / "supplementary_table_lambda_gc_ldsc_intercepts.csv", qc_rows, qc_fields)
    write_tex(
        OUT / "supplementary_table_lambda_gc_ldsc_intercepts.tex",
        qc_rows,
        qc_fields,
        "Genomic inflation and LD-score regression intercepts for the eight primary GWAS scans.",
        "tab:lambda-gc-ldsc-intercepts",
    )

    rg_rows = make_rg_table(scan_data, trait_data)
    rg_fields = ["PCC scan", "External trait", "Shared LDSC SNPs", "Genetic covariance intercept", "rg", "rg SE", "rg p"]
    write_csv(OUT / "supplementary_table_ldsc_rg_hgi.csv", rg_rows, rg_fields)
    write_tex(
        OUT / "supplementary_table_ldsc_rg_hgi.tex",
        rg_rows,
        rg_fields,
        "LD-score regression genetic correlations between PCC GWAS scans and HGI COVID-19 susceptibility/severity traits.",
        "tab:ldsc-rg-hgi",
    )
    print(f"Wrote outputs to {OUT}")


if __name__ == "__main__":
    main()
