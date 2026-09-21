#!/usr/bin/env python3
"""Run LDSC after assigning LC rsIDs from HGI chr:pos:allele matches.

The LC WGS summary statistics do not contain rsIDs. For LDSC, the least lossy
local harmonisation is to use the HGI C2 summary statistics as the rsID/build
reference, then align LC and HGI B2 alleles to that same HGI C2 REF/ALT pair.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Draft" / "ldsc_hgi_position_outputs"
OUT.mkdir(parents=True, exist_ok=True)

LD_PREFIX = str(ROOT / "Draft" / "ldsc" / "eur_w_ld_chr") + "/"
HGI_C2 = ROOT / "GWAS" / "COVID19_HGI_C2_ALL_eur_leave23andme_20220403.tsv.gz"
HGI_B2 = ROOT / "GWAS" / "COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv.gz"

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


def fmt(x, digits=4):
    if x is None or x == "":
        return ""
    try:
        x = float(x)
    except ValueError:
        return str(x)
    if not math.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_p(x):
    if x is None:
        return ""
    x = float(x)
    if not math.isfinite(x):
        return ""
    return f"{x:.2e}" if x < 1e-3 else f"{x:.4f}"


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def latex_escape(text: object) -> str:
    out = "" if text is None else str(text)
    for raw, repl in {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_"}.items():
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


def load_ld_snps() -> set[str]:
    snps = set()
    for chrom in range(1, 23):
        with gzip.open(Path(LD_PREFIX) / f"{chrom}.l2.ldscore.gz", "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            snp_i = header.index("SNP")
            for line in f:
                snps.add(line.rstrip("\n").split("\t")[snp_i])
    return snps


def load_hgi_c2_reference(ld_snps: set[str]) -> pd.DataFrame:
    pieces = []
    for chunk in pd.read_csv(HGI_C2, sep="\t", usecols=["#CHR", "POS", "REF", "ALT", "rsid"], chunksize=1_000_000):
        chunk = chunk.rename(columns={"#CHR": "chr", "POS": "pos"})
        chunk["REF"] = chunk["REF"].astype(str).str.upper()
        chunk["ALT"] = chunk["ALT"].astype(str).str.upper()
        chunk = chunk[
            chunk["rsid"].isin(ld_snps)
            & chunk["REF"].isin(["A", "C", "G", "T"])
            & chunk["ALT"].isin(["A", "C", "G", "T"])
        ]
        pieces.append(chunk)
    ref = pd.concat(pieces, ignore_index=True)
    ref = ref.drop_duplicates(["chr", "pos", "REF", "ALT", "rsid"])
    ref = ref.sort_values("rsid").drop_duplicates(["chr", "pos"], keep="first")
    return ref


def count_rows(path: Path) -> int:
    with gzip.open(path, "rt") as f:
        return max(sum(1 for _ in f) - 1, 0)


def write_sumstats(out_path: Path, chunks) -> int:
    n = 0
    with gzip.open(out_path, "wt") as f:
        f.write("SNP\tA1\tA2\tZ\tP\tN\n")
        for frame in chunks:
            if frame.empty:
                continue
            frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["SNP", "A1", "A2", "Z", "P", "N"])
            frame = frame.sort_values("P").drop_duplicates("SNP", keep="first")
            frame[["SNP", "A1", "A2", "Z", "P", "N"]].to_csv(f, sep="\t", header=False, index=False)
            n += len(frame)
    return n


def munge_lc(label: str, path: Path, short_name: str, hgi_ref: pd.DataFrame) -> tuple[Path, int]:
    out_path = OUT / f"{short_name}.hgi_position.sumstats.gz"
    if out_path.exists():
        return out_path, count_rows(out_path)

    ref = hgi_ref.rename(columns={"chr": "CHROM", "pos": "POSITION"})

    def chunks():
        usecols = ["CHROM", "POSITION", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "BETA", "SE", "P", "N"]
        for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=1_000_000):
            chunk["EFFECT_ALLELE"] = chunk["EFFECT_ALLELE"].astype(str).str.upper()
            chunk["NON_EFFECT_ALLELE"] = chunk["NON_EFFECT_ALLELE"].astype(str).str.upper()
            merged = chunk.merge(ref, on=["CHROM", "POSITION"], how="inner")
            same = (merged["EFFECT_ALLELE"] == merged["ALT"]) & (merged["NON_EFFECT_ALLELE"] == merged["REF"])
            flip = (merged["EFFECT_ALLELE"] == merged["REF"]) & (merged["NON_EFFECT_ALLELE"] == merged["ALT"])
            merged = merged[same | flip].copy()
            beta = pd.to_numeric(merged["BETA"], errors="coerce")
            se = pd.to_numeric(merged["SE"], errors="coerce")
            z = beta / se
            z.loc[flip.loc[merged.index]] = -z.loc[flip.loc[merged.index]]
            yield pd.DataFrame({
                "SNP": merged["rsid"],
                "A1": merged["ALT"],
                "A2": merged["REF"],
                "Z": z,
                "P": pd.to_numeric(merged["P"], errors="coerce"),
                "N": pd.to_numeric(merged["N"], errors="coerce"),
            })

    print(f"Munging {label}")
    return out_path, write_sumstats(out_path, chunks())


def munge_hgi(path: Path, label: str, short_name: str, hgi_ref: pd.DataFrame) -> tuple[Path, int]:
    out_path = OUT / f"{short_name}.hgi_position.sumstats.gz"
    if out_path.exists():
        return out_path, count_rows(out_path)

    ref = hgi_ref[["rsid", "REF", "ALT"]].drop_duplicates("rsid")

    def chunks():
        usecols = ["rsid", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p", "all_inv_var_meta_effective"]
        for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=1_000_000):
            chunk["REF"] = chunk["REF"].astype(str).str.upper()
            chunk["ALT"] = chunk["ALT"].astype(str).str.upper()
            merged = chunk.merge(ref, on="rsid", how="inner", suffixes=("", "_ref"))
            same = (merged["REF"] == merged["REF_ref"]) & (merged["ALT"] == merged["ALT_ref"])
            flip = (merged["REF"] == merged["ALT_ref"]) & (merged["ALT"] == merged["REF_ref"])
            merged = merged[same | flip].copy()
            beta = pd.to_numeric(merged["all_inv_var_meta_beta"], errors="coerce")
            se = pd.to_numeric(merged["all_inv_var_meta_sebeta"], errors="coerce")
            z = beta / se
            z.loc[flip.loc[merged.index]] = -z.loc[flip.loc[merged.index]]
            yield pd.DataFrame({
                "SNP": merged["rsid"],
                "A1": merged["ALT_ref"],
                "A2": merged["REF_ref"],
                "Z": z,
                "P": pd.to_numeric(merged["all_inv_var_meta_p"], errors="coerce"),
                "N": pd.to_numeric(merged["all_inv_var_meta_effective"], errors="coerce"),
            })

    print(f"Munging {label}")
    return out_path, write_sumstats(out_path, chunks())


def lambda_gc_from_full(path: Path) -> tuple[int, float, float]:
    values = []
    n = 0
    for chunk in pd.read_csv(path, sep="\t", usecols=["CHISQ"], chunksize=2_000_000):
        arr = pd.to_numeric(chunk["CHISQ"], errors="coerce").dropna().to_numpy(float)
        arr = arr[np.isfinite(arr)]
        n += arr.size
        values.append(arr)
    all_values = np.concatenate(values)
    median_chisq = float(np.median(all_values))
    return n, median_chisq, median_chisq / 0.454936423119572


def run_ldsc(args: list[str], out_prefix: Path) -> str:
    cmd = [
        "uv", "tool", "run", "--python", "3.11",
        "--with", "pandas<2", "--with", "numpy<2",
        "--from", "ldsc", "ldsc.py",
        *args, "--out", str(out_prefix),
    ]
    subprocess.run(cmd, cwd=ROOT, check=False)
    return out_prefix.with_suffix(".log").read_text()


def parse_h2_log(text: str) -> dict:
    n = re.search(r"After merging with regression SNP LD, ([0-9]+) SNPs remain", text)
    intercept = re.search(r"Intercept: ([^ ]+) \(([^)]+)\)", text)
    return {
        "n": int(n.group(1)) if n else None,
        "intercept": float(intercept.group(1)) if intercept else None,
        "intercept_se": float(intercept.group(2)) if intercept else None,
    }


def parse_optional_float(value: str):
    if value in {"NA", "nan", "NaN", ""}:
        return None
    try:
        value = float(value)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def parse_rg_log(text: str) -> list[dict]:
    marker = "Summary of Genetic Correlation Results"
    if marker not in text:
        return []
    rows = []
    lines = text.split(marker, 1)[1].strip().splitlines()
    if len(lines) < 2:
        return rows
    header = re.split(r"\s+", lines[0].strip())
    for line in lines[1:]:
        if not line.strip() or line.startswith("Analysis finished"):
            break
        parts = re.split(r"\s+", line.strip())
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))
    return rows


def main() -> None:
    ld_snps = load_ld_snps()
    hgi_ref = load_hgi_c2_reference(ld_snps)
    print(f"HGI C2 LD-score reference SNPs: {len(hgi_ref):,}")

    sumstats = {}
    for label, path, short in GWAS_SCANS:
        path_out, n = munge_lc(label, path, short, hgi_ref)
        sumstats[short] = path_out
        print(f"  {short}: {n:,}")

    for label, path, short in [
        ("HGI C2 susceptibility", HGI_C2, "HGI_C2_susceptibility"),
        ("HGI B2 hospitalisation/severity", HGI_B2, "HGI_B2_severity"),
    ]:
        path_out, n = munge_hgi(path, label, short, hgi_ref)
        sumstats[short] = path_out
        print(f"  {short}: {n:,}")

    qc_rows = []
    for label, path, short in GWAS_SCANS:
        log = run_ldsc([
            "--h2", str(sumstats[short]),
            "--ref-ld-chr", LD_PREFIX,
            "--w-ld-chr", LD_PREFIX,
        ], OUT / f"ldsc_h2_{short}")
        parsed = parse_h2_log(log)
        full_n, med, lam = lambda_gc_from_full(path)
        qc_rows.append({
            "Scan": label,
            "Full-summary variants": f"{full_n:,}",
            "Median chi-square": fmt(med),
            "lambda_GC": fmt(lam),
            "LDSC SNPs": f"{parsed['n']:,}" if parsed["n"] else "",
            "LDSC intercept": fmt(parsed["intercept"]),
            "LDSC intercept SE": fmt(parsed["intercept_se"]),
        })

    qc_fields = ["Scan", "Full-summary variants", "Median chi-square", "lambda_GC", "LDSC SNPs", "LDSC intercept", "LDSC intercept SE"]
    write_csv(OUT / "supplementary_table_lambda_gc_ldsc_intercepts.csv", qc_rows, qc_fields)
    write_tex(OUT / "supplementary_table_lambda_gc_ldsc_intercepts.tex", qc_rows, qc_fields, "Genomic inflation and LD-score regression intercepts for the eight primary GWAS scans after HGI-position rsID harmonisation.", "tab:lambda-gc-ldsc-intercepts")

    rg_rows = []
    for label, _, short in GWAS_SCANS:
        log = run_ldsc([
            "--rg", ",".join([str(sumstats[short]), str(sumstats["HGI_C2_susceptibility"]), str(sumstats["HGI_B2_severity"])]),
            "--ref-ld-chr", LD_PREFIX,
            "--w-ld-chr", LD_PREFIX,
        ], OUT / f"ldsc_rg_{short}_vs_HGI")
        rows = parse_rg_log(log)
        if not rows:
            for hgi in ["HGI C2 susceptibility", "HGI B2 hospitalisation/severity"]:
                rg_rows.append({"PCC scan": label, "External trait": hgi, "rg": "", "rg SE": "", "z": "", "rg p": "", "Status": "Unavailable"})
            continue
        for row in rows:
            p2 = Path(row["p2"]).name.replace(".hgi_position.sumstats.gz", "")
            hgi = {
                "HGI_C2_susceptibility": "HGI C2 susceptibility",
                "HGI_B2_severity": "HGI B2 hospitalisation/severity",
            }.get(p2, p2)
            rg = parse_optional_float(row["rg"])
            se = parse_optional_float(row["se"])
            z = parse_optional_float(row["z"])
            p = parse_optional_float(row["p"])
            rg_rows.append({
                "PCC scan": label,
                "External trait": hgi,
                "rg": fmt(rg),
                "rg SE": fmt(se),
                "z": fmt(z, 3),
                "rg p": fmt_p(p),
                "Status": "OK" if rg is not None and abs(rg) <= 1 else "Unstable/unavailable",
            })
    rg_fields = ["PCC scan", "External trait", "rg", "rg SE", "z", "rg p", "Status"]
    write_csv(OUT / "supplementary_table_ldsc_rg_hgi.csv", rg_rows, rg_fields)
    write_tex(OUT / "supplementary_table_ldsc_rg_hgi.tex", rg_rows, rg_fields, "LD-score regression genetic correlations between PCC GWAS scans and HGI COVID-19 susceptibility/severity traits after HGI-position rsID harmonisation.", "tab:ldsc-rg-hgi")
    print(f"Wrote outputs to {OUT}")


if __name__ == "__main__":
    main()
