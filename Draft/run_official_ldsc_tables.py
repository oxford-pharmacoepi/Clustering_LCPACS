#!/usr/bin/env python3
"""Munge PCC GWAS and HGI summary statistics, run LDSC, and write tables."""

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
OUT = ROOT / "Draft" / "ldsc_official_outputs"
OUT.mkdir(parents=True, exist_ok=True)

LD_PREFIX = str(ROOT / "Draft" / "ldsc" / "eur_w_ld_chr") + "/"
HM3_MAP = ROOT / "Draft" / "ldsc" / "hm3_rsid_position_allele_map_full.tsv"

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
    return f"{x:.2e}" if x < 1e-3 else f"{x:.4f}"


def parse_optional_float(value: str) -> float | None:
    if value in {"NA", "nan", "NaN", ""}:
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    return out if math.isfinite(out) else None


def latex_escape(text: object) -> str:
    out = "" if text is None else str(text)
    for raw, repl in {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_"}.items():
        out = out.replace(raw, repl)
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


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


def load_map() -> pd.DataFrame:
    mapping = pd.read_csv(HM3_MAP, sep="\t")
    mapping = mapping.rename(columns={"POSITION": "position", "RSID": "SNP", "REF": "ref", "ALT": "alt"})
    mapping["CHR_REFSEQ"] = mapping["CHR_REFSEQ"].astype(str)
    mapping["position"] = pd.to_numeric(mapping["position"], errors="coerce").astype("Int64")
    mapping["ref"] = mapping["ref"].astype(str).str.upper()
    mapping["alt"] = mapping["alt"].astype(str).str.upper()
    mapping["allele_pair"] = mapping.apply(lambda row: "/".join(sorted([row["ref"], row["alt"]])), axis=1)
    pair_counts = mapping[["SNP", "allele_pair"]].drop_duplicates().groupby("SNP").size()
    unique_pair_snps = set(pair_counts[pair_counts == 1].index)
    mapping = mapping[mapping["SNP"].isin(unique_pair_snps)].drop(columns=["allele_pair"])
    return mapping


def load_ld_snps() -> set[str]:
    snps: set[str] = set()
    for chrom in range(1, 23):
        path = Path(LD_PREFIX) / f"{chrom}.l2.ldscore.gz"
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            snp_i = header.index("SNP")
            for line in f:
                snps.add(line.rstrip("\n").split("\t")[snp_i])
    return snps


def write_sumstats(path: Path, frame_iter) -> int:
    tmp_path = path.with_suffix("")
    n = 0
    with gzip.open(path, "wt") as f:
        f.write("SNP\tA1\tA2\tZ\tP\tN\n")
        for frame in frame_iter:
            if frame.empty:
                continue
            frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["SNP", "A1", "A2", "Z", "P", "N"])
            frame = frame.sort_values("P").drop_duplicates("SNP", keep="first")
            frame[["SNP", "A1", "A2", "Z", "P", "N"]].to_csv(f, sep="\t", header=False, index=False)
            n += len(frame)
    tmp_path.unlink(missing_ok=True)
    return n


def munge_gwas(label: str, path: Path, short_name: str, mapping: pd.DataFrame) -> tuple[Path, int]:
    out_path = OUT / f"{short_name}.unique_allele.sumstats.gz"
    if out_path.exists():
        return out_path, count_gzip_rows(out_path)

    def chunks():
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
            yield pd.DataFrame({
                "SNP": merged["SNP"],
                "A1": merged["alt"],
                "A2": merged["ref"],
                "Z": z,
                "P": pd.to_numeric(merged["P"], errors="coerce"),
                "N": pd.to_numeric(merged["N"], errors="coerce"),
            })

    print(f"Munging {label}")
    return out_path, write_sumstats(out_path, chunks())


def munge_hgi(label: str, path: Path, short_name: str, ld_snps: set[str], mapping: pd.DataFrame) -> tuple[Path, int]:
    out_path = OUT / f"{short_name}.hm3_unique_allele.sumstats.gz"
    if out_path.exists():
        return out_path, count_gzip_rows(out_path)

    allele_map = mapping[["SNP", "ref", "alt"]].drop_duplicates("SNP")

    def chunks():
        usecols = ["rsid", "REF", "ALT", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p", "all_inv_var_meta_effective"]
        for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=1_000_000):
            chunk["REF"] = chunk["REF"].astype(str).str.upper()
            chunk["ALT"] = chunk["ALT"].astype(str).str.upper()
            chunk = chunk[
                chunk["rsid"].isin(ld_snps)
                & chunk["REF"].isin(["A", "C", "G", "T"])
                & chunk["ALT"].isin(["A", "C", "G", "T"])
            ]
            chunk = chunk.merge(allele_map, left_on="rsid", right_on="SNP", how="inner")
            same = (chunk["REF"] == chunk["ref"]) & (chunk["ALT"] == chunk["alt"])
            flip = (chunk["REF"] == chunk["alt"]) & (chunk["ALT"] == chunk["ref"])
            chunk = chunk[same | flip]
            beta = pd.to_numeric(chunk["all_inv_var_meta_beta"], errors="coerce")
            se = pd.to_numeric(chunk["all_inv_var_meta_sebeta"], errors="coerce")
            yield pd.DataFrame({
                "SNP": chunk["rsid"],
                "A1": chunk["ALT"],
                "A2": chunk["REF"],
                "Z": beta / se,
                "P": pd.to_numeric(chunk["all_inv_var_meta_p"], errors="coerce"),
                "N": pd.to_numeric(chunk["all_inv_var_meta_effective"], errors="coerce"),
            })

    print(f"Munging {label}")
    return out_path, write_sumstats(out_path, chunks())


def count_gzip_rows(path: Path) -> int:
    with gzip.open(path, "rt") as f:
        return max(sum(1 for _ in f) - 1, 0)


def run_ldsc(args: list[str], log_prefix: Path) -> str:
    cmd = [
        "uv", "tool", "run", "--python", "3.11",
        "--with", "pandas<2",
        "--with", "numpy<2",
        "--from", "ldsc", "ldsc.py",
        *args,
        "--out", str(log_prefix),
    ]
    subprocess.run(cmd, cwd=ROOT, check=False)
    return log_prefix.with_suffix(".log").read_text()


def parse_h2_log(text: str) -> dict[str, float | int]:
    patterns = {
        "ldsc_snp_count": r"After merging with regression SNP LD, ([0-9]+) SNPs remain",
        "intercept": r"Intercept: ([^ ]+) \(([^)]+)\)",
    }
    n_match = re.search(patterns["ldsc_snp_count"], text)
    i_match = re.search(patterns["intercept"], text)
    return {
        "ldsc_snp_count": int(n_match.group(1)) if n_match else math.nan,
        "intercept": float(i_match.group(1)) if i_match else math.nan,
        "intercept_se": float(i_match.group(2)) if i_match else math.nan,
    }


def parse_rg_log(text: str) -> list[dict[str, str]]:
    rows = []
    marker = "Summary of Genetic Correlation Results"
    if marker not in text:
        return rows
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
    if not HM3_MAP.exists():
        raise SystemExit(f"Missing {HM3_MAP}. Run Draft/build_full_hm3_map.R first.")

    mapping = load_map()
    ld_snps = load_ld_snps()

    sumstats: dict[str, tuple[str, Path, int]] = {}
    for label, path, short_name in GWAS_SCANS:
        out, n = munge_gwas(label, path, short_name, mapping)
        sumstats[short_name] = (label, out, n)
        print(f"  {short_name}: {n:,} HapMap3 SNPs")

    for label, path, short_name in HGI_TRAITS:
        out, n = munge_hgi(label, path, short_name, ld_snps, mapping)
        sumstats[short_name] = (label, out, n)
        print(f"  {short_name}: {n:,} SNPs before LDSC merge")

    qc_rows = []
    for label, path, short_name in GWAS_SCANS:
        print(f"Running LDSC h2 for {label}")
        log = run_ldsc(
            [
                "--h2", str(sumstats[short_name][1]),
                "--ref-ld-chr", LD_PREFIX,
                "--w-ld-chr", LD_PREFIX,
            ],
            OUT / f"ldsc_h2_{short_name}",
        )
        parsed = parse_h2_log(log)
        full_n, median_chisq, lambda_gc = lambda_gc_from_full(path)
        qc_rows.append({
            "Scan": label,
            "Full-summary variants": f"{full_n:,}",
            "Median chi-square": fmt(median_chisq, 4),
            "lambda_GC": fmt(lambda_gc, 4),
            "LDSC SNPs": f"{int(parsed['ldsc_snp_count']):,}",
            "LDSC intercept": fmt(float(parsed["intercept"]), 4),
            "LDSC intercept SE": fmt(float(parsed["intercept_se"]), 4),
        })

    qc_fields = ["Scan", "Full-summary variants", "Median chi-square", "lambda_GC", "LDSC SNPs", "LDSC intercept", "LDSC intercept SE"]
    write_csv(OUT / "supplementary_table_lambda_gc_ldsc_intercepts.csv", qc_rows, qc_fields)
    write_tex(
        OUT / "supplementary_table_lambda_gc_ldsc_intercepts.tex",
        qc_rows,
        qc_fields,
        "Genomic inflation and LD-score regression intercepts for the eight primary GWAS scans.",
        "tab:lambda-gc-ldsc-intercepts",
    )

    rg_rows = []
    for label, _, short_name in GWAS_SCANS:
        print(f"Running LDSC rg for {label}")
        rg_arg = ",".join([
            str(sumstats[short_name][1]),
            str(sumstats["HGI_C2_susceptibility"][1]),
            str(sumstats["HGI_B2_severity"][1]),
        ])
        log = run_ldsc(
            [
                "--rg", rg_arg,
                "--ref-ld-chr", LD_PREFIX,
                "--w-ld-chr", LD_PREFIX,
            ],
            OUT / f"ldsc_rg_{short_name}_vs_HGI",
        )
        parsed_rows = parse_rg_log(log)
        if not parsed_rows:
            for trait_label in ["HGI C2 susceptibility", "HGI B2 hospitalisation/severity"]:
                rg_rows.append({
                    "PCC scan": label,
                    "External trait": trait_label,
                    "rg": "",
                    "rg SE": "",
                    "z": "",
                    "rg p": "",
                    "Status": "Unavailable: LDSC rg failed, likely because h2 was out of bounds or unstable",
                })
            continue
        for row in parsed_rows:
            p2 = Path(row.get("p2", "")).name.replace(".sumstats.gz", "")
            trait_label = {
                "HGI_C2_susceptibility": "HGI C2 susceptibility",
                "HGI_B2_severity": "HGI B2 hospitalisation/severity",
                "HGI_C2_susceptibility.hm3_unique_allele": "HGI C2 susceptibility",
                "HGI_B2_severity.hm3_unique_allele": "HGI B2 hospitalisation/severity",
            }.get(p2, p2)
            rg = parse_optional_float(row["rg"])
            se = parse_optional_float(row["se"])
            z = parse_optional_float(row["z"])
            p = parse_optional_float(row["p"])
            status = "OK"
            if rg is None:
                status = "Unavailable: LDSC h2 out of bounds"
            elif abs(rg) > 1:
                status = "Unstable: LDSC reported rg out of bounds"
            rg_rows.append({
                "PCC scan": label,
                "External trait": trait_label,
                "rg": fmt(rg, 4),
                "rg SE": fmt(se, 4),
                "z": fmt(z, 3),
                "rg p": fmt_p(p),
                "Status": status,
            })

    rg_fields = ["PCC scan", "External trait", "rg", "rg SE", "z", "rg p", "Status"]
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
