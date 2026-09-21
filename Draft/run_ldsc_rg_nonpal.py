#!/usr/bin/env python3
"""Run LDSC rg on non-palindromic HGI-position harmonised sumstats."""

from __future__ import annotations

import csv
import gzip
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Draft" / "ldsc_hgi_position_outputs"
LD_PREFIX = str(ROOT / "Draft" / "ldsc" / "eur_w_ld_chr") + "/"

GWAS_SCANS = [
    ("PCC vs responder controls", "PCC_responder"),
    ("Pauci-symptomatic PCC vs responder controls", "C1_responder"),
    ("Fatigue/sleep-dominated PCC vs responder controls", "C2_responder"),
    ("Multi-symptomatic PCC vs responder controls", "C3_responder"),
    ("PCC vs infected non-PCC controls", "PCC_infected_nonPCC"),
    ("Pauci-symptomatic PCC vs infected non-PCC controls", "C1_infected_nonPCC"),
    ("Fatigue/sleep-dominated PCC vs infected non-PCC controls", "C2_infected_nonPCC"),
    ("Multi-symptomatic PCC vs infected non-PCC controls", "C3_infected_nonPCC"),
]

HGI_TRAITS = [
    ("HGI C2 susceptibility", "HGI_C2_susceptibility"),
    ("HGI B2 hospitalisation/severity", "HGI_B2_severity"),
]


def is_pal(a1: str, a2: str) -> bool:
    return {a1, a2} in [{"A", "T"}, {"C", "G"}]


def nonpal_file(short: str) -> Path:
    return OUT / f"{short}.hgi_position.nonpal.sumstats.gz"


def filter_nonpal(short: str) -> tuple[Path, int]:
    src = OUT / f"{short}.hgi_position.sumstats.gz"
    dst = nonpal_file(short)
    if dst.exists():
        with gzip.open(dst, "rt") as f:
            return dst, max(sum(1 for _ in f) - 1, 0)
    n = 0
    with gzip.open(src, "rt") as inp, gzip.open(dst, "wt") as out:
        header = inp.readline().rstrip("\n").split("\t")
        out.write("\t".join(header) + "\n")
        a1_i = header.index("A1")
        a2_i = header.index("A2")
        for line in inp:
            parts = line.rstrip("\n").split("\t")
            if not is_pal(parts[a1_i], parts[a2_i]):
                out.write(line)
                n += 1
    return dst, n


def run_ldsc_rg(short: str) -> str:
    rg_arg = ",".join([
        str(nonpal_file(short)),
        str(nonpal_file("HGI_C2_susceptibility")),
        str(nonpal_file("HGI_B2_severity")),
    ])
    out_prefix = OUT / f"ldsc_rg_nonpal_{short}_vs_HGI"
    cmd = [
        "uv", "tool", "run", "--python", "3.11",
        "--with", "pandas<2", "--with", "numpy<2",
        "--from", "ldsc", "ldsc.py",
        "--rg", rg_arg,
        "--ref-ld-chr", LD_PREFIX,
        "--w-ld-chr", LD_PREFIX,
        "--out", str(out_prefix),
    ]
    subprocess.run(cmd, cwd=ROOT, check=False)
    return out_prefix.with_suffix(".log").read_text()


def parse_optional_float(value: str):
    if value in {"NA", "nan", "NaN", ""}:
        return None
    try:
        value = float(value)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def fmt(value, digits=4):
    if value is None:
        return ""
    try:
        value = float(value)
    except ValueError:
        return str(value)
    return f"{value:.{digits}f}" if math.isfinite(value) else ""


def fmt_p(value):
    if value is None:
        return ""
    value = float(value)
    return f"{value:.2e}" if value < 1e-3 else f"{value:.4f}"


def parse_rg_log(text: str) -> list[dict]:
    marker = "Summary of Genetic Correlation Results"
    if marker not in text:
        return []
    lines = text.split(marker, 1)[1].strip().splitlines()
    if len(lines) < 2:
        return []
    header = re.split(r"\s+", lines[0].strip())
    rows = []
    for line in lines[1:]:
        if not line.strip() or line.startswith("Analysis finished"):
            break
        parts = re.split(r"\s+", line.strip())
        if len(parts) == len(header):
            rows.append(dict(zip(header, parts)))
    return rows


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


def write_tex(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w") as f:
        f.write("\\begin{table}[htbp]\n\\centering\n")
        f.write("\\caption{LD-score regression genetic correlations after HGI-position harmonisation and palindromic SNP removal.}\n")
        f.write("\\label{tab:ldsc-rg-hgi-nonpal}\n")
        f.write("\\resizebox{\\linewidth}{!}{%\n")
        f.write("\\begin{tabular}{" + "l" * len(fields) + "}\n\\toprule\n")
        f.write(" & ".join(latex_escape(x) for x in fields) + " \\\\\n\\midrule\n")
        for row in rows:
            f.write(" & ".join(latex_escape(row.get(x, "")) for x in fields) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}%\n}\n\\end{table}\n")


def main() -> None:
    counts = {}
    for _, short in GWAS_SCANS:
        _, counts[short] = filter_nonpal(short)
    for _, short in HGI_TRAITS:
        _, counts[short] = filter_nonpal(short)
    print(counts)

    rows = []
    for label, short in GWAS_SCANS:
        log = run_ldsc_rg(short)
        parsed = parse_rg_log(log)
        if not parsed:
            for hgi_label, _ in HGI_TRAITS:
                rows.append({
                    "PCC scan": label,
                    "External trait": hgi_label,
                    "Shared LDSC SNPs": "",
                    "rg": "",
                    "rg SE": "",
                    "z": "",
                    "rg p": "",
                    "Status": "Unavailable",
                })
            continue
        for row in parsed:
            p2 = Path(row["p2"]).name.replace(".hgi_position.nonpal.sumstats.gz", "")
            hgi_label = dict((short, label) for label, short in HGI_TRAITS).get(p2, p2)
            rg = parse_optional_float(row["rg"])
            se = parse_optional_float(row["se"])
            z = parse_optional_float(row["z"])
            p = parse_optional_float(row["p"])
            status = "OK"
            if rg is None:
                status = "Unavailable"
            elif abs(rg) > 1:
                status = "Unstable: LDSC reported rg out of bounds"
            rows.append({
                "PCC scan": label,
                "External trait": hgi_label,
                "Shared LDSC SNPs": f"{min(counts[short], counts[p2]):,}",
                "rg": fmt(rg),
                "rg SE": fmt(se),
                "z": fmt(z, 3),
                "rg p": fmt_p(p),
                "Status": status,
            })
    fields = ["PCC scan", "External trait", "Shared LDSC SNPs", "rg", "rg SE", "z", "rg p", "Status"]
    write_csv(OUT / "supplementary_table_ldsc_rg_hgi_nonpal.csv", rows, fields)
    write_tex(OUT / "supplementary_table_ldsc_rg_hgi_nonpal.tex", rows, fields)
    print(f"Wrote non-palindromic rg table to {OUT}")


if __name__ == "__main__":
    main()
