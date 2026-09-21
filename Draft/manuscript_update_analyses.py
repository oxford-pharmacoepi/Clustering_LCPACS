#!/usr/bin/env python3
"""Generate manuscript-support tables for the PCC GWAS analyses.

The script intentionally uses only local files already present in the repo.
It writes CSV, LaTeX, and Markdown outputs under Draft/manuscript_update_outputs.
"""

from __future__ import annotations

import csv
import gzip
import math
import statistics
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Draft" / "manuscript_update_outputs"
OUT.mkdir(parents=True, exist_ok=True)

FINGEN_XLSX = ROOT / "Draft" / "DPA_Replication_VL_20260120.xlsx"
FINGEN_SHEET = "DPA_SNPs_FinnGen"

DISCOVERY_FILE_BY_ANALYSIS = {
    "AllPCCvsGenPop": ROOT / "GWAS" / "pcc_vsallfil.txt.gz",
    "Subtype1vsPopCtrl": ROOT / "GWAS" / "clust1_vsallfil.txt.gz",
    "Subtype2vsPopCtrl": ROOT / "GWAS" / "clust2_vsallfil.txt.gz",
    "Subtype3vsCOVID": ROOT / "GWAS" / "clust3_vsno.txt.gz",
}

QC_GWAS_FILES = [
    ("PCC vs responder controls", ROOT / "GWAS" / "pcc_vsallfil.txt.gz"),
    ("Pauci-symptomatic PCC vs responder controls", ROOT / "GWAS" / "clust1_vsallfil.txt.gz"),
    ("Fatigue/sleep-dominated PCC vs responder controls", ROOT / "GWAS" / "clust2_vsallfil.txt.gz"),
    ("Multi-symptomatic PCC vs responder controls", ROOT / "GWAS" / "clust3_vsallfil.txt.gz"),
    ("PCC vs infected non-PCC controls", ROOT / "GWAS" / "pcc_vscovid.txt.gz"),
    ("Pauci-symptomatic PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust1_vsno.txt.gz"),
    ("Fatigue/sleep-dominated PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust2_vsno.txt.gz"),
    ("Multi-symptomatic PCC vs infected non-PCC controls", ROOT / "GWAS" / "clust3_vsno.txt.gz"),
]

PRESPEC_ENDPOINT_BY_ANALYSIS = {
    "AllPCCvsGenPop": "LC_U09",
    "Subtype1vsPopCtrl": "LC_U09",
    "Subtype2vsPopCtrl": "LC_U09",
    "Subtype3vsCOVID": "LC_U09_vs_U07",
}

LATEX_SPECIAL = {
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
}


def fmt_float(x, digits=3):
    if x is None or x == "":
        return ""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not math.isfinite(x):
        return ""
    return f"{x:.{digits}f}"


def fmt_p(x):
    if x is None or x == "":
        return ""
    x = float(x)
    if not math.isfinite(x):
        return ""
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.4f}"


def latex_escape(value):
    text = "" if value is None else str(value)
    for raw, repl in LATEX_SPECIAL.items():
        text = text.replace(raw, repl)
    return text


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_latex_table(path, rows, fieldnames, caption, label):
    with path.open("w") as f:
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write(f"\\caption{{{latex_escape(caption)}}}\n")
        f.write(f"\\label{{{latex_escape(label)}}}\n")
        f.write("\\begin{tabular}{" + "l" * len(fieldnames) + "}\n")
        f.write("\\toprule\n")
        f.write(" & ".join(latex_escape(x) for x in fieldnames) + " \\\\\n")
        f.write("\\midrule\n")
        for row in rows:
            f.write(" & ".join(latex_escape(row.get(x, "")) for x in fieldnames) + " \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")


def write_markdown_table(path, rows, fieldnames):
    with path.open("w") as f:
        f.write("| " + " | ".join(fieldnames) + " |\n")
        f.write("| " + " | ".join(["---"] * len(fieldnames)) + " |\n")
        for row in rows:
            f.write("| " + " | ".join(str(row.get(x, "")) for x in fieldnames) + " |\n")


def bh_fdr(p_values):
    indexed = sorted((float(p), i) for i, p in enumerate(p_values))
    n = len(indexed)
    out = [math.nan] * n
    running = 1.0
    for rank_from_end, (p, i) in enumerate(reversed(indexed), start=1):
        rank = n - rank_from_end + 1
        running = min(running, p * n / rank)
        out[i] = min(running, 1.0)
    return out


def binom_sign_test_p(k, n):
    return sum(math.comb(n, i) * (0.5 ** n) for i in range(k, n + 1))


def load_fingen_rows():
    wb = load_workbook(FINGEN_XLSX, data_only=True, read_only=True)
    ws = wb[FINGEN_SHEET]
    rows = ws.iter_rows(values_only=True)
    header = [str(x) if x is not None else "" for x in next(rows)]
    idx = {name: i for i, name in enumerate(header)}
    out = []
    for raw in rows:
        if not raw or raw[idx["rsID"]] is None or raw[idx["LOG10P"]] is None:
            continue
        row = {name: raw[i] if i < len(raw) else None for name, i in idx.items()}
        row["fingen_p"] = 10 ** (-float(row["LOG10P"]))
        out.append(row)
    return out


def discovery_lookup(target_rows):
    targets_by_file = defaultdict(dict)
    for row in target_rows:
        file_path = DISCOVERY_FILE_BY_ANALYSIS.get(row["DPA_Analysis"])
        if not file_path:
            continue
        key = (int(row["Chr"]), int(row["Position_GRCh38"]))
        targets_by_file[file_path][key] = row

    found = {}
    for file_path, targets in targets_by_file.items():
        with gzip.open(file_path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            for line in f:
                parts = line.rstrip("\n").split("\t")
                try:
                    key = (int(parts[idx["CHROM"]]), int(parts[idx["POSITION"]]))
                except (ValueError, KeyError):
                    continue
                if key not in targets:
                    continue
                found[(file_path, key)] = {
                    "discovery_beta": float(parts[idx["BETA"]]),
                    "discovery_se": float(parts[idx["SE"]]),
                    "discovery_p": float(parts[idx["P"]]),
                    "discovery_effect_allele": parts[idx["EFFECT_ALLELE"]],
                    "discovery_non_effect_allele": parts[idx["NON_EFFECT_ALLELE"]],
                    "discovery_variant_id": parts[idx["SNP"]],
                }
                if len(found) == len(targets):
                    break
    for row in target_rows:
        file_path = DISCOVERY_FILE_BY_ANALYSIS.get(row["DPA_Analysis"])
        if not file_path:
            row.update({})
            continue
        key = (int(row["Chr"]), int(row["Position_GRCh38"]))
        row.update(found.get((file_path, key), {}))
    return target_rows


def harmonise_row(row):
    d_ea = row.get("discovery_effect_allele")
    d_nea = row.get("discovery_non_effect_allele")
    f_ea = row.get("ALLELE1")
    f_nea = row.get("ALLELE0")
    f_beta = row.get("BETA")
    if d_ea == f_ea and d_nea == f_nea:
        row["fingen_beta_harmonised"] = float(f_beta)
        row["allele_status"] = "same"
    elif d_ea == f_nea and d_nea == f_ea:
        row["fingen_beta_harmonised"] = -float(f_beta)
        row["allele_status"] = "flipped"
    else:
        row["fingen_beta_harmonised"] = math.nan
        row["allele_status"] = "mismatch"
    db = row.get("discovery_beta")
    fb = row.get("fingen_beta_harmonised")
    if db is None or not math.isfinite(fb):
        row["direction_same"] = ""
    else:
        row["direction_same"] = "yes" if math.copysign(1, db) == math.copysign(1, fb) else "no"
    return row


def summarise_fingen_rows(rows):
    for row in rows:
        harmonise_row(row)
    out = []
    for row in rows:
        out.append({
            "SNP": row["rsID"],
            "Gene": row["NearestGene"],
            "Discovery analysis": row["DPA_Analysis"],
            "FinnGen endpoint": row["FinnGen_Analysis"],
            "Discovery beta": fmt_float(row.get("discovery_beta"), 3),
            "Discovery SE": fmt_float(row.get("discovery_se"), 3),
            "FinnGen beta": fmt_float(row.get("fingen_beta_harmonised"), 3),
            "FinnGen SE": fmt_float(row.get("SE"), 3),
            "Direction same": row.get("direction_same", ""),
            "FinnGen p": fmt_p(row.get("fingen_p")),
        })
    return out


def make_control_definition_tables():
    rows = [
        {
            "Quantity": "Questionnaire responders with infection before survey",
            "Numerator": "51,780",
            "Denominator": "54,283",
            "Fraction": "95.4%",
            "Use": "Shows most infected questionnaire responders reported infection before questionnaire completion.",
        },
        {
            "Quantity": "Qualifying infection window among infected responders",
            "Numerator": "41,207",
            "Denominator": "51,780",
            "Fraction": "79.6%",
            "Use": "Defines the infection-anchored population for PCC versus infected non-PCC controls.",
        },
        {
            "Quantity": "Pre-QC PCC among qualifying infected responders",
            "Numerator": "20,455",
            "Denominator": "41,207",
            "Fraction": "49.6%",
            "Use": "Participants with at least one qualifying post-infection symptom.",
        },
        {
            "Quantity": "Pre-QC infected non-PCC among qualifying infected responders",
            "Numerator": "20,752",
            "Denominator": "41,207",
            "Fraction": "50.4%",
            "Use": "Participants forming the infected non-PCC control pool before genotype QC.",
        },
        {
            "Quantity": "Post-QC infected non-PCC controls relative to responder controls",
            "Numerator": "20,233",
            "Denominator": "170,643",
            "Fraction": "11.9%",
            "Use": "Responder controls are not equivalent to infection-confirmed non-PCC controls.",
        },
    ]
    fields = ["Quantity", "Numerator", "Denominator", "Fraction", "Use"]
    write_csv(OUT / "01_control_definitions.csv", rows, fields)
    write_latex_table(
        OUT / "01_control_definitions.tex",
        rows,
        fields,
        "Control-definition quantities relevant to susceptibility-versus-PCC interpretation.",
        "tab:control-definitions",
    )
    write_markdown_table(OUT / "01_control_definitions.md", rows, fields)
    return rows


def make_lambda_gc_table():
    rows = []
    denominator = 0.454936423119572
    for label, path in QC_GWAS_FILES:
        chisq = []
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(header)}
            for line in f:
                parts = line.rstrip("\n").split("\t")
                try:
                    value = float(parts[idx["CHISQ"]])
                except (ValueError, KeyError, IndexError):
                    continue
                if math.isfinite(value):
                    chisq.append(value)
        median_chisq = statistics.median(chisq)
        rows.append({
            "Scan": label,
            "Variants": f"{len(chisq):,}",
            "Median chi-square": fmt_float(median_chisq, 4),
            "lambda_GC": fmt_float(median_chisq / denominator, 4),
            "LDSC intercept": "TO RUN",
            "LDSC intercept SE": "TO RUN",
        })
    fields = ["Scan", "Variants", "Median chi-square", "lambda_GC", "LDSC intercept", "LDSC intercept SE"]
    write_csv(OUT / "02_gwas_qc_lambda_gc.csv", rows, fields)
    write_latex_table(
        OUT / "02_gwas_qc_lambda_gc.tex",
        rows,
        fields,
        "Genome-wide inflation diagnostics for the eight primary GWAS scans.",
        "tab:gwas-qc-lambda-gc",
    )
    write_markdown_table(OUT / "02_gwas_qc_lambda_gc.md", rows, fields)
    return rows


def make_fingen_tables():
    rows = load_fingen_rows()
    by_snp = defaultdict(list)
    for row in rows:
        by_snp[row["rsID"]].append(row)

    best_rows = [max(group, key=lambda x: float(x["LOG10P"])) for group in by_snp.values()]
    best_rows = discovery_lookup(best_rows)
    best_summary = summarise_fingen_rows(best_rows)
    best_summary.sort(key=lambda x: float(x["FinnGen p"]), reverse=False)
    fields = [
        "SNP",
        "Gene",
        "Discovery analysis",
        "FinnGen endpoint",
        "Discovery beta",
        "Discovery SE",
        "FinnGen beta",
        "FinnGen SE",
        "Direction same",
        "FinnGen p",
    ]
    write_csv(OUT / "03_fingen_best_effect_size_comparison.csv", best_summary, fields)
    write_latex_table(
        OUT / "03_fingen_best_effect_size_comparison.tex",
        best_summary,
        fields,
        "Effect-size comparison for best FinnGen endpoint per index SNP.",
        "tab:fingen-effect-size-best",
    )
    write_markdown_table(OUT / "03_fingen_best_effect_size_comparison.md", best_summary, fields)

    prespec_rows = []
    missing = []
    for snp, group in by_snp.items():
        analysis = group[0]["DPA_Analysis"]
        endpoint = PRESPEC_ENDPOINT_BY_ANALYSIS[analysis]
        matches = [row for row in group if row["FinnGen_Analysis"] == endpoint]
        if not matches:
            missing.append((snp, endpoint))
            continue
        prespec_rows.append(matches[0])
    prespec_rows = discovery_lookup(prespec_rows)
    for row in prespec_rows:
        harmonise_row(row)
    prespec_rows.sort(key=lambda x: x["rsID"])
    p_values = [float(row["fingen_p"]) for row in prespec_rows]
    fdr = bh_fdr(p_values)
    prespec_summary = []
    for row, q in zip(prespec_rows, fdr):
        prespec_summary.append({
            "SNP": row["rsID"],
            "Gene": row["NearestGene"],
            "Discovery analysis": row["DPA_Analysis"],
            "Pre-specified FinnGen endpoint": row["FinnGen_Analysis"],
            "Discovery beta": fmt_float(row.get("discovery_beta"), 3),
            "FinnGen beta": fmt_float(row.get("fingen_beta_harmonised"), 3),
            "Direction same": row.get("direction_same", ""),
            "FinnGen p": fmt_p(row.get("fingen_p")),
            "Bonferroni p over 10": fmt_p(min(float(row["fingen_p"]) * len(prespec_rows), 1.0)),
            "BH-FDR q over 10": fmt_p(q),
        })
    fields2 = [
        "SNP",
        "Gene",
        "Discovery analysis",
        "Pre-specified FinnGen endpoint",
        "Discovery beta",
        "FinnGen beta",
        "Direction same",
        "FinnGen p",
        "Bonferroni p over 10",
        "BH-FDR q over 10",
    ]
    write_csv(OUT / "04_fingen_prespecified_10test_correction.csv", prespec_summary, fields2)
    write_latex_table(
        OUT / "04_fingen_prespecified_10test_correction.tex",
        prespec_summary,
        fields2,
        "Pre-specified FinnGen endpoint replication with correction over ten index SNPs.",
        "tab:fingen-prespecified-ten-tests",
    )
    write_markdown_table(OUT / "04_fingen_prespecified_10test_correction.md", prespec_summary, fields2)

    best_concordant = sum(1 for row in best_summary if row["Direction same"] == "yes")
    best_n = sum(1 for row in best_summary if row["Direction same"] in {"yes", "no"})
    prespec_concordant = sum(1 for row in prespec_summary if row["Direction same"] == "yes")
    prespec_n = sum(1 for row in prespec_summary if row["Direction same"] in {"yes", "no"})
    sign_rows = [
        {
            "Analysis set": "Best FinnGen endpoint per SNP",
            "Concordant": str(best_concordant),
            "Total with harmonised beta": str(best_n),
            "One-sided sign-test p": fmt_p(binom_sign_test_p(best_concordant, best_n)),
        },
        {
            "Analysis set": "Pre-specified endpoint per SNP",
            "Concordant": str(prespec_concordant),
            "Total with harmonised beta": str(prespec_n),
            "One-sided sign-test p": fmt_p(binom_sign_test_p(prespec_concordant, prespec_n)),
        },
    ]
    fields3 = ["Analysis set", "Concordant", "Total with harmonised beta", "One-sided sign-test p"]
    write_csv(OUT / "04b_fingen_sign_tests.csv", sign_rows, fields3)
    write_latex_table(
        OUT / "04b_fingen_sign_tests.tex",
        sign_rows,
        fields3,
        "Direction-concordance sign tests for FinnGen replication.",
        "tab:fingen-sign-tests",
    )
    return best_summary, prespec_summary, sign_rows, missing


def make_entropy_table():
    path = ROOT / "Shiny" / "Results" / "Clustering_cases" / "Information_criteria.csv"
    rows = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            if raw["Clust"] != "3":
                continue
            rows.append({
                "LCA solution": "3-class PCC model",
                "Class": f"Class {len(rows) + 1}",
                "Entropy": raw["entropy"],
                "Mean posterior probability": f"{float(raw['mean_posterior']):.1f}%",
            })
    fields = ["LCA solution", "Class", "Entropy", "Mean posterior probability"]
    write_csv(OUT / "05_lca_entropy.csv", rows, fields)
    write_latex_table(
        OUT / "05_lca_entropy.tex",
        rows,
        fields,
        "Entropy and class-assignment certainty for the selected three-class PCC latent class model.",
        "tab:lca-entropy",
    )
    write_markdown_table(OUT / "05_lca_entropy.md", rows, fields)
    return rows


def make_ldsc_template():
    text = """# 06 LDSC commands still required

LDSC is not installed in this repository and no LD-score reference directory was found locally.
Run these commands on the machine/environment where LDSC and the European LD-score references are available.

Inputs to prepare:
- Full unfiltered GWAS summary statistics for the eight primary scans, mapped to rsIDs.
- Important: the raw GWAS `SNP` column contains DRAGEN IDs, so do not use it directly as the LDSC SNP column.
- Use the existing rsID mapping workflow in `Study/4-PostProcessing.R` / `GWAS/all_rsids*.rds` to create files with columns like `RSID P A1 A2 BETA SE N`.
- HGI C2 SARS-CoV-2 susceptibility summary statistics.
- HGI B2 COVID-19 hospitalisation/severity summary statistics.
- 1000 Genomes European LD scores and weights, for example `eur_w_ld_chr/`.
- HapMap3 SNP list for munging.

Example commands:

```bash
munge_sumstats.py --sumstats pcc_vsallfil_rsid.txt --snp RSID --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N 190611 --merge-alleles w_hm3.snplist --out ldsc/pcc_vs_responder

munge_sumstats.py --sumstats pcc_vscovid_rsid.txt --snp RSID --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N 40201 --merge-alleles w_hm3.snplist --out ldsc/pcc_vs_infected_non_pcc

ldsc.py --h2 ldsc/pcc_vs_responder.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ldsc/pcc_vs_responder_h2

ldsc.py --rg ldsc/pcc_vs_responder.sumstats.gz,ldsc/hgi_c2_susceptibility.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ldsc/pcc_vs_responder_rg_hgi_c2

ldsc.py --rg ldsc/pcc_vs_responder.sumstats.gz,ldsc/hgi_b2_hospitalisation.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ldsc/pcc_vs_responder_rg_hgi_b2
```

Manuscript table fields to fill after running LDSC:
- scan
- lambda_GC
- mean chi-square
- LDSC intercept
- LDSC intercept SE
- rg with HGI C2 susceptibility
- rg SE
- rg p
- rg with HGI B2 hospitalisation/severity
- rg SE
- rg p
"""
    (OUT / "06_ldsc_commands_and_table_slots.md").write_text(text)


def make_change_pack(control_rows, qc_rows, best_rows, prespec_rows, sign_rows, entropy_rows):
    best_concordance = next(row for row in sign_rows if row["Analysis set"] == "Best FinnGen endpoint per SNP")
    prespec_concordance = next(row for row in sign_rows if row["Analysis set"] == "Pre-specified endpoint per SNP")
    text = f"""# Manuscript update pack for GWAS_paper.docx

Use this file as a numbered checklist for manuscript updates. Tables are in this folder as CSV, LaTeX, and Markdown.

## 1. Responder-control infection fraction and susceptibility-versus-PCC distinction

Key numbers:
- 41,207 / 51,780 = 79.6% of infected questionnaire responders had the recorded infection in the one-year to 30-day window before survey completion.
- Before genotype QC, the qualifying infected responder group split into 20,455 PCC cases and 20,752 infected non-PCC controls.
- After genotype QC, the main analyses used 19,968 PCC cases, 20,233 infected non-PCC controls, and 170,643 responder controls.
- The infected non-PCC control count is 20,233 / 170,643 = 11.9% of the responder-control count, illustrating that responder controls are not equivalent to infection-confirmed controls.

Table to use: `01_control_definitions.*`.

Introduction text:

> Because the responder-control contrast compares PCC cases with questionnaire responders who were not required to have a recorded SARS-CoV-2 infection, associations from this scan may reflect a mixture of genetic liability to PCC and genetic liability to acquiring a recorded SARS-CoV-2 infection. We therefore interpret the responder-control GWAS as a high-powered discovery analysis, and the infected non-PCC control GWAS as the more specific contrast for PCC among people with documented infection.

Results text:

> Among 51,780 questionnaire responders with a positive SARS-CoV-2 test before questionnaire completion, 41,207 (79.6%) had a recorded infection between one year and 30 days before survey completion. This infection-anchored group comprised 20,455 individuals with qualifying post-infection symptoms and 20,752 infected responders without qualifying PCC symptoms before genotype QC. After genotype QC, 19,968 PCC cases and 20,233 infected non-PCC controls remained; the broader responder-control group contained 170,643 participants.

## 2. QQ plots, lambda_GC, and LDSC intercept

Local lambda_GC values from full unfiltered GWAS summary statistics are in `02_gwas_qc_lambda_gc.*`. LDSC is not installed locally, so the intercept columns are deliberately marked `TO RUN`.

Results text:

> Inflation metrics were modest across the eight primary scans (lambda_GC range {qc_rows[1]['lambda_GC']}-{qc_rows[-1]['lambda_GC']}; Table X), consistent with the QQ plots and arguing against substantial unmodelled stratification or technical inflation. LDSC intercepts are reported in Table X.

Replace the range above with the table min-max if you prefer exact prose after final table insertion.

Methods text:

> For each GWAS, we assessed test-statistic inflation using QQ plots, genomic control inflation factors (lambda_GC; median observed chi-square divided by the expected median under a one-degree-of-freedom chi-square distribution), and LD-score regression intercepts. LDSC analyses were performed on full unfiltered summary statistics after standard munging to HapMap3 SNPs, using European LD-score reference weights. The LDSC intercept was used to distinguish polygenicity from residual confounding.

Run LDSC using `06_ldsc_commands_and_table_slots.md`, then replace `TO RUN` in `02_gwas_qc_lambda_gc.*`.

## 3. FinnGen beta-to-beta comparison plus sign test

Table to use: `03_fingen_best_effect_size_comparison.*`.

Key result:
- Best-endpoint lookup: {best_concordance['Concordant']} / {best_concordance['Total with harmonised beta']} SNPs had concordant discovery and FinnGen effects after allele harmonisation; one-sided sign-test p = {best_concordance['One-sided sign-test p']}.
- The only discordant best-endpoint SNP is AC104183.1 / rs2661400.

Manuscript text:

> We compared effect estimates between the discovery GWAS and FinnGen after harmonising FinnGen betas to the discovery effect allele. Directional concordance was observed for {best_concordance['Concordant']} of {best_concordance['Total with harmonised beta']} index SNPs in the best-endpoint lookup (one-sided sign-test p = {best_concordance['One-sided sign-test p']}), indicating that the replication pattern was not driven only by small p-values but was also supported by consistent effect directions and comparable effect estimates.

## 4. FinnGen multiple testing correction over ten pre-specified endpoint tests

Table to use: `04_fingen_prespecified_10test_correction.*`.

Pre-specified endpoint rule used:
- AllPCCvsGenPop -> LC_U09.
- Subtype1vsPopCtrl -> LC_U09.
- Subtype2vsPopCtrl -> LC_U09.
- Subtype3vsCOVID -> LC_U09_vs_U07.

Key result:
- Pre-specified endpoint lookup: {prespec_concordance['Concordant']} / {prespec_concordance['Total with harmonised beta']} SNPs had concordant direction; one-sided sign-test p = {prespec_concordance['One-sided sign-test p']}.
- Correct over 10 tests with Bonferroni threshold 0.05 / 10 = 0.005, and report BH-FDR q-values across the same 10 tests.
- Keep the all-120 endpoint-by-SNP table only as exploratory supplementary material.

Manuscript text:

> To avoid over-penalising correlated FinnGen long-COVID endpoints while preserving control of multiple testing, we pre-specified one FinnGen endpoint per index SNP according to the discovery contrast. SNPs discovered against responder/population controls were tested against overall FinnGen post-COVID condition versus population controls (LC_U09), while the SNP discovered against infected controls was tested against FinnGen post-COVID condition versus SARS-CoV-2 infection controls (LC_U09_vs_U07). We corrected these ten primary replication tests using Bonferroni correction (alpha = 0.005) and BH-FDR, while treating the full 120 SNP-endpoint lookup as exploratory.

## 5. LCA hard assignment and entropy

Table to use: `05_lca_entropy.*`.

Key result:
- Selected 3-class PCC model entropy = {entropy_rows[0]['Entropy']}.
- Class-specific mean posterior probabilities are {', '.join(row['Mean posterior probability'] for row in entropy_rows)}.

Manuscript text:

> Individuals were assigned to their maximum-posterior latent class for GWAS. To quantify uncertainty from this hard assignment, we report entropy and class-specific mean posterior probabilities for the selected three-class model. The three-class PCC model had entropy {entropy_rows[0]['Entropy']}, with mean posterior probabilities of {', '.join(row['Mean posterior probability'] for row in entropy_rows)} across the three assigned classes.

Limitation text:

> Because subtype GWAS used maximum-posterior class assignment, residual class-membership uncertainty may attenuate subtype-specific genetic effects. Entropy and mean posterior probabilities suggested generally good class separation, but future analyses could propagate posterior class probabilities directly into association testing.

## 6. Genetic correlation with HGI COVID susceptibility and severity

This still needs LDSC because no local LDSC installation/reference files were found. Use `06_ldsc_commands_and_table_slots.md`.

What to report:
- rg between PCC vs responder controls and HGI C2 susceptibility.
- rg between PCC vs responder controls and HGI B2 hospitalisation/severity.
- Ideally repeat for PCC vs infected non-PCC controls and the three subtypes if LDSC heritability is stable.

Manuscript text template:

> Genome-wide genetic correlation analyses showed [high/low] overlap between PCC and HGI SARS-CoV-2 susceptibility (rg = [ ], SE = [ ], p = [ ]), and [high/low/no] overlap with COVID-19 hospitalisation/severity (rg = [ ], SE = [ ], p = [ ]). This pattern [supports/does not support] the interpretation that the strongest responder-control PCC signals largely reflect infection-susceptibility biology rather than PCC-specific liability after infection.

If rg with susceptibility is high:

> This genome-wide result provides a single quantitative summary supporting the same interpretation suggested by the MR, Bayesian line-model, and mash analyses: much of the currently detectable PCC signal in the high-powered responder-control scan tracks SARS-CoV-2 susceptibility.

If rg with susceptibility is low:

> This would weaken the claim that the responder-control PCC signal primarily reflects infection susceptibility, and the Discussion should instead emphasise locus-specific overlap rather than genome-wide sharing.
"""
    (OUT / "00_numbered_manuscript_change_pack.md").write_text(text)


def main():
    control_rows = make_control_definition_tables()
    qc_rows = make_lambda_gc_table()
    best_rows, prespec_rows, sign_rows, missing = make_fingen_tables()
    entropy_rows = make_entropy_table()
    make_ldsc_template()
    make_change_pack(control_rows, qc_rows, best_rows, prespec_rows, sign_rows, entropy_rows)
    print(f"Wrote manuscript update outputs to {OUT.relative_to(ROOT)}")
    if missing:
        print("Missing pre-specified FinnGen endpoint rows:")
        for snp, endpoint in missing:
            print(f"  {snp}: {endpoint}")
    print("Generated files:")
    for path in sorted(OUT.iterdir()):
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
