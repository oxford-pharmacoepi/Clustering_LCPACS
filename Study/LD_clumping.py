import pandas as pd
import subprocess
from pathlib import Path
import os

# ---- paths ----
# Assuming script is run from the 'Study' folder or root, we try to resolve paths relative to the repo root
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT  = SCRIPT_DIR.parent

# Input files (HGI Sumstats)
HGI_C2_PATH = REPO_ROOT / "Bayes/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz"
HGI_B2_PATH = REPO_ROOT / "Bayes/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz"

# Reference data
# NOTE: You must unzip '1000G_EUR_Phase3_plink.zip' into this folder in the root
LDREF_DIR = REPO_ROOT / "1000G_EUR_Phase3_plink"

# Output directory
OUT_DIR = SCRIPT_DIR / "clumped_results"
OUT_DIR.mkdir(exist_ok=True)

# PLINK executable
# We look for 'plink' in the repo root first, then fall back to system PATH
PLINK_BIN_PATH = REPO_ROOT / "plink"
if PLINK_BIN_PATH.exists():
    PLINK_BIN = str(PLINK_BIN_PATH)
else:
    PLINK_BIN = "plink"

def run_clumping(sumstats_file, output_tag):
    print(f"\n--- Processing {output_tag} ---")
    if not sumstats_file.exists():
        print(f"Error: Input file not found: {sumstats_file}")
        return

    # 1. Prepare a PLINK-style sumstats file: SNP, CHR, BP, P
    print(f"Reading {sumstats_file.name}...")
    df = pd.read_csv(sumstats_file, sep="\t")

    # HGI columns: #CHR, POS, rsid, all_inv_var_meta_p
    plink_df = pd.DataFrame({
        "SNP": df["rsid"],
        "CHR": df["#CHR"],
        "BP":  df["POS"],
        "P":   df["all_inv_var_meta_p"]
    })

    # Drop variants without rsid or p-value
    plink_df = plink_df.dropna(subset=["SNP", "P"])

    temp_sumstats = OUT_DIR / f"{output_tag}_sumstats_for_plink.txt"
    plink_df.to_csv(temp_sumstats, sep="\t", index=False)
    print(f"Created temporary sumstats: {temp_sumstats}")

    # 2. Run PLINK clumping per chromosome
    clumped_snps = []

    if not LDREF_DIR.exists():
        print(f"CRITICAL ERROR: Reference directory {LDREF_DIR} does not exist.")
        print("Please unzip '1000G_EUR_Phase3_plink.zip' into the root of the workspace.")
        return

    print("Running PLINK clumping per chromosome...")
    for chrom in range(1, 23):
        # Expected reference file pattern: 1000G.EUR.QC.{chrom}.bed/.bim/.fam
        # We construct the prefix for PLINK
        bfile_prefix = LDREF_DIR / f"1000G.EUR.QC.{chrom}"

        # Check if reference file exists.
        # Note: We construct the full path to the .bed file explicitly to avoid pathlib replacing .1, .2 etc as extensions
        bed_file = LDREF_DIR / f"1000G.EUR.QC.{chrom}.bed"

        if not bed_file.exists():
            # Fallback check for different naming if needed, or just warn
            print(f"  Warning: Reference file for chr{chrom} not found at {bed_file}")
            continue

        out_prefix = OUT_DIR / f"{output_tag}_chr{chrom}"

        # Parameters matching Lammi et al. (Nature Genetics):
        # P < 5e-8, r2 < 0.001, kb = 10000 (10Mb)
        cmd = [
            PLINK_BIN,
            "--bfile", str(bfile_prefix),
            "--clump", str(temp_sumstats),
            "--clump-p1", "5e-8",
            "--clump-p2", "1e-6",
            "--clump-r2", "0.001",
            "--clump-kb", "10000",
            "--chr", str(chrom),
            "--out", str(out_prefix),
            "--memory", "4000" # Limit memory usage if needed
        ]

        # Run PLINK silently
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check for output file
        clumped_file = out_prefix.with_suffix(".clumped")
        if clumped_file.exists():
            try:
                c_df = pd.read_csv(clumped_file, delim_whitespace=True)
                if not c_df.empty:
                    snps = c_df["SNP"].tolist()
                    clumped_snps.extend(snps)
                    # print(f"  Chr{chrom}: {len(snps)} SNPs")
            except Exception as e:
                print(f"  Error reading results for chr{chrom}: {e}")
        else:
            # If PLINK failed or found no variants, it might not write .clumped
            # We can check result.stderr if needed, but usually it just means no significant hits
            pass

    # 3. Collect and save
    clumped_snps = sorted(list(set(clumped_snps)))
    print(f"Total independent instruments for {output_tag}: {len(clumped_snps)}")

    out_file = OUT_DIR / f"{output_tag}_clumped_snps.txt"
    pd.Series(clumped_snps, name="SNP").to_csv(out_file, index=False)
    print(f"Saved list to: {out_file}")

if __name__ == "__main__":
    print("Starting LD Clumping Pipeline...")
    run_clumping(HGI_C2_PATH, "c2")
    run_clumping(HGI_B2_PATH, "b2")
    print("\nDone. Now run the R script.")
