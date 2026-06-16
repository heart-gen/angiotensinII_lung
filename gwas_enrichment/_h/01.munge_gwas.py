#!/usr/bin/env python
"""
01.munge_gwas.py — harmonize staged GWAS to rsID-keyed inputs for MAGMA and LDSC.

Each trait is parsed by its native format and emitted as:
  magma/<trait>.pval     columns: SNP  P  N      (MAGMA --pval ... use=SNP,P ncol=N)
  ldsc_raw/<trait>.tsv   columns: SNP  A1  A2  N  Z   (input to LDSC munge_sumstats.py)

Everything is keyed on **rsID** so the pipeline can stay in hg19 (g1000_eur + 1000G Phase3
panels) regardless of each GWAS's native build. Traits lacking a signed effect (COPDGene
reports only |beta|) get a MAGMA file but no LDSC file (LDSC needs signed Z).

Run one or more traits with --traits (default: all in the registry).
"""
import os
import argparse
import numpy as np
import pandas as pd
import session_info

IN = "../../inputs/gwas/_m"

# trait -> dict(file, fmt, [n])  ; fmt selects the parser below
REG = {
    "COPD_COPDGene": dict(file="COPD_COPDGene.phs000179.pha004496.txt", fmt="copdgene"),
    "COPD_UKB":      dict(file="COPD_UKB.imputed_v3.both_sexes.tsv.bgz", fmt="neale"),
    "IPF_UKB":       dict(file="IPF_UKB.imputed_v3.both_sexes.tsv.bgz",  fmt="neale"),
    "ASTHMA":        dict(file="ASTHMA_TAGC_EUR.txt.gz",                 fmt="imputed"),
    "SMKINIT":       dict(file="SMKINIT_GSCAN2_EUR.txt",                 fmt="gscan"),
    "CIGDAY":        dict(file="CIGDAY_GSCAN2_EUR.txt",                  fmt="gscan"),
    "SBP":           dict(file="SBP_ICBP.txt.gz",                        fmt="imputed"),
    "DBP":           dict(file="DBP_ICBP.txt.gz",                        fmt="imputed"),
    "CAD":           dict(file="CAD_CARDIOGRAM.txt.gz",                  fmt="imputed"),
    "EDU":           dict(file="EDU_SSGAC.txt.gz",                       fmt="imputed"),
    "BODYFAT":       dict(file="BODYFAT_UKB.txt.gz",                     fmt="imputed"),
    # External lung-function (present only after the Shrine download in 00.stage_gwas.sh);
    # GWAS-Catalog harmonised build37 files (PMID 30804560, N=400,102 EUR, signed beta).
    "FEV1":          dict(file="FEV1_Shrine2019.build37.tsv.gz",        fmt="shrine", n=400102),
    "FVC":           dict(file="FVC_Shrine2019.build37.tsv.gz",         fmt="shrine", n=400102),
    "FEV1FVC":       dict(file="FEV1FVC_Shrine2019.build37.tsv.gz",     fmt="shrine", n=400102),
}

_VARIANTS_CACHE = {}


def load_neale_variants():
    """UKB Neale variant annotation: variant -> rsid, ref, alt (GRCh37)."""
    if "v" in _VARIANTS_CACHE:
        return _VARIANTS_CACHE["v"]
    path = os.path.join(IN, "UKB_variants.tsv.bgz")
    v = pd.read_csv(path, sep="\t", compression="gzip",
                    usecols=["variant", "rsid", "ref", "alt"], dtype=str)
    v = v.dropna(subset=["rsid"])
    v = v[v["rsid"].str.startswith("rs")]
    _VARIANTS_CACHE["v"] = v
    return v


def parse_imputed(path):
    """imputed_gwas_hg38: variant_id (rsid), effect/non_effect_allele, zscore, pvalue, sample_size."""
    df = pd.read_csv(path, sep="\t", compression="infer",
                     usecols=lambda c: c in {"variant_id", "effect_allele", "non_effect_allele",
                                             "zscore", "pvalue", "sample_size"})
    df = df.rename(columns={"variant_id": "SNP", "effect_allele": "A1",
                            "non_effect_allele": "A2", "zscore": "Z",
                            "pvalue": "P", "sample_size": "N"})
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    return df[["SNP", "A1", "A2", "Z", "P", "N"]]


def parse_shrine(path, n=400102):
    """GWAS-Catalog harmonised build37 (.f.tsv.gz): hm_rsid + harmonised alleles/beta.

    Columns vary slightly across releases, so resolve rsID / allele / effect / p columns
    by falling back from the harmonised (hm_*) names to the raw names. Shrine 2019 has no
    per-SNP N, so the published discovery N (400,102 EUR) is assigned (documented fallback).
    """
    head = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    cols = set(head.columns)

    def pick(*cands):
        for c in cands:
            if c in cols:
                return c
        return None

    c_snp = pick("hm_rsid", "rsid", "variant_id")
    c_a1  = pick("hm_effect_allele", "effect_allele")
    c_a2  = pick("hm_other_allele", "other_allele")
    c_b   = pick("hm_beta", "beta")
    c_se  = pick("standard_error", "se")
    c_p   = pick("p_value", "pvalue", "p")
    c_n   = pick("n", "sample_size", "N")
    use = [c for c in {c_snp, c_a1, c_a2, c_b, c_se, c_p, c_n} if c]
    df = pd.read_csv(path, sep="\t", compression="infer", usecols=use)
    df = df.rename(columns={c_snp: "SNP", c_a1: "A1", c_a2: "A2",
                            c_b: "beta", c_se: "se", c_p: "P"})
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    df["Z"] = pd.to_numeric(df["beta"], errors="coerce") / pd.to_numeric(df["se"], errors="coerce")
    df["N"] = pd.to_numeric(df[c_n], errors="coerce") if c_n else n   # per-SNP N if present
    return df[["SNP", "A1", "A2", "Z", "P", "N"]]


def parse_gscan(path):
    """GSCAN2: CHR POS RSID EFFECT_ALLELE OTHER_ALLELE AF_1000G BETA SE P N (whitespace)."""
    df = pd.read_csv(path, sep=r"\s+",
                     usecols=["RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P", "N"])
    df = df.rename(columns={"RSID": "SNP", "EFFECT_ALLELE": "A1", "OTHER_ALLELE": "A2"})
    df["Z"] = df["BETA"] / df["SE"]
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    return df[["SNP", "A1", "A2", "Z", "P", "N"]]


def parse_neale(path):
    """UKB Neale v3: variant beta se pval n_complete_samples; join variants for rsid/alleles."""
    df = pd.read_csv(path, sep="\t", compression="gzip",
                     usecols=["variant", "beta", "se", "pval", "n_complete_samples",
                              "low_confidence_variant"])
    df = df[df["low_confidence_variant"].astype(str).str.lower() != "true"]
    v = load_neale_variants()
    df = df.merge(v, on="variant", how="inner")
    df = df.rename(columns={"rsid": "SNP", "alt": "A1", "ref": "A2",
                            "pval": "P", "n_complete_samples": "N"})
    df["Z"] = df["beta"] / df["se"]
    return df[["SNP", "A1", "A2", "Z", "P", "N"]]


def parse_copdgene(path):
    """COPDGene dbGaP: '#'-commented preamble then header; |beta| is unsigned -> MAGMA only."""
    hdr = None
    with open(path) as fh:
        for i, line in enumerate(fh):
            if line.startswith("ID\t") or "SNP ID" in line.split("\t")[:3]:
                hdr = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=hdr)
    df = df.rename(columns={"SNP ID": "SNP", "P-value": "P", "Sample size": "N",
                            "Coded Allele": "A1", "Allele1": "_a1", "Allele2": "_a2"})
    df = df[df["SNP"].astype(str).str.startswith("rs")]
    # N column is sparsely populated in this dbGaP export; fall back to published N.
    n = pd.to_numeric(df.get("N"), errors="coerce")
    if n.notna().sum() == 0 or np.nanmedian(n) < 1000:
        df["N"] = 20066  # COPDGene COPD-Primary meta-analysis (documented fallback)
    return df[["SNP", "P", "N"]], None   # no signed Z available


def write_outputs(trait, parsed, outdir):
    magma_dir = os.path.join(outdir, "magma"); os.makedirs(magma_dir, exist_ok=True)
    ldsc_dir = os.path.join(outdir, "ldsc_raw"); os.makedirs(ldsc_dir, exist_ok=True)
    if isinstance(parsed, tuple):                 # MAGMA-only (no Z)
        mp, _ = parsed
        mp = mp.dropna(subset=["SNP", "P"]).drop_duplicates("SNP")
        mp.to_csv(os.path.join(magma_dir, f"{trait}.pval"), sep="\t", index=False)
        print(f"  {trait}: MAGMA {mp.shape[0]} SNPs (no LDSC — unsigned effect)")
        return
    df = parsed.dropna(subset=["SNP", "P"]).drop_duplicates("SNP")
    df["A1"] = df["A1"].astype(str).str.upper()
    df["A2"] = df["A2"].astype(str).str.upper()
    df[["SNP", "P", "N"]].to_csv(os.path.join(magma_dir, f"{trait}.pval"), sep="\t", index=False)
    df[["SNP", "A1", "A2", "N", "Z"]].dropna().to_csv(
        os.path.join(ldsc_dir, f"{trait}.tsv"), sep="\t", index=False)
    print(f"  {trait}: MAGMA + LDSC {df.shape[0]} SNPs")


PARSERS = {"imputed": parse_imputed, "gscan": parse_gscan, "shrine": parse_shrine,
           "neale": parse_neale, "copdgene": parse_copdgene}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traits", nargs="*", default=list(REG))
    ap.add_argument("--outdir", default="../_m")
    args = ap.parse_args()
    for trait in args.traits:
        if trait not in REG:
            print(f"  skip unknown trait {trait}"); continue
        path = os.path.join(IN, REG[trait]["file"])
        if not os.path.exists(path):
            print(f"  skip {trait}: file not staged ({REG[trait]['file']})"); continue
        print(f"Munging {trait} [{REG[trait]['fmt']}] ...", flush=True)
        parsed = PARSERS[REG[trait]["fmt"]](path)
        write_outputs(trait, parsed, args.outdir)
    session_info.show()


if __name__ == "__main__":
    main()
