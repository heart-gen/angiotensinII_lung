"""
Build the GTEx phenotype table used by the lung age-correlation analysis.

Inputs (copied into inputs/gtex/_m so this repo is independent of the source):
  - GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt  (public v11 sample attrs:
        SAMPID, SMTS, SMTSD, SMRIN, SMTSISCH)
  - phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz  (protected
        subject phenotypes with *continuous* AGE, SEX, RACE, DTHHRDY keyed by SUBJID;
        SUBJID is stable across GTEx versions)

SUBJID is derived from SAMPID (first two dash-delimited fields) to join the two tables.
"""
import os
import pandas as pd

MDIR = os.path.join(os.path.dirname(__file__), "..", "_m")

SAMPLE_ATTR = os.path.join(MDIR, "GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")
SUBJECT     = os.path.join(
    MDIR, "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz")

SAMPLE_COLS  = ["SAMPID", "SMTS", "SMTSD", "SMRIN", "SMTSISCH"]
SUBJECT_COLS = ["SUBJID", "SEX", "AGE", "RACE", "DTHHRDY"]


def subjid_from_sampid(sampid):
    return "-".join(str(sampid).split("-")[:2])


def main():
    attr = pd.read_csv(SAMPLE_ATTR, sep="\t", low_memory=False)[SAMPLE_COLS]
    attr["SUBJID"] = attr["SAMPID"].map(subjid_from_sampid)

    # The protected subject file has leading comment (#) lines and a blank line
    # before the header; comment + skip_blank_lines handle both.
    subj = pd.read_csv(SUBJECT, sep="\t", comment="#",
                       skip_blank_lines=True)[SUBJECT_COLS]

    merged = attr.merge(subj, on="SUBJID", how="inner")
    merged.to_csv(os.path.join(MDIR, "gtex_v11_sample_data.tsv"),
                  sep="\t", index=False)
    merged[merged["SMTS"] == "Lung"].to_csv(
        os.path.join(MDIR, "gtex_v11_lung_sample_data.tsv"), sep="\t", index=False)
    print("wrote gtex_v11_sample_data.tsv: %d samples (%d lung)"
          % (len(merged), int((merged["SMTS"] == "Lung").sum())))


if __name__ == "__main__":
    main()
