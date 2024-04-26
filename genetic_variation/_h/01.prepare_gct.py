"""
This script is to prepare GCT files based on GTEx pipeline.
"""
import pandas as pd
from pyhere import here
from functools import lru_cache

def to_gct(filename, df):
    description_df = pd.DataFrame({'Description': df.index.values},
                                  index=df.index)
    dfo = pd.concat([description_df, df], axis=1)
    dfo.index.name = 'Names'
    with open(filename, "wt") as out:
        print("#1.2", file=out)
        print(df.shape[0], df.shape[1], sep="\t", file=out)
        dfo.to_csv(out, sep="\t")


@lru_cache()
def get_pheno():
    return pd.read_csv("gtex_phenotypes_lung.tsv", sep='\t')


@lru_cache()
def get_fam():
    fam_file = here("inputs/gtex/genotypes/_m/gtex_v8.fam")
    return pd.read_csv(fam_file, sep="\t", header=None,
                       names=["subjid","ID","V2","V3","V4","V5"])


@lru_cache()
def get_celltypes():
    return pd.read_csv("gtex_cell_proportions.tsv", sep="\t",
                       index_col=0)


@lru_cache()
def load_data():
    pheno_df = get_pheno()
    pheno_df["ids"] = pheno_df.sampid
    pheno_df.set_index("ids", inplace=True)
    tpm_file = here("inputs/gtex/_m",
                    "GTEx_Analysis_2017-06-05_v8"+\
                    "_RNASeQCv1.1.9_gene_tpm.gct.gz")
    tpm_df = pd.read_csv(tpm_file, index_col=0, skiprows=2, sep='\t')
    counts_file = here("inputs/gtex/_m/genes_gtex_v8_counts.txt.gz")
    counts_df = pd.read_csv(counts_file, sep="\t",
                            index_col=0, skiprows=2)
    samples = list(set(counts_df.columns).intersection(set(pheno_df["sampid"])))
    return pheno_df.loc[samples,:], tpm_df.loc[:,samples], counts_df.loc[:,samples]


def select_idv(pheno_df, counts_df):
    samples = list(set(pheno_df.loc[counts_df.columns,:].subjid)\
                   .intersection(set(get_fam().subjid)))
    new_fam = get_fam()[(get_fam()["subjid"].isin(samples))]\
        .drop_duplicates(subset="subjid")
    new_fam.to_csv("keepFam.txt", sep='\t', index=False, header=False)
    pheno_df[["sampid", "subjid"] + list(get_celltypes().columns.values)]\
        .reset_index().set_index("subjid").loc[new_fam.subjid]\
        .to_csv("celltypes_interaction_list.txt", sep='\t')
    new_pheno = pheno_df.loc[:, ["sampid", "subjid"]]\
                        .reset_index().set_index("subjid")\
                        .loc[new_fam.subjid].reset_index().set_index("ids")
    return new_pheno


def main():
    pheno_df, tpm_df, counts_df = load_data()
    new_pheno = select_idv(pheno_df, counts_df)
    genes = list(set(counts_df.index).intersection(set(tpm_df.index)))
    to_gct("counts.gct", counts_df.loc[genes,new_pheno.index])
    to_gct("tpm.gct", tpm_df.loc[genes,new_pheno.index])
    new_pheno.loc[:, ["sampid", "subjid"]]\
             .to_csv("sample_id_to_subject_id.tsv", sep="\t", index=False)
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('vcf_chr_list.txt', header=False, index=None)


if __name__=='__main__':
    main()
