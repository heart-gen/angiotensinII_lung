"""
This script generates GTEx v8 phenotype data.
"""
import pandas as pd
from functools import lru_cache

@lru_cache()
def get_samples():
    samples_file = "/ceph/tmp/gtex/download/_m/phs000424.v8.pht002741.v8."+\
        "p2.GTEx_Sample.MULTI.txt.gz"
    return pd.read_csv(samples_file, sep='\t', comment="#")


@lru_cache()
def get_subjects():
    subject_file = "/ceph/tmp/gtex/download/_m/phs000424.v8.pht002742.v8."+\
        "p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"
    return pd.read_csv(subject_file, sep='\t', comment="#")


@lru_cache()
def get_attr():
    attr_file = "/ceph/tmp/gtex/download/_m/phs000424.v8.pht002743.v8."+\
        "p2.c1.GTEx_Sample_Attributes.GRU.txt.gz"
    return pd.read_csv(attr_file, sep='\t', comment="#", low_memory=False)


@lru_cache()
def generate_data():
    return get_samples().merge(get_attr(), on=["dbGaP_Sample_ID","SAMPID"])\
                        .merge(get_subjects(), on=["dbGaP_Subject_ID","SUBJID"])


@lru_cache()
def clean_data():
    return pd.concat([generate_data().dropna(axis=1),
                      generate_data().SMRIN], axis=1)


def main():
    clean_data().to_csv("gtex_v8_sample_data.tsv", sep='\t', index=False)
    clean_data()[(clean_data()["SMTS"] == "Lung")]\
        .to_csv("gtex_v8_lung_sample_data.tsv", sep='\t', index=False)


if __name__ == "__main__":
    main()
