
"""
This script runs tensorQTL in python.
"""
import session_info
import pandas as pd
from functools import lru_cache
from tensorqtl import susie, read_phenotype_bed, genotypeio

@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(), variant_df


@lru_cache()
def get_eGenes():
    df = pd.read_csv("../_m/GTEx_v8.genes.txt.gz", sep="\t")\
           .loc[:, ["phenotype_id", "qval", "pval_nominal_threshold"]]
    return df[(df["qval"] <= 0.05)].copy()


@lru_cache()
def get_eqtl():
    df0 = pd.read_csv("../_m/GTEx_v8.nominal.txt.gz", sep="\t", nrows=500,
                      usecols=[0,1,6,7])
    df1 = pd.read_csv("../_m/GTEx_v8.nominal.txt.gz", sep="\t",
                      usecols=[0,1,6,7], dtype=df0.dtypes.to_dict(),
                      engine="c", compression="gzip")
    df1 = df1[(df1["phenotype_id"].isin(get_eGenes().phenotype_id))].copy()
    df = df1.merge(get_eGenes(), on="phenotype_id")
    return df[(df["pval_nominal"] <= df["pval_nominal_threshold"])]\
        .drop(["qval", "pval_nominal_threshold"], axis=1)


@lru_cache()
def get_covars(feature="genes"):
    covar_file = f"../../covariates/_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype(feature="genes"):
    expr_bed = f"../../normalize_expression/_m/{feature}.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


def main():
    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype()
    genotype_df, variant_df = get_genotypes()
    covariates_df = get_covars()
    eqtl_df = get_eqtl()
    prefix = "GTEx_v8"

    # Filter for cis-eQTL (eGenes)
    phenotype_df = phenotype_df.loc[get_eGenes().phenotype_id,:]
    phenotype_pos_df = phenotype_pos_df.loc[get_eGenes().phenotype_id,:]
    genotype_df = genotype_df.loc[eqtl_df.variant_id.unique(), :]
    variant_df = variant_df.loc[eqtl_df.variant_id.unique(), :]
    
    # Run SuSiE fine mapping by chromosome
    chroms = sorted(phenotype_pos_df.chr.unique())
    susie_df = []
    for k,chrom in enumerate(chroms, 1):
        print(f'  * processing chr. {k}/{len(chroms)}', flush=True)
        dfx = susie.map(genotype_df, variant_df,
                        phenotype_df.loc[phenotype_pos_df['chr']==chrom],
                        phenotype_pos_df.loc[phenotype_pos_df['chr']==chrom],
                        covariates_df=covariates_df, max_iter=200,
                        maf_threshold=0.05, window=500000)
        susie_df.append(dfx)
    pd.concat(susie_df, axis=0)\
      .to_csv(f"{prefix}.susie.txt.gz", sep='\t', index=False)


if __name__ == "__main__":
    main()
