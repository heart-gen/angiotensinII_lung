"""
This script runs tensorQTL in python.
"""
import session_info
import pandas as pd
from functools import lru_cache
from tensorqtl import genotypeio, cis
from tensorqtl import read_phenotype_bed, calculate_qvalues

@lru_cache()
def get_covars():
    covar_file = "../../covariates/_m/genes.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0).T


@lru_cache()
def get_phenotype():
    expr_bed = "../../normalize_expression/_m/genes.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(), variant_df


def main():
    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype()
    genotype_df, variant_df = get_genotypes()
    prefix = "GTEx_v8"

    # Nominal
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                    prefix, covariates_df=get_covars(), maf_threshold=0.05,
                    window=1000000, output_dir='.', verbose=True)

    # Permutation
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df,phenotype_pos_df,
                         covariates_df=get_covars(), maf_threshold=0.05,
                         window=1000000, seed=220705, verbose=True)
    calculate_qvalues(cis_df, fdr=0.05)
    cis_df.to_csv(f"{prefix}.genes.txt.gz", sep='\t')
    
    # Conditional
    indep_df = cis.map_independent(genotype_df, variant_df, cis_df,phenotype_df,
                                   phenotype_pos_df, covariates_df=get_covars(),
                                   maf_threshold=0.05, fdr=0.05, fdr_col='qval',
                                   window=1000000, seed=220705, verbose=True)
    indep_df.to_csv(f"{prefix}.conditional.txt.gz", sep='\t')

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
