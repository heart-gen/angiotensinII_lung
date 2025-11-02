"""
This script runs tensorQTL in python.
"""

import pandas as pd
import argparse, session_info
from functools import lru_cache
from tensorqtl.post import get_significant_pairs

@lru_cache()
def get_permutation_results(prefix):
    return pd.read_csv(f"{prefix}.genes.txt.gz", sep="\t", index_col=0)


@lru_cache()
def get_eqtl(prefix):
    return get_significant_pairs(get_permutation_results(prefix), prefix)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default="AA")
    args=parser.parse_args()

    # Load and save eQTL results
    get_eqtl(args.prefix)\
        .to_csv(f"{args.prefix}.signif_variants.txt.gz", sep='\t', index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
