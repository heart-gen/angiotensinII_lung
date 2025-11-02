## This script examines AGTR2 expression in AT2b cells

import session_info
import polars as pl

def load_data():
    fn = "../_m/at2_phate.normalized_expression.tsv.gz"
    return pl.read_csv(fn, separator="\t")\
             .with_columns(pl.when(pl.col("AGTR2") > 0)\
                           .then(1).otherwise(0)\
                           .alias("present"))


def main():
    ## Review within AT2 cells
    norm_df = load_data()
    print(norm_df.shape)
    print(norm_df.group_by("PHATE")\
          .agg([pl.count().alias("N"),
	        pl.sum("present").alias("AGTR2_positive")])\
          .with_columns((pl.col("AGTR2_positive")/pl.col("N"))\
                        .alias("percent_AGTR2")))

    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()



