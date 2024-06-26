## Running BayesPrism with multi-cores

library(argparse)
library(BayesPrism)

#### Main
parser <- ArgumentParser()
parser$add_argument('--select', type="integer", default=1,
                    help="the Prism object to select [default %(default)s]")
args <- parser$parse_args()
intx <- args$select[1]
num  <- stringr::str_pad(intx, 2, pad="0")

                                        # Load Prism objects
load("lung_prop_celltypes.RDat")
prism_lt <- list(lungPrism0, lungPrism, lungPrismX)
outfile  <- paste("prism_results", num,
                  "lung_GTEx.RDat", sep=".")
                                        # Run Prism
bp.res <- run.prism(prism=prism_lt[[intx]], n.cores=20)
save(bp.res, file=outfile)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
