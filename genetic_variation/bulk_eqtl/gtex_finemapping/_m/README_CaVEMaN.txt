The procedures for running CaVEMaN can be found here: https://funpopgen.github.io/CaVEMaN/

For the GTEx analysis we used the same genotype, covariate and expression files as the eQTL mapping.
Firstly, the --single-signal flag and the set of conditional eQTL were passed to produce a set of
expression phenotypes which only capture one eQTL signal. Then CaVEMaN was run with the standard
parameters and 10,000 bootstrap samples, before using the --best flag on the results output, to
extract the posterior probabilities for the lead expression variants.
