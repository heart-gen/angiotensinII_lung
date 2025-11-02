***CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF.txt.gz  —> is a single file for all GTEx tissues and all eGenes where we reported
the CPP (Causal Posterior Probability). Sample header file:

TISSUE	GENE	eQTL	CHROM	POS	Probability
Brain_Caudate_basal_ganglia	ENSG00000248485.1	1_161274374	1	161274374	0.157456
Brain_Caudate_basal_ganglia	ENSG00000248485.1	1_161687648	1	161687648	0.211998
Brain_Caudate_basal_ganglia	ENSG00000248485.1	1_161319639	1	161319639	0.119975

***CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt.gz  —> is a single file for all GTEx tissues and all eGenes where we reported
the CPP (Causal Posterior Probability). Each eQTL contains the allele information. Sample header file:

TISSUE	GENE	eQTL	CHROM	POS	Probability
Brain_Caudate_basal_ganglia	ENSG00000248485.1	chr1_161274374_G_A_b38	1	161274374	0.157456


***CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz --> is a single file for all GTEx tissues and all eGene where we report
all the high causal variants (variants that have posterior probability of > 0.1).


***CAVIAR_Results_v8_GTEx_LD_ALL.tar.gz --> contains all the CAVIAR output.
Each tissue has a folder and the result is divided by chromosome. We have
3 files for each tissue/chromosome:
 
1) Set --> 95%-Causal Set
This file contains the set of SNPs that capture all the causal variants
with probability at least 0.95. 
2) Post --> causal posterior probability (CPP)
This file contains the causal posterior probability (3rd column) for set of SNPs for a given
eGene in a tissue
3) Hist  --> Allelic heterogeneity
This file contains the probability of having i-th independent causal variants.
One line which the first number indicates the probability of having no-causal
variants, the second number indicates the probability of having one causal    
variant, the third number indicates the probability of having two causal variants.


***GTEx-v8-Annot-Enr-GTEx-LD.out --contains the 10 top enrichment for fine-mapped eQTL
with RoadMap/ENCODE annotations.


The following command was run to produce the data:
"${CAVIAR_PATH}/CAVIAR " + \
                         " -o " + ${OUTPUT_FILE} + \
                         " -z " + ${Z_FILE} + \
                         " -l " + ${LD_FILE} + \
                         " -g 0.001 " +\
                         " -c 5 " + \
                         " -f 1 -r 0.95";
