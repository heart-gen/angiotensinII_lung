## SNP vcf annotation file format

The GTEx_v8_finemapping_DAPG.vcf.gz file follows the standard VCF format with a single annotation field. Each SNP annotation field contains a set of records separated by `|`, e.g.,

```
ENSG00000269981:1@Spleen=2.00812e-01[9.997e-01:4]|ENSG00000238009:3@Testis=9.68696e-02[9.998e-01:4]|ENSG00000268903:1@Spleen=8.07834e-02[9.919e-01:3]|
```

These records provide the complete eQTL information on the target SNP with respect to all genes and all tissues.  Each record has the following format

```
signal_id@tissue_name= PIP[SPIP:size_of_cluster]
```

More specifically,

 + signal_id: ID of a signal cluster. It consists of the name of the gene and the signal index separated by ":". e.g., ENSG00000238009:3 indicates the signal is the third signal from gene ENSG00000238009
 + tissue_name: name of the tissue where the SNP is investigated
 + PIP: SNP posterior inclusion probability. Higer PIP value indicates the SNP is more likely to be the casual eQTL.
 + SPIP: signal-level posterior_inclusion probability (sum of the PIPs from all members of the signal cluster)
 + size_of_cluster: number of SNPs in the signal cluster. These member SNPs are in LD, all represent the same underlying association signal




## Credible set information file

The GTEx_v8_finemapping_DAPG.CS95.txt.gz file is also in VCF format, but only includes signal clusters with SPIP > 0.95. The SNPs listed in the file are members of 95% credible sets of strong eQTL signals.

The annotation field contains a set of records separated by `|`, e.g.,

```
ENSG00000117620.14:1@Brain_Cortex=3.810e-02|ENSG00000117620.14:1@Lung=6.536e-03
```

Each record has the following format
```
signal_id@tissue_name=PIP
```

More specifically,

+ signal_id:  ID of a signal cluster. It consists of the name of the gene and the signal index separated by ":". Note that the signal clusters in this file all have SPIP > 0.95
+ tissue_name: name of the tissue where the eQTL signal is evaluated
+ PIP: SNP posterior inclusion probability. Higer PIP value indicates the SNP is more likely to be the casual eQTL.
