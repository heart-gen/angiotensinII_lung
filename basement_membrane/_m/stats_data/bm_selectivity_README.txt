Basement-membrane selectivity -- generated summary
Units (>=5 cells): 2329; cell types: 22; donors: 220
Cohort: Healthy donors only

Variance components, primary endpoint (bm_minus_fibrillar):
       grp        var1 var2       vcov      sdcor
1 donor_id (Intercept) <NA> 0.00904089 0.09508359
2    study (Intercept) <NA> 0.01692599 0.13009993
3 Residual        <NA> <NA> 0.05992558 0.24479701

Pre-specified method sanity check:
  LAMA3/LAMA5 NOT pericyte-top3 (expected TRUE): TRUE
  fraction of mural-expected genes in pericyte top-5: 0.80
  COL11A1 NOT pericyte-top5 (expected TRUE): TRUE

Per-gene tau / pericyte rank (BM panel):
       gene       tau            top_group pericyte_rank
     <char>     <num>               <fctr>         <int>
 1:  COL4A1 0.8719076            Pericytes             1
 2:  COL4A2 0.8423271            Pericytes             1
 3: COL18A1 0.8384303            Pericytes             1
 4:   LAMB1 0.8778667            Pericytes             1
 5:    NID1 0.8915455            Pericytes             1
 6:    NID2 0.8030317            Pericytes             1
 7:   LAMA4 0.8016409 Alveolar fibroblasts             3
 8:    AGRN 0.8806290                  AT1             5
 9:   LAMC1 0.7675729 Alveolar fibroblasts             8
10:   LAMA5 0.8138139                  AT1             9
11:   HSPG2 0.7940798  EC venous pulmonary            10
12:   LAMB2 0.5510528                  AT1            13
13:   LAMA3 0.9645018                  AT1            17
    log2_pericyte_over_next
                      <num>
 1:               1.0927463
 2:               0.9921461
 3:               0.8466301
 4:               0.5433248
 5:               1.2165636
 6:               0.3368719
 7:              -0.9314474
 8:              -2.3147346
 9:              -1.5036006
10:              -2.7781505
11:              -2.1388791
12:              -1.2134886
13:              -6.2235693

Per-gene tau / pericyte rank (fibrillar core + minor):
      gene           block       tau                 top_group pericyte_rank
    <char>          <char>     <num>                    <fctr>         <int>
1:  COL5A3 fibrillar_minor 0.9505790                 Pericytes             1
2:  COL5A2 fibrillar_minor 0.8321868 Peribronchial fibroblasts             3
3:  COL1A2  fibrillar_core 0.8281953 Peribronchial fibroblasts             6
4:  COL5A1 fibrillar_minor 0.8384635 Peribronchial fibroblasts             6
5:  COL1A1  fibrillar_core 0.8877511    Subpleural fibroblasts             7
6:  COL3A1  fibrillar_core 0.8451812 Peribronchial fibroblasts             7
7: COL11A2 fibrillar_minor 0.8548388            Myofibroblasts            12
8: COL11A1 fibrillar_minor 0.8032993            Myofibroblasts            22
   log2_pericyte_over_next
                     <num>
1:               1.8823477
2:              -0.4956259
3:              -1.3013682
4:              -1.9165633
5:              -5.0766635
6:              -3.0652024
7:              -3.7741463
8:                    -Inf
