Basement-membrane selectivity -- generated summary
Units (>=5 cells): 2329; cell types: 22; donors: 220
Cohort: Healthy donors only

Variance components, primary endpoint (bm_minus_fibrillar):
       grp        var1 var2       vcov     sdcor
1 donor_id (Intercept) <NA> 0.01044342 0.1021931
2    study (Intercept) <NA> 0.01688290 0.1299342
3 Residual        <NA> <NA> 0.05424585 0.2329074

Pre-specified method sanity check:
  LAMA3/LAMA5 NOT pericyte-top3 (expected TRUE): TRUE
  fraction of mural-expected genes in pericyte top-5: 0.80
  COL11A1 NOT pericyte-top5 (expected TRUE): TRUE

Per-gene tau / pericyte rank (BM panel):
       gene       tau                 top_group pericyte_rank
     <char>     <num>                    <fctr>         <int>
 1:  COL4A1 0.8719076                 Pericytes             1
 2:  COL4A2 0.8423271                 Pericytes             1
 3: COL18A1 0.8384303                 Pericytes             1
 4:   LAMB1 0.8778667                 Pericytes             1
 5:   LAMB4 0.8648643                 Pericytes             1
 6:   LAMC3 0.9935237                 Pericytes             1
 7:    NID1 0.8915455                 Pericytes             1
 8:    NID2 0.8030317                 Pericytes             1
 9:   LAMA4 0.8016409      Alveolar fibroblasts             3
10:   LAMA2 0.8300583   Adventitial fibroblasts             5
11:    AGRN 0.8806290                       AT1             5
12:   LAMA1 0.6388169      Alveolar macrophages             7
13:   LAMB3 0.9306873                       AT1             8
14:   LAMC1 0.7675729      Alveolar fibroblasts             8
15:   LAMA5 0.8138139                       AT1             9
16: COL15A1 0.8811074 Peribronchial fibroblasts            10
17:   HSPG2 0.7940798       EC venous pulmonary            10
18:   LAMB2 0.5510528                       AT1            13
19:   LAMA3 0.9645018                       AT1            17
20:   LAMC2 0.9567766                       AT1            20
       gene       tau                 top_group pericyte_rank
    log2_pericyte_over_next
                      <num>
 1:               1.0927463
 2:               0.9921461
 3:               0.8466301
 4:               0.5433248
 5:               1.0674290
 6:               4.7465630
 7:               1.2165636
 8:               0.3368719
 9:              -0.9314474
10:              -0.6425932
11:              -2.3147346
12:              -0.9875248
13:              -4.6848095
14:              -1.5036006
15:              -2.7781505
16:              -4.8006950
17:              -2.1388791
18:              -1.2134886
19:              -6.2235693
20:              -5.9447401
    log2_pericyte_over_next

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
