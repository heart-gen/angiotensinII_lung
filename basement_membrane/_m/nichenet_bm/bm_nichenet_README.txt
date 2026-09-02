BM-target ligand ranking vs the frozen ECM-program ranking
Spearman rank correlation: 0.227
Top-20 overlap: 8/20
TGF-beta ligands in the BM ranking:
# A tibble: 3 x 5
  test_ligand  rank aupr_corrected perm_z perm_p_BH
  <chr>       <int>          <dbl>  <dbl>     <dbl>
1 TGFB2           5        0.0152  2.94       0.235
2 TGFB1           8        0.0109  2.07       0.235
3 TGFB3         266        0.00123 0.0312     0.269

Top 15 ligands toward basement membrane targets:
# A tibble: 15 x 5
   test_ligand  rank aupr_corrected perm_z perm_p_BH
   <chr>       <int>          <dbl>  <dbl>     <dbl>
 1 MMP14           1        0.137    24.6     0.0321
 2 TIMP2           2        0.0709   13.4     0.235 
 3 MXRA5           3        0.0235    4.22    0.235 
 4 COL5A1          4        0.0166    3.00    0.235 
 5 TGFB2           5        0.0152    2.94    0.235 
 6 COL1A1          6        0.0141    2.66    0.235 
 7 COPA            7        0.0135    2.41    0.235 
 8 TGFB1           8        0.0109    2.07    0.235 
 9 FGF7            9        0.0102    1.69    0.235 
10 CALM3          10        0.0101    1.56    0.235 
11 LAMA2          11        0.00952   1.59    0.235 
12 POSTN          12        0.00950   1.55    0.235 
13 RELN           13        0.00950   1.77    0.235 
14 FGF2           14        0.00941   1.70    0.235 
15 VEGFC          15        0.00868   1.47    0.235 
