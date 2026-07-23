BM-target ligand ranking vs the frozen ECM-program ranking
Spearman rank correlation: 0.223
Top-20 overlap: 8/20
TGF-beta ligands in the BM ranking:
# A tibble: 3 x 5
  test_ligand  rank aupr_corrected perm_z perm_p_BH
  <chr>       <int>          <dbl>  <dbl>     <dbl>
1 TGFB2           5        0.0182  3.20       0.189
2 TGFB1          12        0.0103  1.70       0.189
3 TGFB3         274        0.00124 0.0257     0.241

Top 15 ligands toward BM targets:
# A tibble: 15 x 5
   test_ligand  rank aupr_corrected perm_z perm_p_BH
   <chr>       <int>          <dbl>  <dbl>     <dbl>
 1 MMP14           1        0.164    27.1     0.0321
 2 TIMP2           2        0.0853   14.4     0.189 
 3 MXRA5           3        0.0285    5.22    0.189 
 4 COL5A1          4        0.0197    3.19    0.189 
 5 TGFB2           5        0.0182    3.20    0.189 
 6 COL1A1          6        0.0168    2.97    0.189 
 7 COPA            7        0.0162    2.63    0.189 
 8 FGF7            8        0.0118    1.97    0.189 
 9 CALM3           9        0.0116    1.72    0.189 
10 RELN           10        0.0110    1.94    0.189 
11 POSTN          11        0.0108    1.65    0.189 
12 TGFB1          12        0.0103    1.70    0.189 
13 IGF2           13        0.00989   1.43    0.189 
14 FGF2           14        0.00917   1.53    0.189 
15 COL10A1        15        0.00872   1.34    0.189 
