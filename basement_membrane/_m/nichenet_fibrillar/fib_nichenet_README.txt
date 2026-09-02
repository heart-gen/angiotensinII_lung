FIB-target ligand ranking vs the frozen ECM-program ranking
Spearman rank correlation: 0.550
Top-20 overlap: 10/20
TGF-beta ligands in the FIB ranking:
# A tibble: 3 x 5
  test_ligand  rank aupr_corrected perm_z perm_p_BH
  <chr>       <int>          <dbl>  <dbl>     <dbl>
1 TGFB2           3        0.176   23.2     0.00642
2 TGFB1          16        0.0226   2.89    0.110  
3 TGFB3          84        0.00421  0.362   0.128  

Top 15 ligands toward fibrillar collagen targets:
# A tibble: 15 x 5
   test_ligand  rank aupr_corrected perm_z perm_p_BH
   <chr>       <int>          <dbl>  <dbl>     <dbl>
 1 COL4A1          1         0.337   40.1    0.00642
 2 COPA            2         0.267   32.1    0.00642
 3 TGFB2           3         0.176   23.2    0.00642
 4 COL1A1          4         0.174   23.6    0.00917
 5 COMP            5         0.172   22.5    0.00642
 6 TIMP2           6         0.172   20.4    0.00642
 7 COL5A1          7         0.170   20.2    0.00917
 8 RELN            8         0.116   17.0    0.0481 
 9 COL10A1         9         0.0883  11.4    0.0525 
10 BMP4           10         0.0466   7.13   0.0481 
11 INHBA          11         0.0458   6.77   0.0481 
12 MMP14          12         0.0453   5.55   0.0588 
13 CCN2           13         0.0315   4.06   0.0940 
14 CALM3          14         0.0306   3.53   0.103  
15 PDGFA          15         0.0303   4.31   0.0790 
