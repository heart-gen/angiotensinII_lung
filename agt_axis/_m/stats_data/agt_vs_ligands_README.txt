AGT versus other pericyte ligands -- generated summary
AGT point-estimate rank in the frozen NicheNet run: 11
Rank bootstrap run: TRUE

NOTE: (A) rank and (C) target convergence inherit the NicheNet prior network
and are hypothesis-generating. (B) co-expression is measured in this dataset.

   ligand n_agt_targets n_other_targets n_shared   jaccard     hyper_p universe
   <char>         <int>           <int>    <int>     <num>       <num>    <int>
1:  TGFB1            12              13        6 0.3157895 0.793175394       24
2:  TGFB2            12              11        5 0.2777778 0.793175394       24
3:   CCN2            12              12       10 0.7142857 0.001664475       24
    hyper_p_BH
         <num>
1: 0.793175394
2: 0.793175394
3: 0.004993425

                   ccc_group partner partial_rho p_value n_donors  p_BH
                      <char>  <char>       <num>   <num>    <num> <num>
 1:   Subpleural fibroblasts   PDGFB          NA      NA       12    NA
 2:    EC aerocyte capillary   TGFB3  -0.6788461       0      137     0
 3:      Classical monocytes   TGFB3   0.7281710       0      339     0
 4:      Classical monocytes   PDGFB  -0.7663005       0      339     0
 5:                      DC2   TGFB3   0.7744921       0      353     0
 6:                      DC2    CCN2   0.4490036       0      353     0
 7: Interstitial macrophages   TGFB3  -0.5544078       0      344     0
 8: Interstitial macrophages    CCN2  -0.4364656       0      344     0
 9:               Mast cells   TGFB2  -0.7244475       0      238     0
10:               Mast cells   TGFB3  -0.6656859       0      238     0
