AGT versus other pericyte ligands -- generated summary
AGT point-estimate rank in the frozen NicheNet run: 11
Rank bootstrap run: TRUE

NOTE: (A) rank and (C) target overlap inherit the NicheNet prior network and
are hypothesis-generating. (B) co-expression is the only block measured in
this dataset, so (B) and (C) are NOT two independent lines of evidence.
In (C), read pair_pctile and degree_p; hyper_p_MISCALIBRATED assumes uniform
target draws, which this 24-gene shortlist badly violates.

   ligand n_agt_targets n_other_targets n_shared   jaccard universe
   <char>         <int>           <int>    <int>     <num>    <int>
1:  TGFB1            12              13        6 0.3157895       24
2:  TGFB2            12              11        5 0.2777778       24
3:   CCN2            12              12       10 0.7142857       24
   n_ligand_pairs jaccard_median_all_pairs pair_pctile  degree_p
            <int>                    <num>       <num>     <num>
1:            435                0.5833333    24.13793 0.9993000
2:            435                0.5833333    18.62069 0.9989501
3:            435                0.5833333    75.17241 0.2041398
   hyper_p_MISCALIBRATED degree_p_BH hyper_p_MISCALIBRATED_BH
                   <num>       <num>                    <num>
1:           0.793175394   0.9993000              0.793175394
2:           0.793175394   0.9993000              0.793175394
3:           0.001664475   0.6124194              0.004993425

                    ccc_group partner partial_rho      p_value n_donors
                       <char>  <char>       <num>        <num>    <num>
 1:      Alveolar fibroblasts   PDGFB   0.4498635 5.589890e-10      176
 2:              AT2_AGTR2det   PDGFB  -0.5212166 2.273468e-09      118
 3: Peribronchial fibroblasts   PDGFB   0.4159109 8.175108e-07      133
 4:              AT2_AGTR2det   TGFB3  -0.3683273 4.601250e-05      118
 5:            Myofibroblasts   PDGFB  -0.4462189 3.591181e-04       61
 6:   Adventitial fibroblasts   TGFB3   0.2566701 4.545572e-04      184
 7:              AT2_AGTR2det   TGFB2  -0.2332599 1.116944e-02      118
 8:   Adventitial fibroblasts   PDGFB   0.1836778 1.266302e-02      184
 9:      Alveolar fibroblasts    CCN2   0.1879045 1.261309e-02      176
10:                 Pericytes    CCN2   0.2490151 1.058561e-02      105
    agt_detect_group n_cells_group p_underflow tested excluded_reason
               <num>         <int>      <lgcl> <lgcl>          <char>
 1:       0.01711990         16180       FALSE   TRUE                
 2:       0.01119751          3215       FALSE   TRUE                
 3:       0.03343166          6102       FALSE   TRUE                
 4:       0.01119751          3215       FALSE   TRUE                
 5:       0.02192982          1824       FALSE   TRUE                
 6:       0.04461756         12708       FALSE   TRUE                
 7:       0.01119751          3215       FALSE   TRUE                
 8:       0.04461756         12708       FALSE   TRUE                
 9:       0.01711990         16180       FALSE   TRUE                
10:       0.01846417          8178       FALSE   TRUE                
            p_BH
           <num>
 1: 1.956462e-08
 2: 3.978570e-08
 3: 9.537626e-06
 4: 4.026094e-04
 5: 2.513827e-03
 6: 2.651583e-03
 7: 4.432058e-02
 8: 4.432058e-02
 9: 4.432058e-02
10: 4.432058e-02
