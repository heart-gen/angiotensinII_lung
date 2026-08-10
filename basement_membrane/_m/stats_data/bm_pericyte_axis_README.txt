BM pericyte axis -- generated summary
Donor x cluster units (>=5 cells): 214; donors: 95
Variance components: study SD 0.000, residual SD 0.907 -> ok

Grouping variable is pericyte_state (panel-independent Leiden clusters).
state_program is NOT used as an outcome grouping: after the gate escalated it
is derived from the BM score, so testing against it would be circular.

BM vs fibrillar (donor x cluster pseudobulk):
                                            comparison pearson_r pearson_ci_lo
                                                <char>     <num>         <num>
1:            basement_membrane vs fibrillar_ecm_score 0.2593153     0.1297099
2: basement_membrane vs fibroblast_like_noCOL4A1_score 0.2584707     0.1288197
3:          basement_membrane vs fibroblast_like_score 0.3903410     0.2703791
   pearson_ci_hi    pearson_p spearman_rho   spearman_p n_units
           <num>        <num>        <num>        <num>   <int>
1:     0.3802088 1.245046e-04    0.2386003 4.453588e-04     214
2:     0.3794341 1.312285e-04    0.2087917 2.178351e-03     214
3:     0.4983674 3.359348e-09    0.3301946 9.006730e-07     214

Direct AGTR1 vs BM (denoised lens is the readout; raw/detection are the
dropout-sensitive comparison, NOT the answer):
                     outcome         lens   estimate         SE       df
                      <char>       <char>      <num>      <num>    <num>
1: basement_membrane_score_z   AGTR1_scvi 0.02180869 0.06670510 190.0879
2: basement_membrane_score_z   AGTR1_expr 0.08881811 0.07089932 210.6897
3: basement_membrane_score_z AGTR1_detect 0.10939161 0.07600907 210.4516
4:                  bm_resid   AGTR1_scvi 0.03170310 0.05816112 212.0000
5:                  bm_resid   AGTR1_expr 0.05727642 0.05806877 212.0000
6:                  bm_resid AGTR1_detect 0.04589860 0.05811643 212.0000
     t_ratio   p_value depth_adjusted n_units n_donors      p_BH
       <num>     <num>         <lgcl>   <int>    <int>     <num>
1: 0.3269419 0.7440716           TRUE     214       95 0.7440716
2: 1.2527357 0.2116905           TRUE     214       95 0.3175358
3: 1.4391916 0.1515819           TRUE     214       95 0.3175358
4: 0.5450909 0.5862640          FALSE     214       95 0.5862640
5: 0.9863550 0.3250833          FALSE     214       95 0.5862640
6: 0.7897698 0.4305447          FALSE     214       95 0.5862640

Within-donor (cell-level Spearman per donor, >=20 cells, one-sample test):
           lens n_donors   median_rho         q25        q75  p_wilcox
         <char>    <int>        <num>       <num>      <num>     <num>
1:   AGTR1_scvi       69 -0.015972445 -0.08905322 0.07967829 0.6324255
2:   AGTR1_expr       68  0.005213482 -0.08801746 0.07508115 0.8092769
3: AGTR1_detect       68  0.034538796 -0.05559512 0.10347627 0.2621815
        p_BH
       <num>
1: 0.8092769
2: 0.8092769
3: 0.7865446
