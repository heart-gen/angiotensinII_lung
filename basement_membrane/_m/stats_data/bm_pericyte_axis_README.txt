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
