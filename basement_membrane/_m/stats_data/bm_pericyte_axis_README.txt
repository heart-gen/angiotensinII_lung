BM pericyte axis -- generated summary
Donor x cluster units (>=5 cells): 214; donors: 95
Variance components: study SD 0.000, residual SD 0.912 -> ok

Grouping variable is pericyte_state (panel-independent Leiden clusters).
state_program is NOT used as an outcome grouping: after the gate escalated it
is derived from the BM score, so testing against it would be circular.

BM vs fibrillar (donor x cluster pseudobulk):
                                            comparison pearson_r pearson_ci_lo
                                                <char>     <num>         <num>
1:            basement_membrane vs fibrillar_ecm_score 0.1687275    0.03541232
2: basement_membrane vs fibroblast_like_noCOL4A1_score 0.1500919    0.01630353
3:          basement_membrane vs fibroblast_like_score 0.2768100    0.14819519
   pearson_ci_hi    pearson_p spearman_rho  spearman_p n_units
           <num>        <num>        <num>       <num>   <int>
1:     0.2961426 1.345247e-02    0.1483040 0.030180023     214
2:     0.2786003 2.814595e-02    0.1052784 0.124636856     214
3:     0.3962171 4.022414e-05    0.2210316 0.001162914     214

Matrix-vs-predictor associations. Outcomes are BOTH matrix categories, their
difference (taken on the raw scores, then standardized) and two
orthogonalizations. Depth is a covariate in every pseudobulk model and is
partialled out of every within-donor rho: all score_genes scores rise with
capture depth, so an unadjusted correlation between two of them is the
expected result under pure technical variation.

WHICH LENS ARBITRATES (changed 2026-09-02). The denoised lens was originally
designated the AGTR1 readout. It is NOT the arbiter any more: it is a model
output that can in principle manufacture or erase an association by borrowing
across genes, and the first version of it was trained on the soupX float
matrix and failed its validity gate (donor rho 0.014, p 0.91). It is now
reported as a CONCORDANT SENSITIVITY. Between-group contrasts are arbitrated
by 10.agtr1_count_models.R -- AGTR1 integer counts from raw/X as the response
of an NB GLMM with offset(log(library size)), so dropout sits inside the
likelihood and nothing is imputed. raw/detection remain the dropout-sensitive
comparison. Note the matrix scores themselves have no denoised counterpart, so
EVERY row here still has one depth-affected side, which is the other reason
the offset-bearing count model carries the claim.
                       outcome                 predictor family    estimate
                        <char>                    <char> <char>       <num>
 1:  basement_membrane_score_z                AGTR1_scvi  agtr1  0.12586104
 2:  basement_membrane_score_z                AGTR1_expr  agtr1  0.16953707
 3:  basement_membrane_score_z              AGTR1_detect  agtr1  0.20942094
 4: fibrillar_collagen_score_z                AGTR1_scvi  agtr1 -0.13445426
 5: fibrillar_collagen_score_z                AGTR1_expr  agtr1 -0.05428682
 6: fibrillar_collagen_score_z              AGTR1_detect  agtr1 -0.02962787
 7:       bm_minus_fibrillar_z                AGTR1_scvi  agtr1  0.31299294
 8:       bm_minus_fibrillar_z                AGTR1_expr  agtr1  0.19656985
 9:       bm_minus_fibrillar_z              AGTR1_detect  agtr1  0.16482494
10:                   bm_resid                AGTR1_scvi  agtr1  0.09869076
11:                   bm_resid                AGTR1_expr  agtr1  0.12474511
12:                   bm_resid              AGTR1_detect  agtr1  0.11806967
13:          bm_resid_collagen                AGTR1_scvi  agtr1  0.17836082
14:          bm_resid_collagen                AGTR1_expr  agtr1  0.16053390
15:          bm_resid_collagen              AGTR1_detect  agtr1  0.15863770
16:  basement_membrane_score_z       tgfb_response_score   tgfb -0.19640313
17:  basement_membrane_score_z tgfb_response_noECM_score   tgfb -0.19399552
18:  basement_membrane_score_z       tgfb_receptor_score   tgfb -0.01481949
19: fibrillar_collagen_score_z       tgfb_response_score   tgfb -0.07285565
20: fibrillar_collagen_score_z tgfb_response_noECM_score   tgfb -0.08991185
21: fibrillar_collagen_score_z       tgfb_receptor_score   tgfb  0.05669635
22:       bm_minus_fibrillar_z       tgfb_response_score   tgfb -0.08402273
23:       bm_minus_fibrillar_z tgfb_response_noECM_score   tgfb -0.06091967
24:       bm_minus_fibrillar_z       tgfb_receptor_score   tgfb -0.07206932
25:                   bm_resid       tgfb_response_score   tgfb -0.17420163
26:                   bm_resid tgfb_response_noECM_score   tgfb -0.16889543
27:                   bm_resid       tgfb_receptor_score   tgfb -0.02187767
28:          bm_resid_collagen       tgfb_response_score   tgfb -0.15411975
29:          bm_resid_collagen tgfb_response_noECM_score   tgfb -0.14557040
30:          bm_resid_collagen       tgfb_receptor_score   tgfb -0.03347436
                       outcome                 predictor family    estimate
            SE       df    t_ratio      p_value depth_adjusted n_units n_donors
         <num>    <num>      <num>        <num>         <lgcl>   <int>    <int>
 1: 0.06766576 194.5246  1.8600402 6.438947e-02           TRUE     214       95
 2: 0.07001725 210.8952  2.4213614 1.630843e-02           TRUE     214       95
 3: 0.07478746 208.4804  2.8002148 5.586839e-03           TRUE     214       95
 4: 0.06690631 210.0517 -2.0095902 4.575405e-02           TRUE     214       95
 5: 0.06872898 205.6773 -0.7898680 4.305145e-01           TRUE     214       95
 6: 0.07444179 210.1228 -0.3980006 6.910338e-01           TRUE     214       95
 7: 0.06528566 210.9628  4.7942063 3.081321e-06           TRUE     214       95
 8: 0.06753085 203.7377  2.9108158 4.005926e-03           TRUE     214       95
 9: 0.07427925 209.8778  2.2189905 2.755831e-02           TRUE     214       95
10: 0.05785751 212.0000  1.7057552 8.951825e-02          FALSE     214       95
11: 0.05761972 212.0000  2.1649725 3.150665e-02          FALSE     214       95
12: 0.05768603 212.0000  2.0467636 4.191505e-02          FALSE     214       95
13: 0.05596710 187.8239  3.1868871 1.684534e-03          FALSE     214       95
14: 0.05552796 211.6891  2.8910460 4.239946e-03          FALSE     214       95
15: 0.05554629 211.9983  2.8559551 4.717458e-03          FALSE     214       95
16: 0.06830555 210.8685 -2.8753615 4.449556e-03           TRUE     214       95
17: 0.06822635 210.9935 -2.8434105 4.902097e-03           TRUE     214       95
18: 0.07030150 210.4357 -0.2107990 8.332481e-01           TRUE     214       95
19: 0.06695423 201.8757 -1.0881412 2.778302e-01           TRUE     214       95
20: 0.06679380 202.1379 -1.3461108 1.797745e-01           TRUE     214       95
21: 0.06753625 200.0888  0.8394951 4.021934e-01           TRUE     214       95
22: 0.06752045 200.3239 -1.2444043 2.148048e-01           TRUE     214       95
23: 0.06753420 200.2648 -0.9020567 3.681104e-01           TRUE     214       95
24: 0.06725174 195.4630 -1.0716350 2.852054e-01           TRUE     214       95
25: 0.05701133 212.0000 -3.0555613 2.534854e-03          FALSE     214       95
26: 0.05708660 212.0000 -2.9585828 3.441593e-03          FALSE     214       95
27: 0.05823381 212.0000 -0.3756867 7.075255e-01          FALSE     214       95
28: 0.05561279 209.5853 -2.7713005 6.085473e-03          FALSE     214       95
29: 0.05572565 208.9176 -2.6122692 9.647434e-03          FALSE     214       95
30: 0.05661169 210.4815 -0.5912976 5.549555e-01          FALSE     214       95
            SE       df    t_ratio      p_value depth_adjusted n_units n_donors
            p_BH
           <num>
 1: 6.438947e-02
 2: 2.446265e-02
 3: 1.676052e-02
 4: 1.372622e-01
 5: 6.457718e-01
 6: 6.910338e-01
 7: 9.243962e-06
 8: 6.008889e-03
 9: 2.755831e-02
10: 8.951825e-02
11: 6.287257e-02
12: 6.287257e-02
13: 4.717458e-03
14: 4.717458e-03
15: 4.717458e-03
16: 7.353145e-03
17: 7.353145e-03
18: 8.332481e-01
19: 4.021934e-01
20: 4.021934e-01
21: 4.021934e-01
22: 3.681104e-01
23: 3.681104e-01
24: 3.681104e-01
25: 5.162389e-03
26: 5.162389e-03
27: 7.075255e-01
28: 1.447115e-02
29: 1.447115e-02
30: 5.549555e-01
            p_BH

Within-donor (cell-level Spearman per donor, >=20 cells, one-sample test):
                     outcome                 predictor family n_donors
                      <char>                    <char> <char>    <int>
 1:  basement_membrane_score                AGTR1_scvi  agtr1       69
 2:  basement_membrane_score                AGTR1_expr  agtr1       68
 3:  basement_membrane_score              AGTR1_detect  agtr1       68
 4: fibrillar_collagen_score                AGTR1_scvi  agtr1       69
 5: fibrillar_collagen_score                AGTR1_expr  agtr1       68
 6: fibrillar_collagen_score              AGTR1_detect  agtr1       68
 7:       bm_minus_fibrillar                AGTR1_scvi  agtr1       69
 8:       bm_minus_fibrillar                AGTR1_expr  agtr1       68
 9:       bm_minus_fibrillar              AGTR1_detect  agtr1       68
10:  basement_membrane_score       tgfb_response_score   tgfb       69
11:  basement_membrane_score tgfb_response_noECM_score   tgfb       69
12:  basement_membrane_score       tgfb_receptor_score   tgfb       69
13: fibrillar_collagen_score       tgfb_response_score   tgfb       69
14: fibrillar_collagen_score tgfb_response_noECM_score   tgfb       69
15: fibrillar_collagen_score       tgfb_receptor_score   tgfb       69
16:       bm_minus_fibrillar       tgfb_response_score   tgfb       69
17:       bm_minus_fibrillar tgfb_response_noECM_score   tgfb       69
18:       bm_minus_fibrillar       tgfb_receptor_score   tgfb       69
      median_rho          q25        q75     p_wilcox n_donors_partial
           <num>        <num>      <num>        <num>            <int>
 1:  0.084976796 -0.027976555 0.17682117 4.321191e-05               69
 2:  0.039638258 -0.035930433 0.13533215 1.789605e-02               68
 3:  0.067008892 -0.020513300 0.14065585 4.891189e-03               68
 4: -0.046589551 -0.145094972 0.07042091 1.130980e-01               69
 5:  0.005040658 -0.067621645 0.07727464 5.432006e-01               68
 6:  0.019995792 -0.075524117 0.10135007 2.298736e-01               68
 7:  0.088140716  0.002651306 0.18785681 2.367066e-05               69
 8:  0.041346391 -0.052837663 0.10121413 1.467217e-01               68
 9:  0.038447788 -0.074749591 0.12291877 2.595902e-01               68
10: -0.036470588 -0.100170940 0.02142857 8.224926e-03               69
11: -0.034569123 -0.110263322 0.03113553 1.733047e-02               69
12: -0.014338575 -0.097755102 0.05071532 2.765208e-01               69
13:  0.004202959 -0.080039526 0.11284932 9.095549e-01               69
14: -0.007508266 -0.088266044 0.09343968 8.109832e-01               69
15:  0.029169526 -0.043478261 0.10684935 3.329524e-02               69
16: -0.023137481 -0.097774448 0.08454151 2.203176e-01               69
17: -0.021691352 -0.093369788 0.07746904 3.540633e-01               69
18: -0.035612245 -0.095183260 0.06684880 1.130980e-01               69
    median_partial_rho p_wilcox_partial      depth_control         p_BH
                 <num>            <num>             <char>        <num>
 1:       0.1105007644     8.112945e-07 log10_total_counts 1.296357e-04
 2:       0.0333808908     1.114415e-01 log10_total_counts 1.789605e-02
 3:       0.0354257196     1.114415e-01 log10_total_counts 7.336784e-03
 4:      -0.0206935519     2.007230e-01 log10_total_counts 3.392940e-01
 5:       0.0008151748     9.196927e-01 log10_total_counts 5.432006e-01
 6:      -0.0012597016     8.235138e-01 log10_total_counts 3.448104e-01
 7:       0.0950449384     1.205511e-05 log10_total_counts 7.101199e-05
 8:       0.0459549183     6.543543e-02 log10_total_counts 2.200825e-01
 9:       0.0431531020     8.653715e-02 log10_total_counts 2.595902e-01
10:      -0.0468964875     1.499292e-03 log10_total_counts 2.467478e-02
11:      -0.0504222067     2.687139e-03 log10_total_counts 2.599570e-02
12:       0.0002973228     8.953507e-01 log10_total_counts 2.765208e-01
13:       0.0021866435     7.832920e-01 log10_total_counts 9.095549e-01
14:      -0.0037888437     8.202683e-01 log10_total_counts 9.095549e-01
15:       0.0430345277     3.638209e-02 log10_total_counts 9.988573e-02
16:      -0.0381524697     9.768777e-02 log10_total_counts 3.304764e-01
17:      -0.0430546586     1.904025e-01 log10_total_counts 3.540633e-01
18:      -0.0147867734     2.158489e-01 log10_total_counts 3.304764e-01
    p_BH_partial
           <num>
 1: 2.433883e-06
 2: 1.114415e-01
 3: 1.114415e-01
 4: 6.021689e-01
 5: 9.196927e-01
 6: 9.196927e-01
 7: 3.616534e-05
 8: 8.653715e-02
 9: 8.653715e-02
10: 4.030708e-03
11: 4.030708e-03
12: 8.953507e-01
13: 8.202683e-01
14: 8.202683e-01
15: 1.091463e-01
16: 2.158489e-01
17: 2.158489e-01
18: 2.158489e-01
