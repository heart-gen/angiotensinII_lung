Fibrillar-collagen ambient controls -- generated summary
Cohort: Healthy donors only | units: 2329 | cell types: 22 | donors: 220
Ambient reference groups: Alveolar macrophages, EC aerocyte capillary, EC general capillary, EC venous pulmonary, EC venous systemic, Interstitial macrophages

NOTE: the soupX layer in pericyte_states.h5ad is bit-identical to counts
(verified 2026-08-10). It carries no ambient correction and is NOT used.

(1) Pericyte minus ambient-floor, log1p(CP10K), expr lens:
       gene           block estimate     p_BH
     <char>          <char>    <num>    <num>
 1:   SFTPB  ambient_tracer   -0.017  7.0e-01
 2: SCGB3A2  ambient_tracer   -0.021  6.6e-01
 3: SCGB1A1  ambient_tracer   -0.049  4.0e-01
 4:   SFTPC  ambient_tracer   -0.064  2.8e-01
 5:   PTPRC  ambient_tracer   -0.370  3.2e-05
 6:  COL1A2  fibrillar_core    1.729 8.4e-198
 7:  COL3A1  fibrillar_core    0.643  3.1e-34
 8:  COL1A1  fibrillar_core    0.245  2.7e-07
 9:  COL5A2 fibrillar_minor    1.192 7.2e-123
10:  COL5A3 fibrillar_minor    1.112 3.2e-184
11:  COL5A1 fibrillar_minor    0.455  2.0e-22
12: COL11A2 fibrillar_minor   -0.014  4.0e-01
13: COL11A1 fibrillar_minor   -0.016  3.8e-01

(2) Ambient-burden slopes (soup_burden = pericyte off-lineage tracer z;
    fibexpr_ = same gene in that donor's fibroblasts; fib_frac = fibroblast
    share of the niche). Positive and significant => consistent with soup.
       gene        term estimate    p_BH
     <char>      <char>    <num>   <num>
 1: COL11A1    fib_frac    0.064 7.7e-02
 2: COL11A1    fibexpr_    0.070 6.1e-03
 3: COL11A1 soup_burden   -0.005 4.5e-01
 4: COL11A2    fib_frac   -0.078 4.6e-01
 5: COL11A2    fibexpr_    0.268 9.3e-02
 6: COL11A2 soup_burden   -0.001 9.0e-01
 7:  COL1A1    fib_frac    2.395 7.7e-02
 8:  COL1A1    fibexpr_    0.395 7.7e-03
 9:  COL1A1 soup_burden    0.070 9.0e-01
10:  COL1A2    fib_frac    1.479 3.0e-01
11:  COL1A2    fibexpr_    1.231 1.9e-06
12:  COL1A2 soup_burden    0.211 4.5e-01
13:  COL3A1    fib_frac    1.051 4.6e-01
14:  COL3A1    fibexpr_    0.869 2.2e-06
15:  COL3A1 soup_burden    0.081 9.0e-01
16:  COL5A1    fib_frac    0.109 8.8e-01
17:  COL5A1    fibexpr_    0.602 2.6e-05
18:  COL5A1 soup_burden   -0.042 9.0e-01
19:  COL5A2    fib_frac   -0.198 8.8e-01
20:  COL5A2    fibexpr_    0.576 4.0e-04
21:  COL5A2 soup_burden   -0.020 9.0e-01
22:  COL5A3    fib_frac    1.526 3.0e-01
23:  COL5A3    fibexpr_    0.792 3.6e-02
24:  COL5A3 soup_burden   -0.048 9.0e-01
       gene        term estimate    p_BH

(3) COL1A1 - COL1A2 within unit, by cell type:
                    ccc_group              role emmean n_units
                       <fctr>            <char>  <num>   <int>
 1:                 Pericytes         Pericytes -1.477      49
 2:      Alveolar fibroblasts       Fibroblasts -0.681      98
 3:   Adventitial fibroblasts       Fibroblasts -0.512      94
 4:    Vascular smooth muscle             Other -0.327      81
 5: Peribronchial fibroblasts       Fibroblasts -0.299      57
 6:            Myofibroblasts       Fibroblasts -0.207      33
 7:       EC venous pulmonary Ambient reference -0.107      26
 8:        EC venous systemic Ambient reference -0.091      36
 9:      Alveolar macrophages Ambient reference -0.076     103
10:  Interstitial macrophages Ambient reference -0.048      73
11:     EC aerocyte capillary Ambient reference -0.037      28
12:      EC general capillary Ambient reference -0.032      45
13:   Non-classical monocytes             Other -0.032      45
14:       Classical monocytes             Other -0.027      63
15:                       DC2             Other -0.015      47
16:              Lymphatic EC             Other -0.003      42
17:                Mast cells             Other  0.024      36
18:                       AT1             Other  0.048      73
19:            AT2_AGTR2undet             Other  0.085     118
20:              AT2_AGTR2det             Other  0.091      45
21:     Transitional Club-AT2             Other  0.104      70
22:    Subpleural fibroblasts       Fibroblasts  0.134       8
                    ccc_group              role emmean n_units
