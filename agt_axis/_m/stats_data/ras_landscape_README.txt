Local RAS landscape -- generated summary
Units (>=5 cells): 4376; cell types: 22; donors: 417
Detection threshold for 'step present': 0.05

Circular units (stratum defined by the gene being scored; excluded before
within-dataset standardization, retained with emmean = NA and the flag
circular_by_construction in ras_celltype_profile.tsv):
  AGTR2: AT2_AGTR2det, AT2_AGTR2undet

Cell types with an autonomous AGT->AngII->AT1R circuit: 0
Maximum REN (renin) detection across all cell types: 0.0229

Top 3 cell types per gene:
      gene                 ccc_group      emmean       detect rank_in_gene
    <char>                    <char>       <num>        <num>        <int>
 1:    ACE     EC aerocyte capillary  2.95862579 0.3513218905            1
 2:    ACE      EC general capillary  1.36186589 0.1441277232            2
 3:    ACE      Alveolar macrophages  0.72703169 0.1962545959            3
 4:   ACE2                 Pericytes  1.61245161 0.0301495177            1
 5:   ACE2              AT2_AGTR2det  0.77569886 0.0219937545            2
 6:   ACE2            AT2_AGTR2undet  0.58797150 0.0105319370            3
 7:    AGT    Vascular smooth muscle  2.62641483 0.0857823871            1
 8:    AGT   Adventitial fibroblasts  0.71162049 0.0358629683            2
 9:    AGT      Alveolar fibroblasts  0.36573471 0.0166316298            3
10:  AGTR1                 Pericytes  3.81033573 0.3417617893            1
11:  AGTR1      Alveolar fibroblasts  1.33328858 0.1007838401            2
12:  AGTR1   Adventitial fibroblasts  1.03325484 0.1068614963            3
13:  AGTR2 Peribronchial fibroblasts  0.58860920 0.0034689508            1
14:  AGTR2   Adventitial fibroblasts  0.14072422 0.0015104041            2
15:  AGTR2      EC general capillary  0.12824954 0.0007253343            3
16:   CCN2            Myofibroblasts  1.78225018 0.6383848922            1
17:   CCN2      Alveolar fibroblasts  1.73643180 0.5830725836            2
18:   CCN2 Peribronchial fibroblasts  1.28382246 0.4803280636            3
19:   CMA1                Mast cells  1.76897775 0.0250321021            1
20:   CMA1            AT2_AGTR2undet -0.05750213 0.0003600370            2
21:   CMA1    Vascular smooth muscle -0.06697232 0.0001860119            3
22:   CTSD      Alveolar macrophages  0.99791812 0.6057259661            1
23:   CTSD  Interstitial macrophages  0.91964736 0.6041899030            2
24:   CTSD                Mast cells  0.62450242 0.4758130500            3
25:   CTSG                Mast cells  2.49497661 0.1086125990            1
26:   CTSG      Alveolar macrophages -0.06257426 0.0086254041            2
27:   CTSG       Classical monocytes -0.07717060 0.0024572255            3
28:  ENPEP      Alveolar fibroblasts  2.41402833 0.1648242110            1
29:  ENPEP                 Pericytes  1.49910352 0.0940872105            2
30:  ENPEP    Vascular smooth muscle  1.25768160 0.0750275171            3
31:   LRP2              AT2_AGTR2det  2.66591867 0.3225349827            1
32:   LRP2            AT2_AGTR2undet  2.40254207 0.1938776190            2
33:   LRP2     Transitional Club-AT2  0.06247345 0.0311507345            3
34:   MAS1              AT2_AGTR2det  0.18101724 0.0025843932            1
35:   MAS1            AT2_AGTR2undet  0.14339423 0.0017090188            2
36:   MAS1  Interstitial macrophages  0.10049395 0.0026780139            3
37:    MME      Alveolar macrophages  2.14801037 0.4215318953            1
38:    MME      Alveolar fibroblasts  1.41137573 0.1643649794            2
39:    MME     Transitional Club-AT2  0.59756310 0.1251539200            3
40:  PDGFB     EC aerocyte capillary  2.48267983 0.2172802656            1
41:  PDGFB      EC general capillary  2.01153471 0.1589581563            2
42:  PDGFB        EC venous systemic  0.77233617 0.1204479730            3
43:    REN              AT2_AGTR2det  1.21904424 0.0229086378            1
44:    REN            AT2_AGTR2undet  1.04349507 0.0149797768            2
45:    REN      Alveolar fibroblasts  0.09598377 0.0058137364            3
46:  TGFB1       Classical monocytes  0.89785323 0.2987569945            1
47:  TGFB1                       DC2  0.71683703 0.4297263797            2
48:  TGFB1  Interstitial macrophages  0.71613617 0.4492535917            3
49:  TGFB2                       AT1  1.41317003 0.1697745601            1
50:  TGFB2            Myofibroblasts  1.20892182 0.1418318244            2
51:  TGFB2                 Pericytes  1.13457227 0.1213685610            3
52:  TGFB3 Peribronchial fibroblasts  2.07758604 0.2060123126            1
53:  TGFB3   Adventitial fibroblasts  1.27995848 0.1328849157            2
54:  TGFB3    Vascular smooth muscle  1.12331064 0.0846411316            3
      gene                 ccc_group      emmean       detect rank_in_gene

Circuit completeness by cell type:
                    ccc_group     ace_step chymase_step receptor_AT1
                       <char>        <num>        <num>        <num>
 1:   Adventitial fibroblasts 0.0450035126 1.456381e-03 1.068615e-01
 2:      Alveolar fibroblasts 0.0031325750 1.224915e-03 1.007838e-01
 3:      Alveolar macrophages 0.1962545959 8.625404e-03 6.935295e-04
 4:     EC aerocyte capillary 0.3513218905 2.320095e-03 3.133926e-04
 5:      EC general capillary 0.1441277232 1.397083e-03 5.649170e-04
 6:        EC venous systemic 0.1073745749 2.073171e-03 1.640481e-03
 7:  Interstitial macrophages 0.1446788199 4.867400e-03 3.622099e-04
 8:                Mast cells 0.0109545523 1.086126e-01 3.368736e-03
 9:                 Pericytes 0.0074979233 1.574045e-03 3.417618e-01
10:    Vascular smooth muscle 0.0096818447 1.447704e-03 4.734600e-02
11:                       AT1 0.0052467330 5.768836e-04 5.505183e-04
12:              AT2_AGTR2det 0.0112382633 2.115412e-03 3.110708e-04
13:            AT2_AGTR2undet 0.0065971290 1.183740e-03 1.053728e-03
14:       Classical monocytes 0.0040984704 2.457226e-03 4.228347e-04
15:                       DC2 0.0178615863 1.982305e-03 3.652174e-04
16:       EC venous pulmonary 0.0268010107 5.338844e-05 9.450082e-05
17:              Lymphatic EC 0.0080148891 9.861547e-04 3.695791e-04
18:            Myofibroblasts 0.0007355158 4.430660e-04 1.567441e-02
19:   Non-classical monocytes 0.0347488173 2.896873e-04 3.771801e-04
20: Peribronchial fibroblasts 0.0121024063 1.973879e-03 4.603838e-02
21:    Subpleural fibroblasts 0.0261806425 0.000000e+00 3.862286e-02
22:     Transitional Club-AT2 0.0229535972 1.257282e-03 7.207752e-04
                    ccc_group     ace_step chymase_step receptor_AT1
      renin_step    substrate n_steps_present has_substrate has_protease
           <num>        <num>           <num>        <lgcl>       <lgcl>
 1: 0.0052503965 0.0358629683               1         FALSE        FALSE
 2: 0.0058137364 0.0166316298               1         FALSE        FALSE
 3: 0.0022366866 0.0025025678               1         FALSE         TRUE
 4: 0.0015735302 0.0003137914               1         FALSE         TRUE
 5: 0.0013021699 0.0001113371               1         FALSE         TRUE
 6: 0.0035643426 0.0002105263               1         FALSE         TRUE
 7: 0.0016458751 0.0007768516               1         FALSE         TRUE
 8: 0.0004815629 0.0001915926               1         FALSE         TRUE
 9: 0.0030646123 0.0154677196               1         FALSE        FALSE
10: 0.0013379008 0.0857823871               1          TRUE        FALSE
11: 0.0028051823 0.0004296456               0         FALSE        FALSE
12: 0.0229086378 0.0077309274               0         FALSE        FALSE
13: 0.0149797768 0.0055193034               0         FALSE        FALSE
14: 0.0005432959 0.0005058451               0         FALSE        FALSE
15: 0.0007747377 0.0003379640               0         FALSE        FALSE
16: 0.0042815429 0.0001002110               0         FALSE        FALSE
17: 0.0024094150 0.0011898089               0         FALSE        FALSE
18: 0.0018271329 0.0104411672               0         FALSE        FALSE
19: 0.0011380641 0.0001374314               0         FALSE        FALSE
20: 0.0019135491 0.0263886846               0         FALSE        FALSE
21: 0.0000000000 0.0156461186               0         FALSE        FALSE
22: 0.0047275648 0.0047157771               0         FALSE        FALSE
      renin_step    substrate n_steps_present has_substrate has_protease
    has_receptor autonomous_circuit
          <lgcl>             <lgcl>
 1:         TRUE              FALSE
 2:         TRUE              FALSE
 3:        FALSE              FALSE
 4:        FALSE              FALSE
 5:        FALSE              FALSE
 6:        FALSE              FALSE
 7:        FALSE              FALSE
 8:        FALSE              FALSE
 9:         TRUE              FALSE
10:        FALSE              FALSE
11:        FALSE              FALSE
12:        FALSE              FALSE
13:        FALSE              FALSE
14:        FALSE              FALSE
15:        FALSE              FALSE
16:        FALSE              FALSE
17:        FALSE              FALSE
18:        FALSE              FALSE
19:        FALSE              FALSE
20:        FALSE              FALSE
21:        FALSE              FALSE
22:        FALSE              FALSE
    has_receptor autonomous_circuit
