# TGF-β response score: specificity pre-specification

**Written 2026-09-02, before any of the models below were fitted.** The decision
rule in §7 is fixed here so that it cannot be chosen after seeing which arm wins.
Results land in `_m/stats_data/tgfb_specificity_*.tsv`; this file is not edited
once those exist, except to record the verdict in a clearly marked postscript.

---

## 1. Biological claim under test

TGF-β signalling in lung pericytes is associated with **lower basement-membrane
output**, and this is a statement about TGF-β/SMAD signalling rather than about a
generic activation or dissociation-stress program.

## 2. Primary hypothesis

The negative association between the TGF-β response score and the
basement-membrane score is carried by the **SMAD-proximal** genes of the panel,
and survives adjustment for the immediate-early/mechano-responsive genes.

Falsifiable, and directional: β_SMAD < 0 after adjustment for the IEG arm.

## 3. Why the current panel cannot answer this

`TGFB_RESPONSE` is 17 genes whose variance is carried almost entirely by six
well-detected genes that are **not** TGF-β-specific:

| Arm | Genes | Mean pericyte detection |
| --- | ----- | ----------------------: |
| IEG / mechano | JUNB, ID1, ID2, ID3, CCN1, CCN2, CDKN1A | **32.2 %** |
| SMAD-proximal | SMAD7, SKIL, SKI, PMEPA1, KLF10, BAMBI, TGIF1 | **9.2 %** |

`sc.tl.score_genes` weights genes by their observed variance, so a panel score is
dominated by its best-detected members. JUNB alone (56.5 %) is detected more
often than the entire SMAD arm combined. AP-1 (JUNB), the ID family (SMAD1/5/8,
i.e. **BMP**, not TGF-β/SMAD2/3), the CCN family (YAP/TAZ, mechanical) and
CDKN1A (p53/stress) all respond to serum, matrix stiffness and **warm tissue
dissociation** — the canonical single-cell artifact.

The existing `noECM` arm drops only TGFBI and does not touch this.

## 4. Arm definitions — assigned by pharmacology, before fitting

**SMAD-proximal (7).** Direct SMAD2/3 transcriptional targets forming the
canonical negative-feedback module: `SMAD7, SKIL, SKI, PMEPA1, KLF10, BAMBI,
TGIF1`. SMAD7/SKIL/SKI/TGIF1 are the inhibitory-feedback core; PMEPA1 and BAMBI
are direct SMAD2/3 targets; KLF10 (TIEG) was cloned as a TGF-β-inducible gene.

**IEG / mechano (7).** `JUNB, ID1, ID2, ID3, CCN1, CCN2, CDKN1A`. TGF-β-inducible
in the literature, but each has a dominant non-TGF-β route: AP-1 (JUNB), BMP/
SMAD1-5-8 (ID1/2/3), YAP/TAZ and strain (CCN1/CCN2), p53 and stress (CDKN1A).

**Unassigned (3), in neither arm.** `TGFBI` (an ECM gene; already the `noECM`
arm), `SERPINE2` and `SNAI1` (both broadly EMT/stress-responsive, and both under
2.2 % detected). Forcing them into either arm would blur the contrast the design
exists to draw.

⚠️ **The arms are confounded with detection by construction** — that is the
concern restated, not a flaw in the split. It is the reason §6 is mandatory
rather than optional.

## 5. Statistical unit and design

Unchanged from the existing TGF-β block in `04.bm_state_stats.R`, so results are
directly comparable to the reported β = −0.196:

- **Unit:** donor × Leiden-cluster pseudobulk, **214 units from 95 donors**.
- **Model:** `outcome_z ~ predictor_z + depth + (1 | study) + (1 | donor_id)`,
  `lmerTest`, Satterthwaite df.
- **Outcomes:** `basement_membrane_score_z` (primary) and
  `bm_minus_fibrillar_z` (secondary).
- **Grouping:** `pericyte_state` (Leiden, panel-independent), never
  `state_program` — the latter is derived from the BM score and would be circular.

## 6. Minimal primary analysis

**Test A — head-to-head, the parameter that answers the question.**

```
bm_score_z ~ tgfb_smad_score_z + tgfb_ieg_score_z + depth + (1|study) + (1|donor_id)
```

> **The parameter answering the scientific question is the coefficient on
> `tgfb_smad_score_z`, adjusted for `tgfb_ieg_score_z`.**

Each arm is also fitted alone, so attenuation between the marginal and adjusted
fits is readable.

**Test B — detection-matched empirical null. Mandatory.**

A null SMAD arm is uninterpretable on its own: the arm is 3.5× sparser, and we do
not otherwise know what effect size a 7-gene panel at ~9 % detection *can*
produce. So for each arm, draw **K = 1,000** random 7-gene panels matched
gene-by-gene on pericyte detection rate (nearest neighbour on `detect_frac`,
sampled without replacement), excluding every gene in any matrix, state, or
TGF-β panel. Score each with the **identical** `sc.tl.score_genes` call and fit
the identical model.

This yields, per arm:

- an **empirical two-sided p** for the observed β against panels of equal sparsity;
- a **power statement** — the fraction of matched-null panels reaching
  |β| ≥ 0.196, i.e. whether a panel this sparse could have detected the reported
  effect at all.

Test B is what converts "the SMAD arm is null" from uninterpretable into either
*null with adequate power* (which refutes the claim) or *underpowered* (which
leaves it unproven, not refuted). Without it the analysis cannot conclude either way.

**Test C — between-study variance per arm.** If the IEG arm is a dissociation
artifact it should partition far more variance to `study` (protocol) than the
SMAD arm does. Report the `study` variance component of each single-arm fit.
This is direct, cheap evidence on the technical hypothesis.

## 7. Decision rule — fixed in advance

| Outcome of Test A + B | Verdict | Consequence for the manuscript |
| --------------------- | ------- | ------------------------------ |
| β_SMAD < 0, outside its matched null, survives IEG adjustment | **TGF-β-specific** | Claim stands as written. |
| β_IEG < 0 survives, β_SMAD inside its null **and** null is well-powered | **Generic activation** | Retract the TGF-β framing; report as an activation/stress-vs-BM association. |
| β_SMAD inside its null, null **underpowered** | **Unproven** | Claim narrows to "TGF-β response *panel*", specificity explicitly undetermined; NicheNet carries the mechanism. |
| Both arms negative and outside their nulls | **Shared, not separable** | Report both; do not claim SMAD-specificity. |

## 8. Sensitivity analyses

1. **Within-donor cell-level arm** (per-donor Spearman + one-sample Wilcoxon),
   mirroring the existing structure. Matters because the marginal correlation
   (r = −0.457) already shrinks to −0.196 adjusted — more than half of the raw
   association is between-donor/between-study covariation.
2. **Leave-one-gene-out** on the full 17-gene panel (17 refits). Tests whether
   any single gene — JUNB above all — carries the reported association.
3. **`bm_minus_fibrillar_z`** as a parallel outcome throughout.
4. **Complexity floor**: refit on cells above the median `n_genes`, where sparse
   panel scores are least fragile.

## 9. Orthogonal validation — already in hand, and independent of this concern

The NicheNet arm (`nichenet_bm/`, `nichenet_fibrillar/`) uses **no response score
at all**: it ranks ligands by prior-network regulatory potential toward each
matrix target set. TGFB2 is perm *z* = 23.2 (BH 0.0064) toward fibrillar collagen
versus 2.94 (BH 0.235) toward BM.

**This matters for how much is at stake.** If the SMAD arm fails, the ligand-side
evidence still stands, and the surviving claim is "TGF-β ligand activity is
directed at fibrillar collagen rather than basement membrane" — a mechanism
statement that never depended on the response panel. Test B decides how much of
the *cell-state* half survives, not whether the whole result does.

## 10. Main figure-worthy result

A forest of β for {full panel, SMAD alone, IEG alone, SMAD | IEG, IEG | SMAD} ×
{BM, BM − fibrillar}, with each arm's **detection-matched null interval shaded
behind it**. The shaded null is the panel's argument: it shows the reader
directly that a sparse panel's null is wide, and whether the observed effect
clears it. Destination: a new panel E in `figureS_bm_associations`.

## 11. Reviewer objections

| Objection | Response |
| --------- | -------- |
| The arms were chosen to give the wanted answer | Assigned by pharmacology and written here before fitting; the file is committed ahead of the results. |
| The SMAD arm is null only because it is sparse | Exactly what Test B measures, and the pre-specified reason a null is not read as a refutation. |
| A generic well-detected panel would also correlate with BM | Test B's IEG-matched null answers this directly. |
| Cells are not independent | Donor × cluster pseudobulk, 214 units / 95 donors, plus a within-donor arm. |
| Multiple testing | 2 arms × 2 outcomes = 4 primary tests, BH within the `tgfb` family, matching the existing `bh_family` convention. |
| Dissociation stress explains everything | Test C: if so, the IEG arm should partition markedly more variance to `study`. |
| Circular gene selection | Neither arm overlaps any matrix or state panel; `_assert_tgfb_disjoint()` enforces this at import. |

## 12. What would fail the hypothesis

β_SMAD ≥ 0, or inside its detection-matched null while β_IEG is negative and
outside its own — with the null shown to be well-powered. That result would mean
the reported association is an activation/stress program, and the TGF-β framing
of the cell-state half of the claim comes out.
