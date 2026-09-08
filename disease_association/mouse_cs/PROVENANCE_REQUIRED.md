# `mouse_cs/` — the only wet-lab data in this repository, and it is unusable as it stands

**Status:** orphaned and **not citable**. Written 2026-09-08 for defect **P2-31**
(`writings/TODO.md`). This is a **bench task, not a code task** — nothing in this
file can be fixed by editing the analysis.

---

## Why this matters enough to have its own file

Every other result in this project comes from the same two atlases. This is the
only orthogonal, non-atlas evidence anywhere in the work, and it points the way
the current thesis predicts:

| Receptor | Room Air | Cigarette smoke | Fold | Welch *t* | df | *P* |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| AT1 | 1.012 (*n* = 5) | 0.797 (*n* = 4) | 0.79× | −1.52 | 5.0 | 0.19 |
| **AT2** | 1.006 (*n* = 5) | **1.653** (*n* = 3) | **1.64×** | +8.74 | 6.0 | **1.2 × 10⁻⁴** |

Recomputed 2026-09-08 directly from `_h/Mouse-receptor-data-with-smoke-exposure.csv`.
That is an **AT1R/AT2R balance shift under an injury exposure**.

It is also, right now, uncitable — and provenance recovery gets harder every
month, which is why this is logged as a live task rather than a footnote.

## Confirmed orphaned

`grep -rn` for `mouse_cs`, `Mouse-receptor`, `smoke_exposure` and `qPCR` across
every `.R`, `.py`, `.sh` and `.md` returns only this module's own files plus
ledger/briefing entries *about* it. No figure script, no table script, no
`figures/_m/figure_panel_manifest.tsv` entry, no `writings/` summary depends on
it. (`tables/_h/01.cohort_mouse.R` is the CELLxGENE mouse scRNA-seq from
`cross_species/` — a different dataset answering a different question.)

**Nothing breaks if this is never recovered.** Nothing improves either.

## The entire record

`_h/Mouse-receptor-data-with-smoke-exposure.csv` is three lines. In full:

```
,Rm Air,Rm Air,Rm Air,Rm Air,Rm Air,CS,CS,CS,CS,CS
AT1,1,1.03,0.8,0.98,1.25,0.46,0.81,1.03,0.89,
AT2,1.13,1.1,0.82,0.91,1.07,1.73,1.58,1.65,,
```

## What must be recovered before any use

1. **Housekeeping gene(s) and ΔΔCt values.** The table holds fold-changes only.
   Room-air values centre on 1.00 **because they are the normaliser**, so that
   group's spread (AT1 sd 0.161, AT2 sd 0.134) is partly a normalisation
   artifact and the Welch test above treats it as biological variance.
2. **Which *Agtr1* paralog was amplified.** Mouse has *Agtr1a* and *Agtr1b*; the
   label is just `AT1`. This is not a detail — the cross-species argument in
   `cross_species/` is specifically about ***Agtr1a***.
3. **Primer sequences**, strain, sex, age, exposure duration and dose.
4. **Biological vs technical replicates.** With *n* = 3–5, whether these are
   animals or wells decides whether there is an experiment here at all.
5. **What the blanks mean.** CS has 5 slots; AT1 fills 4 and AT2 fills 3.
   `01.plotting_cs.R` calls `drop_na()`, so the missing values disappear
   silently and the figure shows no *n*. Were those animals lost, excluded, or
   never run?

## Two limits that recovery cannot remove

- **Whole-lung qPCR cannot attribute anything to pericytes**, which is what this
  project claims. At best this supports a compartment-level AT1R/AT2R balance
  statement in mouse lung under smoke.
- ***n* = 3–5 with one exposure arm** supports a direction, not an effect size.

## Decision required

**Decide now whether this can be documented**, because the answer changes the
figure plan:

- **If recoverable** — it becomes the project's only orthogonal validation and
  is worth a panel. Fix `01.plotting_cs.R` to print *n* and to record rather
  than silently drop the blanks.
- **If not recoverable** — say so explicitly in the ledger and in
  `writings/pi_briefings/disease_association_legacy_cohorts/`, and ensure no
  figure plan assumes it. As of 2026-09-08 none does; that must stay true.

Do not cite the *P* = 1.2 × 10⁻⁴ above in any document until item 1 is resolved.
