Snapshot taken 2026-09-01, immediately BEFORE re-running 00.state_discovery.py
with the expanded 20-gene basement-membrane panel (basement_membrane/_h/bm_panels.py).

These files are the 13-gene-BM-panel state model: every published BM number and
every state_program label in the manuscript as of this date was computed from
them. Kept so the invariant check after the re-run is possible:

  - pericyte_state (Leiden on X_pca_harmony) must be IDENTICAL for all 11,680
    cells -- clustering does not depend on the marker panels.
  - the five original program scores must be bit-identical (max |diff| = 0).
  - only basement_membrane_score may move.

The state gate (basement_membrane/_m/state_gate_summary.tsv, run 2026-09-01
before this snapshot) had already established that the 20-gene panel leaves the
cluster -> program map unchanged: metric A = 0.000% of cells, metric B = 0
clusters flipped.

See also _m_backup_5panel/ (the pre-2026-07-21 five-panel model).
