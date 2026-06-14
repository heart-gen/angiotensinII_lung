# ----------------------------
# Core review function
# ----------------------------
def metadata_review(
    adata,
    outdir=OUTDIR,
    min_cells_per_donor=MIN_CELLS_PER_DONOR,
    id_candidates=ID_CANDIDATES,
    disease_col=DISEASE_COL,
    categorical_cols=CATEGORICAL_COLS,
    numeric_cols=NUMERIC_COLS,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    obs = adata.obs.copy()

    # Pick subject-level ID column
    id_col = choose_existing_column(obs, id_candidates, required=True)
    if disease_col not in obs.columns:
        raise ValueError(f"'{disease_col}' not found in adata.obs")

    print(f"Using donor/subject ID column: {id_col}")
    print(f"Using disease column: {disease_col}")

    # Keep columns that actually exist
    categorical_cols = [c for c in categorical_cols if c in obs.columns]
    numeric_cols = [c for c in numeric_cols if c in obs.columns]

    # Standardize disease and ID as strings
    obs[id_col] = obs[id_col].astype(str)
    obs[disease_col] = obs[disease_col].astype(str)

    # --------------------------------
    # 1. Cell-level summaries
    # --------------------------------
    cells_per_donor = (
        obs.groupby(id_col)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    save_table(cells_per_donor, outdir / "cells_per_donor.csv")

    cells_per_disease = (
        obs.groupby(disease_col)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    save_table(cells_per_disease, outdir / "cells_per_disease.csv")

    # --------------------------------
    # 2. Donor-level metadata table
    # --------------------------------
    donor_meta = pd.DataFrame({id_col: sorted(obs[id_col].unique())})
    donor_disease = obs.groupby(id_col)[disease_col].agg(collapse_unique).reset_index()
    donor_disease_n = (
        obs.groupby(id_col, observed=False)[disease_col]
        .agg(lambda x: len(nonnull_unique(x)))
        .reset_index(name=f"{disease_col}__n_unique")
    )

    donor_meta = donor_meta.merge(donor_disease, on=id_col, how="left")
    donor_meta = donor_meta.merge(donor_disease_n, on=id_col, how="left")

    for col in categorical_cols:
        tmp_val = obs.groupby(id_col)[col].agg(collapse_unique).reset_index()
        tmp_n = (
            obs.groupby(id_col)[col]
            .agg(lambda x: len(nonnull_unique(x)))
            .reset_index(name=f"{col}__n_unique")
        )
        donor_meta = donor_meta.merge(tmp_val, on=id_col, how="left")
        donor_meta = donor_meta.merge(tmp_n, on=id_col, how="left")

    # add numeric donor-level summaries
    for col in numeric_cols:
        vals = extract_numeric(obs[col])
        tmp = (
            pd.DataFrame({id_col: obs[id_col].values, col: vals.values})
            .groupby(id_col)[col]
            .median()
            .reset_index()
            .rename(columns={col: f"{col}__median"})
        )
        donor_meta = donor_meta.merge(tmp, on=id_col, how="left")

    donor_meta = donor_meta.merge(cells_per_donor, on=id_col, how="left")
    donor_meta = donor_meta.sort_values(["n_cells", id_col], ascending=[False, True])

    save_table(donor_meta, outdir / "donor_level_metadata.csv")

    # --------------------------------
    # 3. Donor-level summaries by disease
    # --------------------------------
    donors_per_disease = (
        donor_meta.groupby(disease_col)[id_col]
        .nunique()
        .reset_index(name="n_donors")
        .sort_values("n_donors", ascending=False)
    )
    save_table(donors_per_disease, outdir / "donors_per_disease.csv")

    donor_cell_summary_by_disease = (
        donor_meta.groupby(disease_col)["n_cells"]
        .agg(
            n_donors="size",
            median_cells_per_donor="median",
            mean_cells_per_donor="mean",
            min_cells_per_donor="min",
            max_cells_per_donor="max",
        )
        .reset_index()
        .sort_values("n_donors", ascending=False)
    )
    save_table(donor_cell_summary_by_disease, outdir / "donor_cell_summary_by_disease.csv")

    # numeric summaries by disease
    numeric_summary_tables = {}
    for col in numeric_cols:
        med_col = f"{col}__median"
        if med_col not in donor_meta.columns:
            continue
        tmp = (
            donor_meta.groupby(disease_col)[med_col]
            .agg(["count", "median", "mean", "min", "max"])
            .reset_index()
        )
        numeric_summary_tables[col] = tmp
        save_table(tmp, outdir / f"{col}_summary_by_disease.csv")

    # --------------------------------
    # 4. Ambiguity / consistency checks
    # --------------------------------
    ambiguity_cols = [disease_col] + categorical_cols
    ambiguity_records = []

    for col in ambiguity_cols:
        nuniq_col = f"{col}__n_unique"
        if nuniq_col in donor_meta.columns:
            bad = donor_meta.loc[donor_meta[nuniq_col] > 1, [id_col, col, nuniq_col, "n_cells"]].copy()
            bad["column"] = col
            ambiguity_records.append(bad)

    if len(ambiguity_records) > 0:
        ambiguous_donors = pd.concat(ambiguity_records, axis=0, ignore_index=True)
    else:
        ambiguous_donors = pd.DataFrame(columns=[id_col, "column", "n_cells"])

    save_table(ambiguous_donors, outdir / "ambiguous_donor_metadata.csv")

    # donors with low cell counts
    low_cell_donors = donor_meta.loc[donor_meta["n_cells"] < min_cells_per_donor].copy()
    save_table(low_cell_donors, outdir / f"low_cell_donors_lt_{min_cells_per_donor}.csv")

    # --------------------------------
    # 5. Disease confounding cross-tabs
    # donor-level and cell-level
    # --------------------------------
    available_cat_cols = [c for c in categorical_cols if c in donor_meta.columns]

    for col in available_cat_cols:
        # donor-level crosstab
        donor_ct = pd.crosstab(donor_meta[disease_col], donor_meta[col], dropna=False)
        save_crosstab(donor_ct, outdir / f"donor_crosstab_{disease_col}_by_{col}.csv")
        if donor_ct.shape[0] > 0 and donor_ct.shape[1] > 0:
            plot_heatmap(
                donor_ct,
                title=f"Donors: {disease_col} x {col}",
                outpath=outdir / f"donor_heatmap_{disease_col}_by_{col}.png"
            )

        # cell-level crosstab
        cell_ct = pd.crosstab(obs[disease_col], obs[col], dropna=False)
        save_crosstab(cell_ct, outdir / f"cell_crosstab_{disease_col}_by_{col}.csv")
        if cell_ct.shape[0] > 0 and cell_ct.shape[1] > 0:
            plot_heatmap(
                cell_ct,
                title=f"Cells: {disease_col} x {col}",
                outpath=outdir / f"cell_heatmap_{disease_col}_by_{col}.png"
            )

    # --------------------------------
    # 6. Plots
    # --------------------------------
    # donors per disease
    tmp = donors_per_disease.set_index(disease_col)["n_donors"]
    plot_bar(
        tmp,
        title="Number of donors per disease",
        ylabel="Number of donors",
        outpath=outdir / "bar_donors_per_disease.png"
    )

    # cells per disease
    tmp = cells_per_disease.set_index(disease_col)["n_cells"]
    plot_bar(
        tmp,
        title="Number of mesothelial cells per disease",
        ylabel="Number of cells",
        outpath=outdir / "bar_cells_per_disease.png"
    )

    # cells per donor histogram
    plot_hist(
        donor_meta["n_cells"],
        title="Mesothelial cells per donor",
        xlabel="Number of mesothelial cells",
        outpath=outdir / "hist_cells_per_donor.png"
    )

    # boxplot of cells per donor by disease
    plt.figure(figsize=(9, 4.8))
    order = (
        donor_meta.groupby(disease_col)["n_cells"]
        .median()
        .sort_values(ascending=False)
        .index
    )
    if HAVE_SEABORN:
        sns.boxplot(data=donor_meta, x=disease_col, y="n_cells", order=order)
        sns.stripplot(data=donor_meta, x=disease_col, y="n_cells", order=order, color="black", alpha=0.6, size=4)
    else:
        donor_meta.boxplot(column="n_cells", by=disease_col, rot=45)
        plt.suptitle("")
    plt.xticks(rotation=45, ha="right")
    plt.title("Cells per donor by disease")
    plt.xlabel("")
    plt.ylabel("Mesothelial cells per donor")
    plt.tight_layout()
    plt.savefig(outdir / "box_cells_per_donor_by_disease.png", dpi=200)
    plt.close()

    # --------------------------------
    # 7. Summary report
    # --------------------------------
    report = {
        "n_cells_total": int(adata.n_obs),
        "n_genes_total": int(adata.n_vars),
        "id_col_used": id_col,
        "disease_col_used": disease_col,
        "n_unique_donors": int(donor_meta[id_col].nunique()),
        "n_unique_diseases": int(donor_meta[disease_col].nunique()),
        "min_cells_per_donor_threshold": int(min_cells_per_donor),
        "n_low_cell_donors": int((donor_meta["n_cells"] < min_cells_per_donor).sum()),
        "n_ambiguous_disease_donors": int((donor_meta[f"{disease_col}__n_unique"] > 1).sum()),
        "categorical_cols_reviewed": categorical_cols,
        "numeric_cols_reviewed": numeric_cols,
    }

    with open(outdir / "summary_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("\n=== Metadata review complete ===")
    print(json.dumps(report, indent=2))

    if report["n_ambiguous_disease_donors"] > 0:
        warnings.warn(
            f"{report['n_ambiguous_disease_donors']} donors map to multiple disease labels. "
            "Inspect ambiguous_donor_metadata.csv before disease-state analysis."
        )

    if report["n_low_cell_donors"] > 0:
        warnings.warn(
            f"{report['n_low_cell_donors']} donors have < {min_cells_per_donor} mesothelial cells. "
            "Consider excluding them from donor-level composition / pseudobulk analyses."
        )

    return {
        "obs": obs,
        "donor_meta": donor_meta,
        "cells_per_donor": cells_per_donor,
        "cells_per_disease": cells_per_disease,
        "donors_per_disease": donors_per_disease,
        "donor_cell_summary_by_disease": donor_cell_summary_by_disease,
        "ambiguous_donors": ambiguous_donors,
        "low_cell_donors": low_cell_donors,
        "report": report,
    }
