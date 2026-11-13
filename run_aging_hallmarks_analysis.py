import os
import json
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("results", exist_ok=True)


def preprocess_and_combat(h5ad_path: str):
    """Load raw h5ad, filter missingness, log2-transform, median-normalize, impute, and run ComBat."""
    print("1. Loading and preprocessing proteomics data...")
    adata = ad.read_h5ad(h5ad_path)
    # Remove single mislabeled sample Blue5_F12 ('Drug')
    adata = adata[adata.obs["drug"] != "Drug"].copy()

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.copy()
    zero_rate = (X == 0).mean(axis=0)
    keep_prots = zero_rate <= 0.50
    adata = adata[:, keep_prots].copy()
    X = X[:, keep_prots]
    print(f"   Retained {adata.n_vars} / {len(zero_rate)} proteins with <=50% zeros across {adata.n_obs} samples.")

    with np.errstate(divide="ignore"):
        X_log = np.where(X > 0, np.log2(X), np.nan)

    # Sample-wise median normalization
    sample_medians = np.nanmedian(X_log, axis=1, keepdims=True)
    global_median = np.nanmedian(sample_medians)
    X_norm = X_log - sample_medians + global_median

    # Perseus-style left-censored imputation
    rng = np.random.default_rng(42)
    X_imp = X_norm.copy()
    for i in range(X_imp.shape[0]):
        m = np.isnan(X_imp[i])
        if m.any():
            mu = np.nanmean(X_imp[i])
            sd = np.nanstd(X_imp[i])
            X_imp[i, m] = rng.normal(loc=mu - 1.8 * sd, scale=0.3 * sd, size=m.sum())

    adata_pre = adata.copy()
    adata_pre.X = X_imp.copy()

    adata_post = adata.copy()
    adata_post.X = X_imp.copy()
    sc.pp.combat(adata_post, key="batch")

    return adata_pre, adata_post, X_log


def plot_combat_qc(adata_pre, adata_post, X_log_raw):
    """Plot PCA of QC samples and batch-wise median distributions (boxplot+stripplot) before/after ComBat."""
    print("2. Generating ComBat QC figure...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A & B: PCA of QC samples before and after ComBat
    for ax, ad_obj, title in [
        (axes[0, 0], adata_pre, "QC Samples PCA — Before ComBat"),
        (axes[0, 1], adata_post, "QC Samples PCA — After ComBat"),
    ]:
        qc = ad_obj[ad_obj.obs["sample_type"] == "QC"].copy()
        sc.pp.pca(qc, n_comps=2)
        coords = qc.obsm["X_pca"]
        var_exp = qc.uns["pca"]["variance_ratio"] * 100
        df_pca = pd.DataFrame({
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
            "batch": qc.obs["batch"].values,
        })
        sns.scatterplot(data=df_pca, x="PC1", y="PC2", hue="batch", s=55, alpha=0.85, ax=ax, palette="tab20")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)

    # Panel C & D: Boxplot + stripplot of sample medians per batch (raw log2 vs post-ComBat)
    raw_meds = np.nanmedian(X_log_raw, axis=1)
    post_meds = np.median(adata_post.X, axis=1)
    df_med = pd.DataFrame({
        "batch": adata_post.obs["batch"].values,
        "Raw Log2 Median": raw_meds,
        "Post-ComBat Median": post_meds,
    })
    batch_order = sorted(df_med["batch"].unique())

    for ax, col, title in [
        (axes[1, 0], "Raw Log2 Median", "Sample Median Log2 Intensity — Raw"),
        (axes[1, 1], "Post-ComBat Median", "Sample Median Log2 Intensity — Normalized + ComBat"),
    ]:
        sns.boxplot(data=df_med, x="batch", y=col, order=batch_order, color="lightgray", fliersize=0, ax=ax)
        sns.stripplot(data=df_med, x="batch", y=col, order=batch_order, size=2.5, alpha=0.55, color="#1f77b4", ax=ax)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Batch")
        ax.set_ylabel("Median Log2 Intensity")
        ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig("results/fig1_combat_qc.png", dpi=250, bbox_inches="tight")
    plt.close()


def build_gene_sets(var_names):
    """Define the 12 Hallmarks of Aging gene sets (Lopez-Otin et al., 2023) + 2 aggregate aging sets."""
    print("3. Fetching and curating 14 aging gene sets...")
    hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
    reactome = gp.get_library(name="Reactome_Pathways_2024", organism="Human")
    gobp = gp.get_library(name="GO_Biological_Process_2023", organism="Human")

    gene_sets = {
        "1. Genomic Instability (DNA Repair)": hallmark["DNA Repair"],
        "2. Telomere Attrition (Telomere Ext.)": reactome["Extension of Telomeres"],
        "3. Epigenetic Alterations (Chromatin Org.)": reactome["Chromatin Organization"],
        "4. Loss of Proteostasis (UPR)": hallmark["Unfolded Protein Response"],
        "5. Disabled Macroautophagy": reactome["Macroautophagy"],
        "6. Nutrient-Sensing (mTORC1)": hallmark["mTORC1 Signaling"],
        "7. Mitochondrial Dysfunction (OXPHOS)": hallmark["Oxidative Phosphorylation"],
        "8. Cellular Senescence": reactome["Cellular Senescence"],
        "9. Stem Cell Exhaustion (Maint.)": gobp["Positive Regulation Of Stem Cell Population Maintenance (GO:1902459)"],
        "10. Intercellular Comm. (Cell-Cell)": reactome["Cell-Cell Communication"],
        "11. Chronic Inflammation": hallmark["Inflammatory Response"],
        "12. Dysbiosis (TLR Cascades)": reactome["Toll-like Receptor Cascades"],
    }

    for name in ["DEMAGALHAES_AGING_UP", "DEMAGALHAES_AGING_DN"]:
        with open(f"gene-sets/{name}.v2026.2.Hs.json") as f:
            d = json.load(f)
        gene_sets[name] = d[name]["geneSymbols"]

    prots = set(var_names)
    summary = []
    for k, genes in gene_sets.items():
        ov = len(set(genes) & prots)
        summary.append({"GeneSet": k, "Total_Genes": len(genes), "Proteomics_Overlap": ov})
    df_gs = pd.DataFrame(summary)
    df_gs.to_csv("results/geneset_summary.csv", index=False)
    with open("results/aging_hallmarks_genesets.json", "w") as f:
        json.dump(gene_sets, f, indent=2)
    print(df_gs.to_string(index=False))
    return gene_sets


def run_de_and_gsea(adata_post, gene_sets):
    """Run Welch's t-test DE (vs DMSO) and Pre-Ranked GSEA for all drugs."""
    print("4. Running Differential Expression (BH FDR) and Pre-Ranked GSEA for all treatments...")
    # Exclude QC samples
    ad_trt = adata_post[adata_post.obs["sample_type"] != "QC"].copy()
    ctrl_mask = ad_trt.obs["drug"] == "DMSO"
    ctrl_X = ad_trt[ctrl_mask].X
    prots = ad_trt.var_names

    drugs = sorted([d for d in ad_trt.obs["drug"].unique() if d != "DMSO"])
    print(f"   Testing {len(drugs)} treatments against DMSO (n={ctrl_mask.sum()})...")

    de_summaries = []
    gsea_records = []
    t0 = time.time()

    for idx, drug in enumerate(drugs):
        d_mask = ad_trt.obs["drug"] == drug
        d_X = ad_trt[d_mask].X
        log2fc = d_X.mean(axis=0) - ctrl_X.mean(axis=0)
        t_stat, p_val = stats.ttest_ind(d_X, ctrl_X, axis=0, equal_var=False)
        # Handle any nan t_stat
        t_stat = np.nan_to_num(t_stat, nan=0.0)
        p_val = np.nan_to_num(p_val, nan=1.0)
        _, q_val, _, _ = multipletests(p_val, method="fdr_bh")

        de_summaries.append({
            "drug": drug,
            "n_rep": int(d_mask.sum()),
            "batch": ",".join(sorted(ad_trt.obs.loc[d_mask, "batch"].unique())),
            "n_DE_q05": int((q_val < 0.05).sum()),
            "n_DE_q10": int((q_val < 0.10).sum()),
            "n_DE_q20": int((q_val < 0.20).sum()),
            "n_up_q05": int(((q_val < 0.05) & (log2fc > 0)).sum()),
            "n_dn_q05": int(((q_val < 0.05) & (log2fc < 0)).sum()),
        })

        rnk = pd.Series(t_stat, index=prots).sort_values(ascending=False)
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            min_size=5,
            max_size=500,
            permutation_num=1000,
            seed=42,
            threads=4,
            verbose=False,
        )
        res_df = pre_res.res2d.copy()
        res_df["drug"] = drug
        gsea_records.append(res_df)

        if (idx + 1) % 25 == 0 or (idx + 1) == len(drugs):
            print(f"   Processed {idx + 1}/{len(drugs)} drugs ({time.time() - t0:.1f}s)")

    df_de = pd.DataFrame(de_summaries)
    df_gsea = pd.concat(gsea_records, ignore_index=True)
    df_gsea["NES"] = pd.to_numeric(df_gsea["NES"], errors="coerce")
    df_gsea["FDR q-val"] = pd.to_numeric(df_gsea["FDR q-val"], errors="coerce")

    df_de.to_csv("results/de_summary_per_drug.csv", index=False)
    df_gsea.to_csv("results/gsea_all_drugs_14pathways.csv", index=False)
    return df_de, df_gsea


def merge_metadata_and_plot(df_de, df_gsea, excel_path, gene_sets):
    """Merge Excel drug metadata, filter to significant drugs, and plot NES heatmap + distributions."""
    print("5. Merging drug metadata and generating NES heatmap...")
    df_meta = pd.read_excel(excel_path)
    meta_sub = df_meta[["ProductName", "Category", "PathWay", "Target"]].drop_duplicates("ProductName")

    df_de_m = df_de.merge(meta_sub, left_on="drug", right_on="ProductName", how="left")
    df_de_m.loc[df_de_m["drug"] == "DFO", "Category"] = "Positive_Control"
    df_de_m["Category"] = df_de_m["Category"].fillna("Unannotated")
    df_de_m.to_csv("results/de_summary_with_metadata.csv", index=False)

    pathway_order = list(gene_sets.keys())
    nes_mat = df_gsea.pivot(index="drug", columns="Term", values="NES")[pathway_order]
    fdr_mat = df_gsea.pivot(index="drug", columns="Term", values="FDR q-val")[pathway_order]

    nes_mat.to_csv("results/nes_matrix_all_drugs.csv")
    fdr_mat.to_csv("results/fdr_q_matrix_all_drugs.csv")

    # Filter to drugs with at least one significant pathway change (FDR q-val < 0.05)
    sig_mask_05 = (fdr_mat < 0.05).any(axis=1)
    sig_drugs_05 = fdr_mat.index[sig_mask_05].tolist()
    print(f"   Drugs with >=1 pathway at FDR q < 0.05: {len(sig_drugs_05)} / {len(fdr_mat)}")

    sig_mask_10 = (fdr_mat < 0.10).any(axis=1)
    print(f"   Drugs with >=1 pathway at FDR q < 0.10: {sig_mask_10.sum()} / {len(fdr_mat)}")

    nes_sig = nes_mat.loc[sig_drugs_05].copy()
    fdr_sig = fdr_mat.loc[sig_drugs_05].copy()

    # Map categories for row colors
    cat_map = df_de_m.set_index("drug")["Category"].to_dict()
    row_cats = pd.Series([cat_map.get(d, "Unannotated") for d in nes_sig.index], index=nes_sig.index, name="Category")
    unique_cats = sorted(row_cats.unique())
    cat_palette = dict(zip(unique_cats, sns.color_palette("Set2", n_colors=len(unique_cats))))
    row_colors = row_cats.map(cat_palette)

    # Create annotation matrix (* q<0.05, ** q<0.01)
    annot_mat = fdr_sig.applymap(lambda q: "**" if q < 0.01 else ("*" if q < 0.05 else ""))

    # Clustermap of NES for significant drugs
    fig_h = max(10, len(sig_drugs_05) * 0.24)
    g = sns.clustermap(
        nes_sig,
        cmap="RdBu_r",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        col_cluster=False,
        row_cluster=True,
        row_colors=row_colors,
        annot=annot_mat,
        fmt="",
        annot_kws={"size": 7, "color": "black"},
        linewidths=0.3,
        linecolor="#e0e0e0",
        figsize=(15, fig_h),
        cbar_kws={"label": "Normalized Enrichment Score (NES)"},
        dendrogram_ratio=(0.12, 0.02),
        colors_ratio=0.025,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=40, ha="right", fontsize=10)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
    g.ax_heatmap.set_xlabel("Aging Hallmarks & Aggregate Gene Sets", fontsize=11, fontweight="bold")
    g.ax_heatmap.set_ylabel("Drug Treatments (FDR q < 0.05 in >=1 pathway)", fontsize=11, fontweight="bold")
    g.fig.suptitle(
        f"Aging Hallmarks Pre-Ranked GSEA NES Heatmap ({len(sig_drugs_05)} Significant Drugs, FDR q < 0.05)",
        y=1.01,
        fontsize=13,
        fontweight="bold",
    )

    # Add legend for Drug Category
    handles = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=cat_palette[c], markersize=9, label=c) for c in unique_cats]
    g.fig.legend(handles=handles, title="Drug Category", loc="upper left", bbox_to_anchor=(0.01, 0.99), frameon=True, fontsize=8, title_fontsize=9)

    plt.savefig("results/fig2_aging_hallmarks_nes_heatmap.png", dpi=250, bbox_inches="tight")
    plt.close()

    # Figure 3: Boxplot + Stripplot of NES across all drugs per pathway (No-Barplot Rule)
    fig, ax = plt.subplots(figsize=(14, 7))
    df_plot = df_gsea.copy()
    df_plot["Term"] = pd.Categorical(df_plot["Term"], categories=pathway_order, ordered=True)
    df_plot["Significant (q < 0.05)"] = df_plot["FDR q-val"] < 0.05

    sns.boxplot(data=df_plot, x="Term", y="NES", color="#f0f0f0", fliersize=0, ax=ax)
    sns.stripplot(
        data=df_plot[~df_plot["Significant (q < 0.05)"]],
        x="Term",
        y="NES",
        color="#999999",
        alpha=0.35,
        size=3,
        jitter=0.25,
        ax=ax,
        label="FDR q >= 0.05",
    )
    sns.stripplot(
        data=df_plot[df_plot["Significant (q < 0.05)"]],
        x="Term",
        y="NES",
        color="#d62728",
        alpha=0.8,
        size=4.5,
        jitter=0.25,
        ax=ax,
        label="FDR q < 0.05",
    )
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40, ha="right", fontsize=9.5)
    ax.set_xlabel("Aging Hallmarks & Aggregate Gene Sets", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=11, fontweight="bold")
    ax.set_title("Distribution of Drug-Induced NES Across 14 Aging Pathways (All 173 Treatments)", fontsize=12, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", frameon=True)
    plt.tight_layout()
    plt.savefig("results/fig3_pathway_nes_distributions.png", dpi=250, bbox_inches="tight")
    plt.close()

    # Save top hits table
    sig_long = df_gsea[df_gsea["FDR q-val"] < 0.05].merge(
        df_de_m[["drug", "Category", "Target", "PathWay", "n_DE_q05"]], on="drug", how="left"
    )
    sig_long = sig_long.sort_values(["FDR q-val", "NES"], ascending=[True, False])
    sig_long.to_csv("results/significant_drug_pathway_hits_q05.csv", index=False)
    print("   Saved all results and figures to results/.")


if __name__ == "__main__":
    adata_pre, adata_post, X_log_raw = preprocess_and_combat("data/raw (1).h5ad")
    plot_combat_qc(adata_pre, adata_post, X_log_raw)
    gene_sets = build_gene_sets(adata_post.var_names)
    df_de, df_gsea = run_de_and_gsea(adata_post, gene_sets)
    merge_metadata_and_plot(df_de, df_gsea, "data/drug_details (1).xlsx", gene_sets)
    print("DONE!")
