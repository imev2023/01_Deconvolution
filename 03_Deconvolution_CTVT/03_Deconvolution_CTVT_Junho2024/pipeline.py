#!/usr/bin/env python3
"""
End-to-end single-cell RNA-seq processing pipeline for a folder of .h5ad files.

Implements the steps described by the user:
- Re-align genes to a GRCh38 (Cell Ranger 2020-A) reference gene list
- Per-dataset QC filtering: remove cells with <1000 UMIs, <500 detected genes, >7000 genes
- Doublet detection with Scrublet
- Geometric sketching to downsample while preserving diversity (geosketch)
- Merge datasets and perform batch correction with BBKNN
- Compute UMAP + Leiden
- Optional rule-based cell type annotation from a YAML marker list

Outputs:
- Per-dataset: filtered .h5ad, scrublet scores, optional figures
- Per-dataset: sketched .h5ad (subset used for downstream merging)
- Merged: batch-corrected .h5ad with UMAP/Leiden and optional annotations
- CSVs summarizing QC, doublets, clusters and (if provided) annotations

Usage example:
    python scrna_pipeline.py \
        --input-dir ./h5ad_inputs \
        --output-dir ./outputs \
        --ref-genes ./cellranger_2020A_genes.txt \
        --gene-field index \
        --sketch-n 20000 \
        --expected-doublet-rate 0.06 \
        --markers-yaml ./markers.yaml

# python scrna_pipeline.py \
#   --input-dir ./h5ad_inputs \
#   --output-dir ./outputs \
#   --ref-genes ./cellranger_2020A_genes.txt \
#   --gene-field index \
#   --sketch-n 20000 \
#   --expected-doublet-rate 0.06


Notes:
- The reference genes file should contain a single column of gene symbols (or Ensembl IDs), one per line,
  matching the chosen --gene-field.
- If your .h5ad stores symbols in adata.var['gene_symbols'] (or similar), set --gene-field gene_symbols (or the key name).
- Requires: scanpy>=1.8.2, anndata, bbknn, scrublet, geosketch, pandas, numpy, scipy, yaml, matplotlib

"""

import argparse
import os
import sys
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

# External tools
try:
    import scrublet as scr
except Exception as e:
    scr = None

try:
    from geosketch import gs
    _HAS_GEOSKETCH = True
except Exception:
    _HAS_GEOSKETCH = False

try:
    import bbknn  # noqa: F401 (import registers sc.external.pp.bbknn)
    _HAS_BBKNN = True
except Exception:
    _HAS_BBKNN = False

try:
    import yaml
    _HAS_YAML = True
except Exception:
    _HAS_YAML = False

sc.set_figure_params(figsize=(5, 5), dpi=120)


def parse_args():
    p = argparse.ArgumentParser(description="Single-cell RNA-seq pipeline (Scanpy + Scrublet + BBKNN)")
    p.add_argument("--input-dir", required=True, help="Folder containing .h5ad files")
    p.add_argument("--output-dir", required=True, help="Folder to write outputs")
    p.add_argument("--ref-genes", required=True, help="Text file: one gene per line (Cell Ranger 2020-A order)")
    p.add_argument("--gene-field", default="index", help="Which field holds gene IDs to align on: 'index' or a column in .var (e.g. 'gene_symbols' or 'gene_ids')")
    p.add_argument("--min-umis", type=int, default=1000, help="Cells with fewer UMIs are removed (empty droplets)")
    p.add_argument("--min-genes", type=int, default=500, help="Cells with fewer detected genes are removed (empty droplets)")
    p.add_argument("--max-genes", type=int, default=7000, help="Cells with more detected genes are removed (potential doublets)")
    p.add_argument("--expected-doublet-rate", type=float, default=0.06, help="Expected doublet rate for Scrublet (typical 0.05-0.1)")
    p.add_argument("--doublet-threshold", type=float, default=None, help="Optional fixed Scrublet score threshold. If omitted, Scrublet's automatic threshold is used")
    p.add_argument("--sketch-n", type=int, default=20000, help="Number of cells to keep per dataset via geometric sketch. If <=0, skip sketching")
    p.add_argument("--hvg-n", type=int, default=3000, help="Number of highly variable genes for PCA/sketching")
    p.add_argument("--random-seed", type=int, default=0, help="Random seed for reproducibility")
    p.add_argument("--markers-yaml", default=None, help="Optional YAML mapping of cell_type -> [marker genes]; used for provisional annotation")
    p.add_argument("--umap-min-dist", type=float, default=0.3, help="UMAP min_dist")
    p.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden resolution for clustering")
    return p.parse_args()


def read_ref_genes(path):
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                # split on whitespace or comma if needed
                s = s.split("\t")[0].split(",")[0]
                genes.append(s)
    # ensure unique while preserving order
    seen = set()
    ordered = []
    for g in genes:
        if g not in seen:
            ordered.append(g)
            seen.add(g)
    return ordered


def ensure_var_names(adata: ad.AnnData, gene_field: str):
    """Set adata.var_names to the desired gene field and make unique."""
    if gene_field == "index":
        pass
    else:
        if gene_field not in adata.var.columns:
            raise KeyError(f"gene_field '{gene_field}' not found in adata.var columns: {list(adata.var.columns)[:10]}...")
        adata.var_names = adata.var[gene_field].astype(str)
    adata.var_names_make_unique()
    return adata


def compute_shared_genes(h5ad_paths, ref_genes, gene_field):
    shared = set(ref_genes)
    for p in h5ad_paths:
        a = sc.read_h5ad(p, backed=None)
        ensure_var_names(a, gene_field)
        genes = set(a.var_names.tolist())
        shared &= genes
    # keep the Cell Ranger reference order
    shared_ordered = [g for g in ref_genes if g in shared]
    return shared_ordered


def qc_filter(adata: ad.AnnData, min_umis: int, min_genes: int, max_genes: int):
    # Basic QC metrics
    if "n_genes_by_counts" not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, log1p=False)
    # Filter
    keep = (
        (adata.obs["total_counts"] >= min_umis)
        & (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
    )
    adata._inplace_subset_obs(keep.values)
    return adata


def run_scrublet(adata: ad.AnnData, expected_doublet_rate: float, doublet_threshold: float | None):
    if scr is None:
        warnings.warn("Scrublet not installed; skipping doublet detection")
        adata.obs["doublet_score"] = np.nan
        adata.obs["predicted_doublet"] = False
        return adata
    # Use counts matrix for scrublet; ensure dense or csr
    counts = adata.X
    if sparse.issparse(counts):
        counts = counts.tocsr()
    else:
        counts = sparse.csr_matrix(counts)
    scrub = scr.Scrublet(counts, expected_doublet_rate=expected_doublet_rate, random_state=0)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=None, min_cells=None, synth_doublet_fraction=0.25)
    # Optionally override threshold
    if doublet_threshold is not None:
        predicted_doublets = doublet_scores >= doublet_threshold
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = pd.Series(predicted_doublets, index=adata.obs_names).astype(bool)
    # Remove predicted doublets
    adata = adata[~adata.obs["predicted_doublet"].values].copy()
    return adata


def geometric_sketch(adata: ad.AnnData, sketch_n: int, hvg_n: int, seed: int):
    if sketch_n is None or sketch_n <= 0 or adata.n_obs <= sketch_n:
        return adata.copy()
    # Normalize/log/HVG/PCA for sketch space
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=hvg_n, flavor="seurat_v3", subset=True)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=min(50, a.n_vars - 1), svd_solver="arpack", random_state=seed)

    if _HAS_GEOSKETCH:
        # Use PCA space for farthest-first traversal
        X = a.obsm["X_pca"]
        idx = gs(X, sketch_n, replace=False, seed=seed)
    else:
        warnings.warn("geosketch not installed; falling back to Scanpy's pp.subsample (random)")
        sc.pp.subsample(a, n_obs=sketch_n, random_state=seed)
        idx = a.obs_names
    # Select from the original (un-normalized) adata
    if isinstance(idx, (list, np.ndarray)) and not isinstance(idx, pd.Index):
        selected = adata[np.array(idx, dtype=int)].copy() if np.issubdtype(np.array(idx).dtype, np.integer) else adata[a.obs_names.isin(idx)].copy()
    else:
        selected = adata[adata.obs_names.isin(idx)].copy()
    selected.obs["sketched"] = True
    return selected


def merge_and_bbknn(adatas: list[ad.AnnData], batch_key: str, umap_min_dist: float, leiden_resolution: float, seed: int):
    merged = sc.concat(adatas, label=batch_key, keys=[a.uns.get("_dataset_name", f"ds{i}") for i, a in enumerate(adatas)], join="inner")

    # Standard pipeline
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    sc.pp.highly_variable_genes(merged, flavor="seurat_v3", n_top_genes=3000, subset=True)
    sc.pp.scale(merged, max_value=10)
    sc.tl.pca(merged, n_comps=min(50, merged.n_vars - 1), svd_solver="arpack", random_state=seed)

    if not _HAS_BBKNN:
        warnings.warn("bbknn not installed; using standard neighbors instead of BBKNN")
        sc.pp.neighbors(merged, n_pcs=min(50, merged.obsm["X_pca"].shape[1]))
    else:
        sc.external.pp.bbknn(merged, batch_key=batch_key, n_pcs=min(50, merged.obsm["X_pca"].shape[1]))

    sc.tl.umap(merged, min_dist=umap_min_dist, random_state=seed)
    sc.tl.leiden(merged, resolution=leiden_resolution, key_added="leiden")
    return merged


def provisional_annotation(adata: ad.AnnData, markers_yaml: str, batch_key: str = "batch"):
    if not markers_yaml:
        return None
    if not _HAS_YAML:
        warnings.warn("pyyaml not installed; skipping marker-based annotation")
        return None
    with open(markers_yaml, "r", encoding="utf-8") as f:
        marker_map = yaml.safe_load(f)
    if not isinstance(marker_map, dict):
        warnings.warn("markers YAML must be a mapping of cell_type -> [markers]")
        return None

    # Score each cell type and assign per cluster the best-scoring label
    scores = {}
    for cell_type, genes in marker_map.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) == 0:
            warnings.warn(f"No markers present for {cell_type}; skipping")
            continue
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=f"score_{cell_type}")
        scores[cell_type] = adata.obs[f"score_{cell_type}"]

    if not scores:
        return None

    # Aggregate by cluster
    clust = adata.obs["leiden"].astype(str)
    score_df = pd.DataFrame(scores)
    score_df["cluster"] = clust.values
    cluster_scores = score_df.groupby("cluster").mean()
    anno = cluster_scores.idxmax(axis=1).rename("provisional_annotation")

    # Write back to adata.obs per cell
    adata.obs = adata.obs.join(anno, on="leiden")

    # Store summary
    adata.uns["provisional_annotation"] = {
        "markers": marker_map,
        "cluster_scores": cluster_scores.to_dict(orient="index"),
    }
    return anno


def main():
    args = parse_args()
    np.random.seed(args.random_seed)

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ads = sorted([p for p in in_dir.glob("*.h5ad")])
    if len(h5ads) == 0:
        print(f"No .h5ad files found in {in_dir}", file=sys.stderr)
        sys.exit(1)

    ref_genes = read_ref_genes(args.ref_genes)
    print(f"Loaded {len(ref_genes)} reference genes")

    # Compute the shared gene list across datasets intersected with reference
    shared_genes = compute_shared_genes(h5ads, ref_genes, args.gene_field)
    if len(shared_genes) == 0:
        print("No shared genes across datasets after intersecting with reference. Check --gene-field and ref list.", file=sys.stderr)
        sys.exit(2)
    print(f"Shared genes across datasets (∩ ref): {len(shared_genes)}")

    per_dataset_sketched = []
    qc_summary_rows = []

    for path in h5ads:
        name = path.stem
        print(f"\n=== Processing: {name} ===")
        a = sc.read_h5ad(path)
        a.uns["_dataset_name"] = name
        ensure_var_names(a, args.gene_field)

        # Align to shared genes (consistent order from reference)
        order = [g for g in shared_genes if g in a.var_names]
        a = a[:, order].copy()

        # QC filter
        pre_n = a.n_obs
        sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None, log1p=False)
        a = qc_filter(a, args.min_umis, args.min_genes, args.max_genes)
        post_qc_n = a.n_obs

        # Doublet detection & removal
        a = run_scrublet(a, args.expected_doublet_rate, args.doublet_threshold)
        post_dd_n = a.n_obs

        # Save filtered dataset with scrublet scores
        ds_out_dir = out_dir / name
        ds_out_dir.mkdir(exist_ok=True, parents=True)
        a.write(ds_out_dir / f"{name}.qc_dd_filtered.h5ad")

        # Geometric sketching
        sketched = geometric_sketch(a, args.sketch_n, args.hvg_n, args.random_seed)
        sketched.obs["batch"] = name
        sketched.uns["_dataset_name"] = name
        sketched.write(ds_out_dir / f"{name}.sketched.h5ad")
        per_dataset_sketched.append(sketched)

        # QC summary row
        qc_summary_rows.append({
            "dataset": name,
            "n_cells_pre": pre_n,
            "n_cells_post_qc": post_qc_n,
            "n_cells_post_doublet": post_dd_n,
            "n_cells_sketched": sketched.n_obs,
            "n_genes": sketched.n_vars,
        })

    # Write QC summary
    qc_df = pd.DataFrame(qc_summary_rows)
    qc_df.to_csv(out_dir / "qc_summary.csv", index=False)

    # Merge + BBKNN + UMAP/Leiden
    print("\n=== Merging sketched datasets and running BBKNN/UMAP/Leiden ===")
    merged = merge_and_bbknn(per_dataset_sketched, batch_key="batch", umap_min_dist=args.umap_min_dist, leiden_resolution=args.leiden_resolution, seed=args.random_seed)

    # Optional marker-based provisional annotations
    if args.markers_yaml:
        print("Applying provisional annotation from marker YAML...")
        provisional_annotation(merged, args.markers_yaml, batch_key="batch")

    # Rank genes per cluster (for manual annotation support)
    sc.tl.rank_genes_groups(merged, groupby="leiden", method="wilcoxon")

    # Save merged object
    merged.write(out_dir / "merged_bbknn_umap.h5ad")

    # Export UMAP and cluster assignments
    emb = pd.DataFrame(merged.obsm["X_umap"], index=merged.obs_names, columns=["umap1", "umap2"])\
        .join(merged.obs[["batch", "leiden"]])
    if "provisional_annotation" in merged.obs.columns:
        emb = emb.join(merged.obs[["provisional_annotation"]])
    emb.to_csv(out_dir / "umap_leiden_annotations.csv")

    # Export ranked genes per cluster
    rg = merged.uns["rank_genes_groups"]
    groups = rg["names"].dtype.names
    rgg_rows = []
    for g in groups:
        names = rg["names"][g]
        scores = rg["scores"][g]
        pvals = rg["pvals_adj"][g]
        for rank, (gene, score, p) in enumerate(zip(names, scores, pvals), start=1):
            rgg_rows.append({"cluster": g, "rank": rank, "gene": gene, "score": score, "pval_adj": p})
    pd.DataFrame(rgg_rows).to_csv(out_dir / "rank_genes_groups.csv", index=False)

    # Figures (optional small set)
    print("Saving static figures...")
    sc.pl.umap(merged, color=["batch", "leiden"], wspace=0.4, show=False, save=f"_batch_leiden.png")
    if "provisional_annotation" in merged.obs.columns:
        sc.pl.umap(merged, color=["provisional_annotation"], show=False, save=f"_annotations.png")

    print("\nDone. Key outputs:")
    print(f" - QC summary: {out_dir / 'qc_summary.csv'}")
    print(f" - Merged AnnData: {out_dir / 'merged_bbknn_umap.h5ad'}")
    print(f" - UMAP & clusters: {out_dir / 'umap_leiden_annotations.csv'}")
    print(f" - Ranked markers: {out_dir / 'rank_genes_groups.csv'}")
    print(f" - Per-dataset outputs under: {out_dir}/<dataset>/ ...")


if __name__ == "__main__":
    main()
