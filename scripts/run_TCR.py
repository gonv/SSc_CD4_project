#!/usr/bin/env python3

import scanpy as sc
import scirpy as ir
import muon as mu
import os

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
os.chdir('/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob')

# ------------------------------------------------------------------
# 1. Load AnnData object
# ------------------------------------------------------------------
adata = sc.read('scRNAseq_SSc_CD4_allMerge_20PCs_TCR_it92_2.h5ad')

# ------------------------------------------------------------------
# 2. TCR preprocessing & QC
# ------------------------------------------------------------------
ir.io.upgrade_schema(adata)
ir.pp.index_chains(adata)
ir.tl.chain_qc(adata)

# ------------------------------------------------------------------
# 3. Keep only productive single-pair TCRαβ cells
# ------------------------------------------------------------------
adata = adata[
    (adata.obs["receptor_type"] == "TCR") &
    (adata.obs["chain_pairing"] == "single pair")
].copy()

print("Cells after TCR QC:", adata.n_obs)

# ------------------------------------------------------------------
# 4. Define clonotypes (exact AA match, standard definition)
# ------------------------------------------------------------------
N_JOBS = 24

ir.tl.define_clonotypes(
    adata,
    receptor_arms="all",   # usa VJ y VDJ (α y β)
    dual_ir="any",         # si quedara alguna célula con dual chains (idealmente no)
    n_jobs=N_JOBS
)

ir.pp.ir_dist(
    adata,
    metric="alignment",
    sequence="aa",
    cutoff=10,
    n_jobs=N_JOBS
)

ir.tl.define_clonotype_clusters(
    adata,
    sequence="aa",
    metric="alignment",
    receptor_arms="all",
    dual_ir="primary_only",
    key_added="cc_aa",
    n_jobs=N_JOBS
)

ir.tl.clonotype_network(adata, 
                        min_cells=10, 
                        sequence="aa",
                        random_state=42)

# Quick sanity check
print(
    adata.obs[["clone_id", "clone_id_size"]]
    .value_counts("clone_id_size")
    .head()
)

# ------------------------------------------------------------------
# 5. Export metadata for R
# ------------------------------------------------------------------
cols = [
    "sample",
    "condition",
    "cell_type",
    "clone_id",
    "clone_id_size",
    "cc_aa",
    "cc_aa_size"
]

adata.obs[cols].to_csv("TCR_metadata_for_R.csv")

# ------------------------------------------------------------------
# 6. Save AnnData
# ------------------------------------------------------------------
print("Initial cells:", sc.read('scRNAseq_SSc_CD4_allMerge_20PCs_TCR_it92_2.h5ad').n_obs)
print("Cells after TCR QC:", adata.n_obs)

adata.write('scRNAseq_SSc_CD4_allMerge_20PCs_TCR_clones_it92_2.h5ad')