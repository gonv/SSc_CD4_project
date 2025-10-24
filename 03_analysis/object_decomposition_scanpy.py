#!/usr/bin/env python3 

import scanpy as sc
import numpy as np
import scipy.sparse as sparse
import pandas as pd
from scipy import io

# ============================================================================
# 1. Load AnnData object
# ============================================================================
adata3 = sc.read(
    '/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/scRNAseq_SSc_CD4_allMerge_PCs_it92_2_5_CT_TCR.h5ad'
)
print("Loaded AnnData successfully")

# ============================================================================
# 2. Export raw counts as MatrixMarket (sparse coordinate)
# ============================================================================
print("Exporting raw counts...")

raw_X = adata3.raw.X
if not sparse.issparse(raw_X):
    raw_X = sparse.csr_matrix(raw_X)

io.mmwrite(
    "/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/ssc_CD4_counts_raw.mtx",
    raw_X
)

# ============================================================================
# 3. Export features (genes) and barcodes (cells)
# ============================================================================
print("Exporting features and barcodes...")

adata3.raw.var.to_csv(
    "/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/ssc_CD4_features.tsv",
    sep="\t"
)
adata3.obs.to_csv(
    "/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/ssc_CD4_barcodes.tsv",
    sep="\t"
)

# ============================================================================
# 4. Export graph connectivities and distances (already sparse)
# ============================================================================
print("Exporting features and barcodes...")

if "connectivities" in adata3.obsp:
    print("Exporting connectivities...")
    conn = adata3.obsp["connectivities"]
    if not sparse.issparse(conn):
        conn = sparse.csr_matrix(conn)
    io.mmwrite(
        "/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/ssc_CD4_neighbors_connectivities.mtx",
        conn
    )

if "distances" in adata3.obsp:
    print("Exporting distances...")
    dist = adata3.obsp["distances"]
    if not sparse.issparse(dist):
        dist = sparse.csr_matrix(dist)
    io.mmwrite(
        "/lab-share/IM-Gutierrez-e2/Public/gvillanueva/00_CD4/sc_ob/ssc_CD4_neighbors_distances.mtx",
        dist
    )

print("All exports completed successfully")
