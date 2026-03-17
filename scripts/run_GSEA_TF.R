#############################################################
# GSEA + TF ACTIVITY ANALYSIS
#############################################################

# -------------------------------
# 1. Setup
# -------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(msigdbr)
  library(fgsea)
  library(decoupleR)
  library(dorothea)
  library(OmnipathR)
})

set.seed(42)

dir.create("../results/GSEA", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/TF_activity", recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 2. Load input data
# -------------------------------

# Differential expression results
res2_ann <- read.csv("ShinyCD4/data/dea_pseudobulk_condition_v5.csv")

# Adjust p-values
res2_ann$adj.P.Val <- p.adjust(res2_ann$p_value, method = "fdr")

# Ranking for GSEA (t-statistic)
ranking <- setNames(res2_ann$t_value, res2_ann$gene)


# -------------------------------
# 3. GSEA - REACTOME (C2)
# -------------------------------

c2_gene_sets <- msigdbr(species = "human", category = "C2")

react <- c2_gene_sets %>%
  filter(gs_subcollection == "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

gsea.c2 <- fgseaMultilevel(
  pathways = react,
  stats = ranking,
  minSize = 5,
  maxSize = 500
)

# Format output
gsea.c2_out <- gsea.c2 %>%
  mutate(
    leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
  )

# Save
write.csv(gsea.c2_out,
          "../results/GSEA/2026_SSc_CD4_GSEA_REACTOME.csv",
          row.names = FALSE)


# -------------------------------
# 4. GSEA - Hallmark (H)
# -------------------------------

H_gene_sets <- msigdbr(species = "human", category = "H")

H.ensembl.ls <- H_gene_sets %>%
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(all.genes = list(unique(gene_symbol))) %>%
  deframe()

gsea.H <- fgseaMultilevel(
  pathways = H.ensembl.ls,
  stats = ranking,
  minSize = 10,
  maxSize = 500
)

gsea.H_out <- gsea.H %>%
  mutate(
    leadingEdge = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
  )

# Save
write.csv(gsea.H_out,
          "../results/GSEA/2026_SSc_CD4_GSEA_HALLMARK.csv",
          row.names = FALSE)


# -------------------------------
# 5. TF Activity - per sample
# -------------------------------

# Load normalized pseudobulk counts
counts_con <- readRDS("ShinyCD4/data/2026_SSc_CD4_counts_norm.rds")

counts_df <- as.data.frame(as.matrix(counts_con))

# TF network
net <- decoupleR::get_collectri(
  organism = 'human',
  split_complexes = FALSE
)

# Run ULM (per sample TF activity)
sample_acts <- decoupleR::run_ulm(
  mat = counts_df,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 5
)

# Save
write.csv(sample_acts,
          "../results/TF_activity/2026_SSc_CD4_sample_TF_activity.csv",
          row.names = FALSE)


# -------------------------------
# 6. TF Activity - contrast level
# -------------------------------

# Prepare DEG matrix (t-values)
deg <- res2_ann %>%
  dplyr::select(gene, t_value) %>%
  dplyr::filter(!is.na(t_value)) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Run ULM on contrast
contrast_acts <- decoupleR::run_ulm(
  mat = deg[, 't_value', drop = FALSE],
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'mor',
  minsize = 5
)

# Save
write.csv(contrast_acts,
          "../results/TF_activity/2026_SSc_CD4_contrast_TF_activity.csv",
          row.names = FALSE)


# -------------------------------
# 7. Save ranking (optional, useful)
# -------------------------------

write.csv(
  data.frame(gene = names(ranking), t_value = ranking),
  "../results/GSEA/2026_SSc_CD4_ranking.csv",
  row.names = FALSE
)