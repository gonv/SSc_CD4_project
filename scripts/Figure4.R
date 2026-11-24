#!/usr/bin/env Rscript

# Figure 4 - R panels A-D and donor-level TCR analyses.
# Panels A and B use current annotations and validated cell-level clonotypes.
# Panel D uses the original cell-level antigen annotations across the full cohort.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

pathTO <- "figures/plots4figures/Figure4"
dir.create(pathTO, recursive = TRUE, showWarnings = FALSE)
custom_colors <- c(CTRL = "#D8BFD8", SSc = "#C3E5C3")
analysis_seed <- 42L
set.seed(42)

annotations <- read.csv("../ShinyCD4/data/20260124_SSc_CD4_meta_v3.csv",
                        row.names = 1, colClasses = "character",
                        na.strings = c("", "NA", "nan"))
meta <- read.csv("../data/SSc_CD4_meta_TCR.csv", colClasses = "character",
                 na.strings = c("", "NA", "nan"))
annotation_fields <- c("cell", "sample", "condition", "R09", "cell_type")
tcr_fields <- c("cell_id", "sample", "condition", "clone_id",
                "clone_id_size", "R09", "cell_type", "chain_pairing")

idx <- match(meta$cell_id, annotations$cell)

# Export annotation changes before using the latest labels.
annotation_changes <- data.frame(
  cell_id = meta$cell_id, donor = meta$sample,
  old_cluster = meta$R09, new_cluster = annotations$R09[idx],
  old_subset = meta$cell_type, new_subset = annotations$cell_type[idx]
)

write.csv(annotation_changes, file.path(pathTO, "annotation_comparison.csv"), row.names = FALSE)
meta$R09 <- annotations$R09[idx]
meta$cell_type <- annotations$cell_type[idx]

# Use the stable biological donor identifier for all aggregations.
meta$Sample <- meta$sample

# A clonotype may span several compartments. Classify it in the full donor
# repertoire first and count each observed cell once in its own compartment.
clone_key <- interaction(meta$Sample, meta$clone_id, drop = TRUE)
meta$size <- as.integer(ave(rep(1L, nrow(meta)), clone_key, FUN = sum))

meta$expanded <- as.integer(meta$size >= 2L)
meta$total <- 1L
meta$compartment <- "All"

summarize_counts <- function(data, grouping) {
  aggregate(data[c("expanded", "total")], data[grouping], sum)
}

analyze_resolution <- function(column, resolution) {
  data <- meta
  data$compartment <- data[[column]]
  strata <- c("Sample", "condition", "compartment")
  donors <- summarize_counts(data, strata)
  donors$fraction_expanded <- donors$expanded / donors$total
  donors$resolution <- resolution
  
  # These pooled summaries describe the figure and weight donors by cell yield.
  pooled <- summarize_counts(data, c("condition", "compartment"))
  pooled$fraction_expanded <- pooled$expanded / pooled$total
  pooled$share_of_expanded <- pooled$expanded /
    ave(pooled$expanded, pooled$condition, FUN = sum)
  pooled$resolution <- resolution
  
  # Count each expanded donor-clonotype once per occupied compartment.
  # A clone spanning compartments contributes one incidence to each.
  incidences <- unique(data[data$expanded == 1L, c(strata, "clone_id")])
  donors$expanded_clone_incidences <- 0L
  if (nrow(incidences)) {
    incidences$expanded_clone_incidences <- 1L
    counts <- aggregate(incidences["expanded_clone_incidences"], incidences[strata], sum)
    donors$expanded_clone_incidences <- NULL
    donors <- merge(donors, counts, by = strata, all.x = TRUE)
    donors$expanded_clone_incidences[is.na(donors$expanded_clone_incidences)] <- 0L
  }
  incidence_total <- ave(donors$expanded_clone_incidences, donors$Sample, FUN = sum)
  donors$share_of_expanded_clone_incidences <- ifelse(
    incidence_total > 0L, donors$expanded_clone_incidences / incidence_total, NA_real_
  )
  list(donors = donors, pooled = pooled)
}

# Compare unweighted donor means using a two-sided donor-label permutation
# test. Enumerate all assignments for this cohort, or sample for larger cohorts.
# A donor with zero observed cells in a compartment has no defined fraction.
compare_donors <- function(data, metric) {
  data <- data[is.finite(data[[metric]]), ]
  x <- data[[metric]]
  case <- data$condition == "SSc"
  n_case <- sum(case)
  n_ctrl <- sum(!case)
  result <- data.frame(n_SSc = n_case, n_CTRL = n_ctrl,
                       mean_SSc = if (n_case) mean(x[case]) else NA_real_,
                       mean_CTRL = if (n_ctrl) mean(x[!case]) else NA_real_)
  result$difference_pp <- 100 * (result$mean_SSc - result$mean_CTRL)
  result$p_value <- result$leave_one_out_min_pp <- result$leave_one_out_max_pp <- NA_real_
  if (min(n_case, n_ctrl) < 2L) return(result)
  
  permutation_difference <- function(indices) {
    case_sum <- sum(x[indices])
    case_sum / n_case - (sum(x) - case_sum) / n_ctrl
  }
  cutoff <- abs(result$difference_pp / 100) - 1e-12
  if (choose(length(x), n_case) <= 100000) {
    null <- combn(seq_along(x), n_case, permutation_difference)
    result$p_value <- mean(abs(null) >= cutoff)
  } else {
    null <- replicate(100000, permutation_difference(sample.int(length(x), n_case)))
    result$p_value <- (1 + sum(abs(null) >= cutoff)) / (length(null) + 1)
  }
  
  leave_one_out <- vapply(seq_along(x), function(i) {
    100 * (mean(x[-i][case[-i]]) - mean(x[-i][!case[-i]]))
  }, numeric(1))
  result$leave_one_out_min_pp <- min(leave_one_out)
  result$leave_one_out_max_pp <- max(leave_one_out)
  result
}

set.seed(analysis_seed)
analyses <- Map(analyze_resolution, c("compartment", "R09", "cell_type"),
                c("global", "cluster", "subset"))
donors <- do.call(rbind, lapply(analyses, "[[", "donors"))
pooled <- do.call(rbind, lapply(analyses, "[[", "pooled"))
contrasts <- function(data, metric) {
  strata <- split(data, list(data$resolution, data$compartment), drop = TRUE)
  results <- lapply(strata, function(group) {
    cbind(group[1, c("resolution", "compartment")], compare_donors(group, metric))
  })
  results <- do.call(rbind, results)
  # Adjust within each resolution for this endpoint, using full precision.
  results$p_BH <- ave(results$p_value, results$resolution,
                      FUN = function(p) p.adjust(p, method = "BH", n = length(p)))
  results
}
tests <- contrasts(donors, "fraction_expanded")
composition_tests <- contrasts(subset(donors, resolution == "subset"),
                               "share_of_expanded_clone_incidences")

write.csv(donors, file.path(pathTO, "donor_expansion.csv"), row.names = FALSE)
write.csv(pooled, file.path(pathTO, "pooled_expansion.csv"), row.names = FALSE)
write.csv(tests, file.path(pathTO, "donor_expansion_tests.csv"), row.names = FALSE)
write.csv(composition_tests, file.path(pathTO, "donor_clone_composition_tests.csv"),
          row.names = FALSE)


# ----------------------------------------------------------------------------------------
# Figure 4A - global expansion and expansion across compartments -------------------------
# ----------------------------------------------------------------------------------------

# Each TCR-assigned cell contributes once. These descriptive pooled
# fractions differ from the equal-weight donor means used in the tests.
expansion_global <- pooled %>%
  filter(resolution == "global") %>%
  mutate(Singleton = total - expanded, Expanded = expanded) %>%
  select(condition, Singleton, Expanded) %>%
  pivot_longer(c("Singleton", "Expanded"),
               names_to = "status", values_to = "cells") %>%
  group_by(condition) %>%
  mutate(fraction = cells / sum(cells),
         label = sprintf("%.2f", 100 * fraction)) %>%
  ungroup()
expansion_global$condition <- factor(expansion_global$condition, c("SSc", "CTRL"))

expansion_colors <- c(Expanded = "#D23931", Singleton = "#F7AF3C")
global_expansion_plot <- ggplot(expansion_global,
                                aes(x = 1, y = fraction, fill = status)) +
  geom_col(width = 1, colour = "white") +
  geom_label(aes(label = label), position = position_stack(vjust = 0.5),
             size = 2.5, fill = "white", colour = "black", show.legend = FALSE) +
  coord_polar(theta = "y") +
  facet_wrap(~condition, ncol = 1) +
  scale_fill_manual(values = expansion_colors, guide = "none") +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

cluster_levels <- as.character(sort(unique(as.integer(meta$R09))))
cellt <- c("Tn", "Tcm", "Tem", "Temra", "Treg")

plot_pooled_expansion <- function(resolution_name, compartment_order, title) {
  data <- donors %>%
    filter(resolution == resolution_name) %>%
    group_by(compartment) %>%
    summarise(Expanded = sum(expanded), Singleton = sum(total - expanded),
              .groups = "drop") %>%
    pivot_longer(c("Singleton", "Expanded"),
                 names_to = "status", values_to = "cells") %>%
    group_by(compartment) %>%
    mutate(fraction = cells / sum(cells)) %>%
    ungroup()
  data$compartment <- factor(data$compartment, compartment_order)
  ggplot(data, aes(compartment, fraction, fill = status)) +
    geom_col(width = 0.95, colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = expansion_colors) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = title, x = NULL, y = "Fraction of TCR-assigned cells", fill = NULL) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 9),
          legend.position = "bottom")
}

cluster_expansion_plot <- plot_pooled_expansion(
  "cluster", cluster_levels, "Clonotype expansion by cluster"
)
subset_expansion_plot <- plot_pooled_expansion(
  "subset", cellt, "Clonotype expansion by cell type"
)

figure4A <- global_expansion_plot + cluster_expansion_plot + subset_expansion_plot +
  plot_layout(widths = c(0.65, 2, 1.25), guides = "collect") &
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Figure 4B - descriptive clone sizes by cell type ----------------------------
# -----------------------------------------------------------------------------

# Use one point per donor-clonotype-subset incidence. A clonotype found in
# several subsets appears once in each occupied subset. Its plotted size
# is its full donor repertoire size. This avoids assigning every cell of
# a shared clonotype to its first observed subset.
clone_sizes <- meta %>%
  group_by(sample, condition, clone_id, cell_type) %>%
  summarise(clone_size = first(size), cells_in_subset = n(), .groups = "drop") %>%
  filter(clone_size >= 2)
clone_sizes$cell_type <- factor(clone_sizes$cell_type, cellt)
clone_sizes$condition <- factor(clone_sizes$condition, c("CTRL", "SSc"))

figure4B <- ggplot(clone_sizes, aes(cell_type, clone_size, fill = condition)) +
  geom_boxplot(width = 0.65, alpha = 0.6, outlier.shape = NA) +
  geom_point(shape = 21, size = 0.8, alpha = 0.5,
             position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.75, seed = analysis_seed)) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Clone size", x = NULL, y = "Full donor clonotype size", fill = NULL) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
write.csv(clone_sizes, file.path(pathTO, "clone_sizes_by_subset.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Figure 4C - donor-level primary-chain V and J gene usage
# ------------------------------------------------------------

# MERGED UPDATE C BEGIN - integrated gene usage analysis.
usage_pathTO <- "figures/plots4figures/Figure4C_gene_usage"
dir.create(usage_pathTO, recursive = TRUE, showWarnings = FALSE)

gene_columns <- c(TRAV = "IR_VJ_1_v_call", TRAJ = "IR_VJ_1_j_call",
                  TRBV = "IR_VDJ_1_v_call", TRBJ = "IR_VDJ_1_j_call")
n_bootstrap <- 10000L
n_display <- 8L

# Reuse the inputs loaded for Panels A and B.
tcr <- meta
vdj <- read.csv("../ShinyCD4/20260202_SSc_CD4_TCR_irVDJ.csv",
                row.names = 1, stringsAsFactors = FALSE, colClasses = "character",
                na.strings = c("", "NA", "nan", "None", "selfdefinedNA"))

meta_index <- match(tcr$cell_id, annotations$cell)
vdj_index <- match(tcr$cell_id, vdj$X)

# Clone sizes were already recomputed for the expansion analysis.
cells <- data.frame(donor = tcr$sample, condition = tcr$condition,
                    clone_id = tcr$clone_id, vdj[vdj_index, unname(gene_columns)],
                    check.names = FALSE)
long <- cells %>%
  pivot_longer(all_of(unname(gene_columns)), names_to = "column", values_to = "gene") %>%
  mutate(segment = names(gene_columns)[match(column, gene_columns)],
         gene = ifelse(is.na(gene) | trimws(gene) == "", NA_character_, gene))

# Check that a donor-clonotype has one consistent primary gene call per segment.
# Stop rather than selecting the first call if conflicting annotations occur.
discordance <- long %>%
  filter(!is.na(gene)) %>%
  distinct(donor, clone_id, segment, gene) %>%
  count(donor, clone_id, segment) %>%
  filter(n > 1)
if (nrow(discordance)) stop("Conflicting primary calls within donor-clonotypes.")

# ------------------------------------------------------------
# Coverage and donor frequencies  ----------------------------
# ------------------------------------------------------------

coverage <- long %>%
  group_by(donor, condition, segment) %>%
  summarise(total_cells = n(), assigned_cells = sum(!is.na(gene)),
            missing_cells = sum(is.na(gene)), .groups = "drop")
write.csv(coverage, file.path(usage_pathTO, "gene_call_coverage.csv"), row.names = FALSE)

# The denominator is valid calls in the segment, separately for each donor.
# Include zero frequencies for genes absent from a donor.
frequencies <- function(data, weighting) {
  denominators <- data %>%
    group_by(donor, condition, segment) %>%
    summarise(denominator = n(), .groups = "drop")
  counts <- data %>% count(donor, condition, segment, gene, name = "count")
  gene_set <- data %>% distinct(segment, gene)
  grid <- merge(unique(data[c("donor", "condition")]), gene_set, by = NULL)
  grid %>%
    left_join(counts, by = c("donor", "condition", "segment", "gene")) %>%
    left_join(denominators, by = c("donor", "condition", "segment")) %>%
    mutate(count = replace_na(count, 0L),
           fraction = count / denominator,
           weighting = weighting)
}
valid <- long %>% filter(!is.na(gene))
cell_frequencies <- frequencies(valid, "Cell")
clone_frequencies <- frequencies(
  valid %>% distinct(donor, condition, clone_id, segment, gene), "Clonotype"
)
donor_frequencies <- bind_rows(cell_frequencies, clone_frequencies)
write.csv(donor_frequencies, file.path(usage_pathTO, "donor_gene_frequencies.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------------------
# Exact donor-label permutation tests and bootstrap intervals ----------------------------
# ----------------------------------------------------------------------------------------

usage_donors <- tcr %>% distinct(sample, condition) %>% arrange(condition, sample)
n_case <- sum(usage_donors$condition == "SSc")
n_ctrl <- sum(usage_donors$condition == "CTRL")

assignments <- combn(seq_len(nrow(usage_donors)), n_case)

permutation_weights <- matrix(-1 / n_ctrl, nrow(usage_donors), ncol(assignments))
permutation_weights[cbind(as.vector(assignments),
                          rep(seq_len(ncol(assignments)), each = n_case))] <- 1 / n_case
observed_weights <- ifelse(usage_donors$condition == "SSc", 1 / n_case, -1 / n_ctrl)

# Resample donors independently within each condition. All genes share the
# same donor resamples. Intervals are pointwise 95% percentile intervals.
set.seed(analysis_seed)
case_index <- which(usage_donors$condition == "SSc")
ctrl_index <- which(usage_donors$condition == "CTRL")

bootstrap_weights <- replicate(n_bootstrap, {
  weights <- numeric(nrow(usage_donors))
  weights <- weights + tabulate(sample(case_index, n_case, replace = TRUE),
                                nbins = nrow(usage_donors)) / n_case
  weights <- weights - tabulate(sample(ctrl_index, n_ctrl, replace = TRUE),
                                nbins = nrow(usage_donors)) / n_ctrl
  weights
})

run_usage_tests <- function(data) {
  genes <- data %>% distinct(segment, gene) %>% arrange(segment, gene)
  matrix <- matrix(NA_real_, nrow(genes), nrow(usage_donors))
  for (i in seq_len(nrow(genes))) {
    values <- data %>% filter(segment == genes$segment[i], gene == genes$gene[i])
    matrix[i, ] <- values$fraction[match(usage_donors$sample, values$donor)]
  }
  if (any(!is.finite(matrix))) stop("Missing donor-segment denominators.")
  effect <- as.vector(matrix %*% observed_weights)
  null <- matrix %*% permutation_weights
  p_value <- rowMeans(abs(null) >= abs(effect) - 1e-12)
  bootstrap <- 100 * matrix %*% bootstrap_weights
  intervals <- t(apply(bootstrap, 1, quantile, probs = c(0.025, 0.975)))
  genes %>%
    mutate(
      n_SSc = n_case, n_CTRL = n_ctrl,
      mean_SSc = rowMeans(matrix[, case_index, drop = FALSE]),
      mean_CTRL = rowMeans(matrix[, ctrl_index, drop = FALSE]),
      difference_pp = 100 * effect,
      lower_95_pp = intervals[, 1], upper_95_pp = intervals[, 2],
      p_value = p_value,
      # All observed genes across the four segments form one testing family.
      p_BH = p.adjust(p_value, method = "BH"),
      weighting = first(data$weighting)
    )
}
gene_usage_tests <- bind_rows(run_usage_tests(cell_frequencies),
                              run_usage_tests(clone_frequencies))
write.csv(gene_usage_tests, file.path(usage_pathTO, "gene_usage_tests.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------------------------
# Figure 4C - differences in the eight most-used genes per segment ----------------------------
# ---------------------------------------------------------------------------------------------

# Select displayed genes by mean donor frequency across both conditions,
# without using effect size or significance. Test all observed genes above.
display_genes <- cell_frequencies %>%
  group_by(segment, gene) %>%
  summarise(mean_frequency = mean(fraction), .groups = "drop") %>%
  arrange(segment, desc(mean_frequency), gene) %>%
  group_by(segment) %>%
  slice_head(n = n_display) %>%
  ungroup()

plot_data <- gene_usage_tests %>% inner_join(display_genes, by = c("segment", "gene"))
plot_data$segment <- factor(plot_data$segment, c("TRAV", "TRAJ", "TRBV", "TRBJ"))
plot_data$gene <- factor(plot_data$gene, levels = rev(display_genes$gene))

write.csv(display_genes, file.path(usage_pathTO, "displayed_genes.csv"), row.names = FALSE)

plot_usage <- function(data, title) {
  ggplot(data, aes(difference_pp, gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_segment(aes(x = lower_95_pp, xend = upper_95_pp, yend = gene),
                 colour = "#3B747E", linewidth = 0.5) +
    geom_point(colour = "#174D59", size = 2) +
    facet_wrap(~segment, nrow = 2, scales = "free_y") +
    labs(title = title, x = "SSc - CTRL difference in gene usage (percentage points)",
         y = NULL,
         caption = "Equal-weight donor means; pointwise 95% donor-bootstrap intervals.") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 12, face = "bold"))
}

figure4C <- plot_usage(filter(plot_data, weighting == "Cell"), "Primary-chain V and J gene usage")
sensitivity_plot <- plot_usage(filter(plot_data, weighting == "Clonotype"),
                               "V and J usage with equal weight per clonotype")
# MERGED UPDATE C END

# ----------------------------------------------------------------------------------------
# Donor-level expansion plot supporting the cohort comparisons----------------------------
# ----------------------------------------------------------------------------------------

donor_data <- donors %>% filter(resolution == "cluster")
donor_data$compartment <- factor(donor_data$compartment, cluster_levels)
donor_data$condition <- factor(donor_data$condition, c("CTRL", "SSc"))
donor_expansion_plot <- ggplot(donor_data,
                               aes(compartment, fraction_expanded, fill = condition)) +
  geom_point(shape = 21, colour = "grey25",
             position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.65, seed = analysis_seed)) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 8,
               position = position_dodge(width = 0.65)) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Cluster", y = "Fraction of TCR-assigned cells", fill = NULL) +
  theme_classic()

# ------------------------------------------------------------
# Figure 4D - antigen annotations ----------------------------
# ------------------------------------------------------------

# REVIEW UPDATE D BEGIN - local antigen export and explicit denominators.
antigen_meta <- read.csv("../data/scRNAseq_SScCD4_adataFnl_it92.csv",
                         stringsAsFactors = FALSE,
                         colClasses = "character",
                         na.strings = c("", "NA", "nan"))
antigen_fields <- c("cell_id", "sample", "condition", "antigen.species")

# Match against the latest metadata using cell IDs. The antigen export
# contains the full cohort, including cells outside the single-pair TCR set.
antigen_idx <- match(antigen_meta$cell_id, annotations$cell)

adataEpi <- antigen_meta %>%
  filter(!is.na(antigen.species), trimws(antigen.species) != "")


# Preserve the original category labels, including ambiguous annotations.
# The denominator is the number of annotated cells in a category across
# both cohorts. This is a descriptive pooled plot, weighted by cell yield.
epitope_counts <- adataEpi %>%
  count(antigen.species, condition, name = "cells") %>%
  complete(antigen.species, condition = c("CTRL", "SSc"), fill = list(cells = 0L)) %>%
  pivot_wider(names_from = condition, values_from = cells) %>%
  mutate(
    total = CTRL + SSc,
    pctCTRL = CTRL / total,
    pctSSc = SSc / total,
    shown_in_figure = CTRL > 0 & SSc > 0
  )

# The original Rmd displayed categories observed in both cohorts. Retain
# that selection for Figure 4D and export all categories to the count table.
epitope_plot_data <- epitope_counts %>%
  filter(shown_in_figure) %>%
  select(antigen.species, CTRL, SSc) %>%
  pivot_longer(c("CTRL", "SSc"), names_to = "condition", values_to = "cells") %>%
  group_by(antigen.species) %>%
  mutate(fraction = cells / sum(cells)) %>%
  ungroup()

epitope_plot_data$antigen.species <- factor(
  epitope_plot_data$antigen.species,
  levels = sort(unique(epitope_plot_data$antigen.species))
)
epitope_plot_data$condition <- factor(epitope_plot_data$condition, c("CTRL", "SSc"))

figure4D <- ggplot(epitope_plot_data,
                   aes(antigen.species, fraction, fill = condition)) +
  geom_col(width = 0.95, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::label_percent(),
                     expand = expansion(mult = c(0, 0))) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of annotated cells within category", fill = NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        legend.position = "right")

epitope_donor_counts <- adataEpi %>%
  count(sample, condition, antigen.species, name = "annotated_cells")
epitope_coverage <- antigen_meta %>%
  mutate(annotated = !is.na(antigen.species) & trimws(antigen.species) != "") %>%
  group_by(condition) %>%
  summarise(total_cells = n(), annotated_cells = sum(annotated),
            fraction_annotated = mean(annotated), .groups = "drop")

write.csv(epitope_counts, file.path(pathTO, "epitope_category_counts.csv"), row.names = FALSE)
write.csv(epitope_donor_counts, file.path(pathTO, "epitope_donor_counts.csv"), row.names = FALSE)
write.csv(epitope_coverage, file.path(pathTO, "epitope_annotation_coverage.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Save Figure 4 R panels -------------------------------------
# ------------------------------------------------------------

figure4A
figure4B
figure4C
figure4D

ggsave(plot = figure4A, filename = file.path(pathTO, "figure4A.pdf"),
       width = 11, height = 4)
ggsave(plot = figure4A, filename = file.path(pathTO, "figure4A.png"),
       width = 11, height = 4, dpi = 300, bg = "white")
ggsave(plot = figure4B, filename = file.path(pathTO, "figure4B_clone_size.pdf"),
       width = 5, height = 4)
ggsave(plot = figure4B, filename = file.path(pathTO, "figure4B_clone_size.png"),
       width = 5, height = 4, dpi = 300, bg = "white")
figure4C
ggsave(plot = figure4C, filename = file.path(usage_pathTO, "figure4C_gene_usage.pdf"),
       width = 5.5, height = 6.7)
ggsave(plot = figure4C, filename = file.path(usage_pathTO, "figure4C_gene_usage.png"),
       width = 5.5, height = 6.7, dpi = 300, bg = "white")
ggsave(plot = sensitivity_plot, filename = file.path(usage_pathTO, "gene_usage_clonotype_weighted.pdf"),
       width = 5.5, height = 6.7)

ggsave(plot = figure4D, filename = file.path(pathTO, "figure4D.pdf"),
       width = 6, height = 7)
ggsave(plot = figure4D, filename = file.path(pathTO, "figure4D.png"),
       width = 6, height = 7, dpi = 300, bg = "white")
ggsave(plot = donor_expansion_plot,
       filename = file.path(pathTO, "donor_expansion_cluster.pdf"),
       width = 9, height = 4)

writeLines(c(
  paste("Seed:", analysis_seed),
  paste("C exact donor-label assignments:", ncol(assignments)),
  paste("C bootstrap resamples:", n_bootstrap),
  "C: Primary-chain V/J usage with cell and equal-clonotype weighting.",
  paste("C: BH across", sum(gene_usage_tests$weighting == "Cell"),
        "observed genes, separately per weighting."),
  "C: Bootstrap intervals are pointwise and not multiplicity-adjusted.",
  "A-B: Single-pair TCR cells with updated cluster and subset annotations.",
  "D: Antigen annotations in the full cohort; pooled cell counts.",
  "D denominator: annotated cells in each category across both cohorts.",
  "D display: categories with at least one cell in each cohort.",
  capture.output(sessionInfo())
), file.path(pathTO, "session_info.txt"))
