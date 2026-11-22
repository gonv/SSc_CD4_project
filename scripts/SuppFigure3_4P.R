#!/usr/bin/env Rscript

# Supplementary Figure - TCR clonotype sharing and IFN signatures
# CD4 T cell single-cell RNA-seq analysis in SSc
# Run from the 02_CD4 project root, as for Figure2_4P.R and Figure3_4P.R.
# All analyses start from the original cell metadata and expression matrix.
# Required packages: ggplot2, dplyr, tidyr, Matrix and patchwork.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(patchwork)
})

# ------------------------------------------------------------
# Create output directory and define analysis settings
# ------------------------------------------------------------

pathTO <- "figures/plots4figures/Figure4_supplementary"
dir.create(pathTO, recursive = TRUE, showWarnings = FALSE)

# Use the settings applied to the reviewed analyses.
n_perm <- 4999L
min_cells <- 5L
subsets <- c("Tcm", "Temra")

# Match the cohort colors in Figure2_4P.R and Figure3_4P.R.
custom_colors <- c(
  "CTRL" = "#D8BFD8",
  "SSc" = "#C3E5C3"
)

# ------------------------------------------------------------
# Load and match metadata
# ------------------------------------------------------------

# Use the latest metadata in ShinyCD4/data, as requested.
meta <- read.csv("ShinyCD4/data/20260124_SSc_CD4_meta_v3.csv",
                 row.names = 1,
                 stringsAsFactors = FALSE,
                 colClasses = "character",
                 na.strings = c("", "NA", "nan"))

tcr <- read.csv("../data/SSc_CD4_meta_TCR.csv",
                stringsAsFactors = FALSE,
                colClasses = "character",
                na.strings = c("", "NA", "nan"))

meta_fields <- c("cell", "sample", "condition", "cell_type", "cell_type_Th")
tcr_fields <- c("cell_id", "sample", "condition", "clone_id",
                "clone_id_size", "chain_pairing")
if (!all(meta_fields %in% names(meta)) || !all(tcr_fields %in% names(tcr))) {
  stop("Required annotation or clonotype columns are missing.")
}
if (!nrow(meta) || !nrow(tcr) || anyNA(meta[meta_fields]) ||
    anyNA(tcr[tcr_fields]) || anyDuplicated(meta$cell) || anyDuplicated(tcr$cell_id)) {
  stop("Check metadata for missing values or duplicate cell IDs.")
}

# Match by cell ID. 'sample' retains the donor IDs across metadata versions.
idx <- match(tcr$cell_id, meta$cell)

if (anyNA(idx) || !all(tcr$sample == meta$sample[idx] &
                      tcr$condition == meta$condition[idx])) {
  stop("Cell, donor or condition mismatch after joining the metadata.")
}

if (!all(tcr$chain_pairing == "single pair") ||
    !setequal(tcr$condition, c("SSc", "CTRL")) ||
    anyDuplicated(unique(tcr[c("sample", "condition")])$sample)) {
  stop("Expected single-pair TCR cells and one condition per donor.")
}

tcr$major <- meta$cell_type[idx]
tcr$helper <- meta$cell_type_Th[idx]

tcr$helper[tcr$helper %in% c("Th1_GZMK", "Th1_GZMB")] <- "Th1"

# ------------------------------------------------------------
# Clonotype sharing within donors - Panels A and B
# ------------------------------------------------------------

# Null model: shuffle state labels within each donor across all TCR cells.
# This preserves clone sizes and state abundance. A clone contributes once
# to each occupied pair of states, regardless of its cell count.
# Donor fractions use all expanded clonotypes in that donor as denominator.
# Cohort statistics give equal weight to each evaluable donor.
by_donor <- split(tcr, tcr$sample)
set.seed(2024)
overlap_tables <- list()

for (resolution in c("major", "helper")) {
  states <- sort(unique(tcr[[resolution]]))
  pairs <- t(combn(seq_along(states), 2))
  state1 <- pairs[, 1]
  state2 <- pairs[, 2]
  keys <- cbind(state1, state2)
  analyses <- lapply(by_donor, function(donor) {
    clone <- match(donor$clone_id, unique(donor$clone_id))
    sizes <- tabulate(clone)
    if (!isTRUE(all(sizes[clone] == as.numeric(donor$clone_id_size)))) stop("Clone size mismatch.")
    expanded <- sizes[clone] >= 2L
    group <- match(clone[expanded], which(sizes >= 2L))
    n_clone <- sum(sizes >= 2L)
    n_state <- length(states)
    labels <- match(donor[[resolution]], states)
    state_counts <- tabulate(labels, n_state)
    counts <- function(label_subset) {
      presence <- matrix(tabulate(group + (label_subset - 1L) * n_clone,
                                  n_clone * n_state), nrow = n_clone, ncol = n_state) > 0L
      crossprod(presence)[keys]
    }
    observed <- counts(labels[expanded])
    # Draw labels without replacement for expanded cells only. The other
    # positions are irrelevant to the statistic, making this a full shuffle marginal.
    simulated <- replicate(n_perm, counts(labels[sample.int(length(labels), length(group))]))
    support <- n_clone > 0 & state_counts[state1] > 0 & state_counts[state2] > 0
    detail <- data.frame(donor = donor$sample[1], condition = donor$condition[1],
                         state1 = states[state1], state2 = states[state2],
                         n_state1 = state_counts[state1], n_state2 = state_counts[state2],
                         n_expanded_clones = n_clone, supported = support,
                         observed = observed, expected = rowMeans(simulated))
    detail$delta_pp <- ifelse(support, 100 * (detail$observed - detail$expected) / n_clone, NA)
    list(detail = detail, simulated = simulated)
  })

  donor_results <- do.call(rbind, lapply(analyses, "[[", "detail"))
  summary <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    eligible <- Filter(function(x) x$detail$supported[i], analyses)
    if (!length(eligible)) return(NULL)
    data <- do.call(rbind, lapply(eligible, function(x) x$detail[i, ]))
    observed <- mean(data$observed / data$n_expanded_clones)
    simulated <- Reduce("+", lapply(eligible, function(x) {
      x$simulated[i, ] / x$detail$n_expanded_clones[i]
    })) / nrow(data)
    # Equal weight for each evaluable donor. Use a conservative two-sided
    # Monte Carlo p value based on the smaller of the two inclusive tails.
    lower <- (1 + sum(simulated <= observed + 1e-12)) / (1 + n_perm)
    upper <- (1 + sum(simulated >= observed - 1e-12)) / (1 + n_perm)
    data.frame(state1 = data$state1[1], state2 = data$state2[1],
               n_donors = nrow(data), n_donors_with_overlap = sum(data$observed > 0),
               observed_clones = sum(data$observed), expected_clones = sum(data$expected),
               observed_pct = 100 * observed, expected_pct = 100 * mean(simulated),
               delta_pp = mean(data$delta_pp), p_value = min(1, 2 * min(lower, upper)))
  }))
  summary$p_BH <- p.adjust(summary$p_value, method = "BH", n = nrow(pairs))
  write.csv(donor_results, file.path(pathTO, paste0(resolution, "_overlap_donors.csv")), row.names = FALSE)
  write.csv(summary, file.path(pathTO, paste0(resolution, "_overlap_summary.csv")), row.names = FALSE)

  overlap_tables[[resolution]] <- list(donors = donor_results, summary = summary)
}

# ------------------------------------------------------------
# Donor-level SSc versus CTRL comparisons - Panel C
# ------------------------------------------------------------

# Primary endpoint: observed minus expected overlap in percentage points.
# Secondary endpoint: observed proportion of expanded clones shared by each pair.
# Enumerate all donor-label assignments. Adjust within each resolution and
# endpoint. Leave-one-donor-out ranges describe sensitivity, not confidence intervals.
compare <- function(data, endpoint) {
  data <- data[data$supported & is.finite(data[[endpoint]]), ]
  if (anyDuplicated(data$donor)) stop("Duplicate donors.")
  x <- data[[endpoint]]
  is_case <- data$condition == "SSc"
  n_case <- sum(is_case)
  n_ctrl <- sum(!is_case)
  result <- data.frame(n_SSc = n_case, n_CTRL = n_ctrl,
                       mean_SSc = mean(x[is_case]), mean_CTRL = mean(x[!is_case]))
  result$difference_pp <- result$mean_SSc - result$mean_CTRL
  result$p_value <- result$loo_min_pp <- result$loo_max_pp <- NA_real_
  if (min(n_case, n_ctrl) < 2L) return(result)
  if (choose(length(x), n_case) > 1e6) stop("Too many assignments for exact enumeration.")
  null <- combn(seq_along(x), n_case, function(i) {
    mean(x[i]) - mean(x[-i])
  })
  result$p_value <- mean(abs(null) >= abs(result$difference_pp) - 1e-12)
  loo <- vapply(seq_along(x), function(i) {
    mean(x[-i][is_case[-i]]) - mean(x[-i][!is_case[-i]])
  }, numeric(1))
  result$loo_min_pp <- min(loo)
  result$loo_max_pp <- max(loo)
  result
}

all_results <- list()

for (resolution in c("major", "helper")) {
  data <- overlap_tables[[resolution]]$donors
  if (!setequal(data$condition, c("SSc", "CTRL"))) stop("Unexpected cohort labels.")
  data$observed_pp <- 100 * data$observed / data$n_expanded_clones
  strata <- split(data, interaction(data$state1, data$state2, drop = TRUE))
  results <- do.call(rbind, lapply(strata, function(group) {
    do.call(rbind, lapply(c("delta_pp", "observed_pp"), function(endpoint) {
      cbind(resolution, group[1, c("state1", "state2")], endpoint,
            compare(group, endpoint))
    }))
  }))
  # BH within each annotation resolution and endpoint across all tested pairs.
  results$p_BH <- ave(results$p_value, results$endpoint,
                       FUN = function(p) p.adjust(p, "BH"))
  all_results[[resolution]] <- results

}
condition_tests <- do.call(rbind, all_results)
write.csv(condition_tests, file.path(pathTO, "condition_tests.csv"), row.names = FALSE)

# ------------------------------------------------------------
# IFN signatures - Panel D
# ------------------------------------------------------------

IFNa_genes <- c(
  "IFITM1","TRAFD1","CD74","CASP1","PLSCR1","PSMB9","PSME1","PSME2",
  "IFI35","UBE2L6","B2M","LY6E","PSMB8","IFI30","IL7","TAP1","RTP4",
  "IFITM2","PARP9","PARP12","HELZ2","USP18","RNF31","GBP2","OAS1",
  "BST2","LAMP3","MVB12A","UBA7","TRIM21","IFI44","LPAR6",
  "LGALS3BP","STAT2","SP110","NUB1","MX1","NMI","DDX60",
  "TDRD7","CASP8"
)

IFNg_genes <- c(
  "RBCK1","TRAFD1","PTPN6","CD74","PSMB10","CASP1","PLSCR1","PSMB9",
  "PSME1","PSME2","IFNAR2","IFI35","UBE2L6","CASP4","MVP","B2M",
  "LY6E","PSMB8","IFI30","PIM1","IL7","ZBP1","TAP1","RTP4",
  "VAMP5","VAMP8","IFITM2","PARP12","TNFSF10","OAS2","HELZ2",
  "USP18","RNF31","MYD88","BST2","JAK2","HLA-B","ITGB7",
  "STAT1","TRIM21","PSMA2","PML","FAS","IFI44","LGALS3BP",
  "STAT3","IL15RA","STAT2","IL18BP","SP110","GCH1","MX1",
  "HLA-DRB1","HLA-A","NMI","DDX60","SRI","TDRD7",
  "PSMB2","CD274","ARID5B","CASP8"
)


# Determine expansion across the entire retained repertoire of each donor.
# Then use the latest cell annotations to assign each cell to a subset.
cells <- tcr %>%
  transmute(cell = cell_id, donor = sample, condition, clone_id,
            reported_size = as.numeric(clone_id_size), subset = tcr$major) %>%
  add_count(donor, clone_id, name = "clone_size")

cells$status <- ifelse(cells$clone_size >= 2, "Expanded", "Singleton")

# Use the same score as Figure 2F: mean normalized expression of detected
# signature genes. Match columns by cell ID before computing the means.

exprs_norm <- readRDS("../ShinyCD4/data/20260124_SSc_CD4_counts_norm.rds")

expr_index <- match(cells$cell, colnames(exprs_norm))

signatures <- list(IFN_alpha = IFNa_genes, IFN_gamma = IFNg_genes)

coverage <- bind_rows(lapply(names(signatures), function(module) {
  data.frame(module, gene = signatures[[module]],
             present = signatures[[module]] %in% rownames(exprs_norm))
}))

for (module in names(signatures)) {
  genes <- coverage$gene[coverage$module == module & coverage$present]
  if (!length(genes)) stop("No genes found for ", module)
  cells[[module]] <- Matrix::colMeans(exprs_norm[genes, expr_index, drop = FALSE])
  if (any(!is.finite(cells[[module]]))) stop("Non-finite scores for ", module)
}

#rm(exprs_norm)

# Use Tcm and Temra as the primary subsets. Summarize both per cell and
# per clonotype. Clonotypes can contribute to multiple subsets.

scored <- cells %>%
  filter(subset %in% subsets) %>%
  pivot_longer(all_of(names(signatures)), names_to = "module", values_to = "score")

by_status <- scored %>%
  group_by(subset, module, donor, condition, status, clone_id) %>%
  summarise(clone_mean = mean(score), n_cells = n(), .groups = "drop_last") %>%
  summarise(cell_mean = weighted.mean(clone_mean, n_cells),
            clone_mean = mean(clone_mean), n_clones = n(),
            n_cells = sum(n_cells), .groups = "drop")

# Keep donors lacking one class visible in the eligibility table.
pairs <- by_status %>%
  complete(nesting(subset, module, donor, condition),
           status = c("Expanded", "Singleton"),
           fill = list(n_cells = 0L, n_clones = 0L)) %>%
  pivot_wider(names_from = status, values_from = c(cell_mean, clone_mean, n_clones, n_cells)) %>%
  mutate(eligible = n_cells_Expanded >= min_cells & n_cells_Singleton >= min_cells,
         cell_delta = cell_mean_Expanded - cell_mean_Singleton,
         clone_delta = clone_mean_Expanded - clone_mean_Singleton)

# Each donor has equal weight. A one-sample t test of donor differences is
# equivalent to a paired t test on the donor means. Intervals are nominal.
# The clone-weighted check uses exactly the same eligible donor pairs.
estimate <- function(x) {
    n <- length(x)
  effect <- if (n) mean(x) else NA_real_
  test <- if (n >= 3L && sd(x) > 0) tryCatch(t.test(x), error = function(e) NULL)
  tibble(n_donors = n, mean_delta = effect,
         lower_95 = if (is.null(test)) NA_real_ else test$conf.int[1],
         upper_95 = if (is.null(test)) NA_real_ else test$conf.int[2],
         p_value = if (is.null(test)) NA_real_ else test$p.value)
}

results <- pairs %>%
  group_by(subset, module) %>%
  group_modify(~ bind_rows(
    mutate(estimate(.x$cell_delta[.x$eligible]), weighting = "Cell"),
    mutate(estimate(.x$clone_delta[.x$eligible]), weighting = "Clonotype")
  )) %>%
  ungroup() %>%
  group_by(weighting) %>%
  mutate(p_BH = p.adjust(p_value, "BH", n = length(subsets) * length(signatures))) %>%
  ungroup()

write.csv(pairs, file.path(pathTO, "donor_IFN_pairs.csv"), row.names = FALSE)
write.csv(results, file.path(pathTO, "donor_IFN_tests.csv"), row.names = FALSE)
write.csv(coverage, file.path(pathTO, "IFN_signature_coverage.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Plotting - Panels A to D
# ------------------------------------------------------------

base <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        plot.margin = margin(8, 10, 8, 8))

heatmap <- function(resolution, states, title) {
  d <- overlap_tables[[resolution]]$summary
  d$state1[d$state1 == "Thf"] <- "Tfh"
  d$state2[d$state2 == "Thf"] <- "Tfh"
  # Match pair orientation to the displayed state ordering.
  reverse <- match(d$state1, states) > match(d$state2, states)
  old <- d$state1
  d$state1[reverse] <- d$state2[reverse]
  d$state2[reverse] <- old[reverse]
  ggplot(d, aes(state1, state2, fill = delta_pp)) +
    geom_tile(colour = "white", linewidth = .8) +
    geom_text(aes(label = sprintf("%+.1f%s", delta_pp,
                                  ifelse(p_BH < .05, "*", ""))), size = 3) +
    scale_x_discrete(limits = states, drop = FALSE) +
    scale_y_discrete(limits = rev(states), drop = FALSE) +
    scale_fill_gradient2(low = "#397BA5", mid = "#FAFAF8", high = "#D46C43",
                         midpoint = 0, limits = c(-23,23),
                         breaks = c(-20,-10,0,10,20),
                         name = "Observed - expected\n(pp)") +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    base + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                 legend.key.height = grid::unit(.4, "cm"))
}
overlap_major_plot <- heatmap("major", c("Tn","Tcm","Tem","Temra","Treg"), "Clonotype sharing across major subsets")
overlap_helper_plot <- heatmap("helper", c("Tn","Th1","Th2","Th17","Tfh","Treg"), "Clonotype sharing across annotated states")

major <- overlap_tables$major$donors
helper <- overlap_tables$helper$donors
major <- subset(major, state1 == "Tem" & state2 == "Temra")
helper <- subset(helper, state1 == "Th17" & state2 == "Th2")
major$pair <- "Tem / Temra"
helper$pair <- "Th17 / Th2"
d <- rbind(major, helper)
d$condition <- factor(d$condition, c("CTRL","SSc"))
annotations <- condition_tests %>%
  filter(endpoint == "delta_pp",
         (resolution == "major" & state1 == "Tem" & state2 == "Temra") |
         (resolution == "helper" & state1 == "Th17" & state2 == "Th2")) %>%
  mutate(pair = paste(state1, state2, sep = " / "),
         label = sprintf("BH P = %.3f", p_BH))
cohort_plot <- ggplot(d, aes(condition, delta_pp)) +
  geom_hline(yintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_point(aes(fill = condition), shape = 21, colour = "grey25", size = 2,
             position = position_jitter(width = .09, height = 0, seed = 42)) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 8) +
  geom_text(data = annotations, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 3) +
  facet_wrap(~pair, nrow = 1) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.08,.18))) +
  labs(title = "Clonotype sharing by cohort", x = NULL,
       y = "Observed - expected overlap (pp)") +
  base + theme(strip.background = element_blank(),
               strip.text = element_text(face = "bold"))

ifn <- pairs
ifn <- subset(ifn, eligible)
ifn$condition <- factor(ifn$condition, c("CTRL","SSc"))
ifn$module <- factor(ifn$module, c("IFN_alpha","IFN_gamma"),
                     labels = c("IFN-alpha","IFN-gamma"))
ifn_plot <- ggplot(ifn, aes(condition, cell_delta)) +
  geom_hline(yintercept = 0, colour = "grey65", linetype = "dashed") +
  geom_point(aes(fill = condition), shape = 21, colour = "grey25", size = 1.8,
             position = position_jitter(width = .09, height = 0, seed = 42)) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 7) +
  facet_grid(module ~ subset, scales = "free_y") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  labs(title = "IFN scores in expanded and singleton cells", x = NULL,
       y = "IFN score difference\n(expanded - singleton)") +
  base + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

ifn_n <- pairs %>%
  filter(eligible) %>%
  distinct(subset, donor) %>%
  count(subset)
ifn_counts <- paste(paste0(ifn_n$subset, ", n = ", ifn_n$n), collapse = "; ")
primary_ifn_tests <- results %>% filter(weighting == "Cell")
ifn_test_caption <- if (all(!is.na(primary_ifn_tests$p_BH)) &&
                        all(primary_ifn_tests$p_BH > 0.05)) {
  "All four IFN comparisons BH P > 0.05."
} else {
  "See donor_IFN_tests.csv for IFN comparisons."
}

figure <- ((overlap_major_plot | overlap_helper_plot) / (cohort_plot | ifn_plot)) +
  plot_layout(heights = c(1,1), widths = c(1,1)) +
  plot_annotation(tag_levels = "A",
    caption = paste("A-B: Equal-weight donor means; * BH-adjusted permutation P < 0.05.",
                    "C-D: Each point represents one donor; black marks show means.",
                    paste0("D: ", ifn_counts, ". Minimum ", min_cells, " cells per category. ", ifn_test_caption),
                    sep = "\n"),
    theme = theme(plot.tag = element_text(face = "bold", size = 15),
                  plot.caption = element_text(size = 9, hjust = 0, lineheight = 1.15)))



# ------------------------------------------------------------
# Save the supplementary figure and individual panels
# ------------------------------------------------------------

# Keep the figure identifier generic until the supplement is numbered.
figure

ggsave(plot = figure,
       filename = file.path(pathTO, "Figure4_supplementary.png"),
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave(plot = figure,
       filename = file.path(pathTO, "Figure4_supplementary.pdf"),
       width = 12, height = 10, bg = "white")

figure_panels <- list(A = overlap_major_plot, B = overlap_helper_plot,
                      C = cohort_plot, D = ifn_plot)
for (panel_name in names(figure_panels)) {
  ggsave(plot = figure_panels[[panel_name]],
         filename = file.path(pathTO, paste0("Figure4_supplementary_", panel_name, ".pdf")),
         width = 6, height = 5, bg = "white")
}

writeLines(c(paste("Within-donor permutations:", n_perm),
             "Permutation seed: 2024",
             paste("Minimum IFN cells per donor, subset and class:", min_cells),
             capture.output(sessionInfo())),
           file.path(pathTO, "session_info.txt"))
