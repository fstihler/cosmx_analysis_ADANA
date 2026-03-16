# PanCK-gated two-pass InSituType cell typing pipeline
# for CosMx 6K breast cancer TMA spatial transcriptomics data
#
# Strategy: Split cells by PanCK protein fluorescence into epithelial vs.
# non-epithelial BEFORE running InSituType. Type each subset separately.
#
# Author: Frederik Stihler

library(Seurat)
library(Matrix)
library(InSituType)
library(ggplot2)
library(data.table)
library(dplyr)


# =============================================================================
# Function 1: determine_panck_threshold
# =============================================================================

determine_panck_threshold <- function(seurat_obj,
                                       panck_col = "Mean.PanCK",
                                       threshold = "auto",
                                       output_dir = "results/celltyping") {

  if (!panck_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", panck_col, "' not found in metadata.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panck_vals <- seurat_obj@meta.data[[panck_col]]

  # Auto-detect threshold from bimodal density minimum
  if (identical(threshold, "auto")) {
    log_vals <- log1p(panck_vals)
    dens <- density(log_vals, n = 1024)
    # Find local minima: where derivative goes from negative to positive
    dy <- diff(dens$y)
    minima_idx <- which(dy[-length(dy)] < 0 & dy[-1] > 0) + 1
    if (length(minima_idx) == 0) {
      # Fallback: use median if no valley found
      threshold <- median(panck_vals, na.rm = TRUE)
      message("No bimodal valley detected, using median as threshold: ", round(threshold))
    } else {
      # Pick the deepest valley (lowest density at the minimum)
      best_idx <- minima_idx[which.min(dens$y[minima_idx])]
      threshold <- expm1(dens$x[best_idx])
      message("Auto-detected bimodal valley at log1p = ", round(dens$x[best_idx], 3),
              " -> threshold = ", round(threshold))
    }
  }

  gate <- ifelse(panck_vals >= threshold, "PanCK_high", "PanCK_low")

  n_high <- sum(gate == "PanCK_high")
  n_low  <- sum(gate == "PanCK_low")
  n_total <- length(gate)

  message("PanCK threshold: ", threshold)
  message("  PanCK_high (epithelial):     ", n_high, " (", round(100 * n_high / n_total, 1), "%)")
  message("  PanCK_low  (non-epithelial): ", n_low,  " (", round(100 * n_low / n_total, 1), "%)")

  # Density plot
  df <- data.frame(log1p_panck = log1p(panck_vals))

  p_density <- ggplot(df, aes(x = log1p_panck)) +
    geom_density(fill = "steelblue", alpha = 0.4) +
    geom_vline(xintercept = log1p(threshold), color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = log1p(threshold), y = Inf, label = paste0("threshold = ", threshold),
             hjust = -0.1, vjust = 2, color = "red", size = 3.5) +
    labs(title = "PanCK Gate: Density of log1p(Mean.PanCK)",
         x = "log1p(Mean.PanCK)", y = "Density") +
    theme_bw()

  ggsave(file.path(output_dir, "panck_gate_density.pdf"), p_density, width = 8, height = 5)
  print(p_density)

  # Violin plot
  df_vln <- data.frame(
    panck = panck_vals,
    gate = factor(gate, levels = c("PanCK_low", "PanCK_high"))
  )

  p_violin <- ggplot(df_vln, aes(x = gate, y = log1p(panck), fill = gate)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_fill_manual(values = c("PanCK_low" = "grey70", "PanCK_high" = "tomato")) +
    labs(title = "PanCK Gate: Violin of log1p(Mean.PanCK) by group",
         x = "Gate", y = "log1p(Mean.PanCK)") +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "panck_gate_violin.pdf"), p_violin, width = 6, height = 5)
  print(p_violin)

  list(
    threshold = threshold,
    gate = gate,
    n_high = n_high,
    n_low = n_low,
    n_total = n_total
  )
}


# =============================================================================
# Function 2: split_by_panck
# =============================================================================

split_by_panck <- function(seurat_obj,
                           panck_col = "Mean.PanCK",
                           threshold = 3500) {

  if (!panck_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", panck_col, "' not found in metadata.", call. = FALSE)
  }

  panck_vals <- seurat_obj@meta.data[[panck_col]]
  seurat_obj$PanCK_gate <- ifelse(panck_vals >= threshold, "PanCK_high", "PanCK_low")

  cells_high <- colnames(seurat_obj)[seurat_obj$PanCK_gate == "PanCK_high"]
  cells_low  <- colnames(seurat_obj)[seurat_obj$PanCK_gate == "PanCK_low"]

  message("Splitting by PanCK gate (threshold = ", threshold, "):")
  message("  PanCK_high (epithelial):     ", length(cells_high), " cells")
  message("  PanCK_low  (non-epithelial): ", length(cells_low), " cells")

  seu_epithelial    <- subset(seurat_obj, cells = cells_high)
  seu_nonepithelial <- subset(seurat_obj, cells = cells_low)

  list(
    nonepithelial = seu_nonepithelial,
    epithelial = seu_epithelial,
    full = seurat_obj
  )
}


# =============================================================================
# Function 3: type_nonepithelial
# =============================================================================

type_nonepithelial <- function(seu_nonepithelial,
                                reference_profiles,
                                neg_probes = NULL,
                                n_free_clusters = 3,
                                min_posterior = 0.5,
                                cohorting = TRUE,
                                n_cohorts = 10,
                                output_dir = "results/celltyping") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message(Sys.time(), " - Typing non-epithelial cells (PanCK-low)...")

  # Extract counts (cells x genes)
  counts <- Matrix::t(seu_nonepithelial[["RNA"]]$counts)

  # Extract or use provided negative probes
  if (is.null(neg_probes)) {
    negmean <- Matrix::colMeans(seu_nonepithelial[["negprobes"]])
  } else {
    negmean <- neg_probes
  }

  # Remove Malignancy from reference
  ref <- reference_profiles
  if ("Malignancy" %in% rownames(ref)) {
    message("  Removing 'Malignancy' from reference profiles")
    ref <- ref[rownames(ref) != "Malignancy", , drop = FALSE]
  }

  # Filter reference to genes present in counts
  ref <- ref[is.element(rownames(ref), colnames(counts)), , drop = FALSE]
  message("  Reference: ", ncol(ref), " cell types, ", nrow(ref), " genes")

  # Cohorting
  coh <- NULL
  if (cohorting) {
    md <- seu_nonepithelial@meta.data
    cohort_cols <- intersect(c("Mean.PanCK", "Mean.CD45"), colnames(md))
    if (length(cohort_cols) > 0) {
      cohort_features <- as.matrix(md[, cohort_cols, drop = FALSE])
      coh <- fastCohorting(cohort_features, gaussian_transform = TRUE, n_cohorts = n_cohorts)
      message("  Cohorting with: ", paste(cohort_cols, collapse = ", "), " (", n_cohorts, " cohorts)")
    } else {
      message("  No cohort columns found, skipping cohorting")
    }
  }

  # Run semi-supervised InSituType
  message("  Running semi-supervised InSituType (", ncol(ref), " reference types + ", n_free_clusters, " free clusters)...")
  result <- insitutype(
    x = counts,
    neg = negmean,
    cohort = coh,
    bg = NULL,
    n_clusts = n_free_clusters,
    reference_profiles = ref,
    update_reference_profiles = FALSE
  )

  # Extract clusters
  clusters <- result$clust
  message("  Found clusters: ", paste(sort(unique(clusters)), collapse = ", "))

  # Compute posteriors from logliks
  logliks <- result$logliks
  posteriors <- exp(logliks - apply(logliks, 1, max))
  posteriors <- posteriors / rowSums(posteriors)
  max_post <- apply(posteriors, 1, max)

  # Flag low-confidence cells
  n_flagged <- sum(max_post < min_posterior)
  if (n_flagged > 0) {
    message("  Flagging ", n_flagged, " cells with max posterior < ", min_posterior, " as Unassigned_NonEpi")
    clusters[max_post < min_posterior] <- "Unassigned_NonEpi"
  }

  # Add to metadata
  cell_ids <- colnames(seu_nonepithelial)
  cluster_metadata <- rep(NA_character_, length(cell_ids))
  names(cluster_metadata) <- cell_ids
  cluster_metadata[names(clusters)] <- clusters

  seu_nonepithelial$celltype_nonepithelial <- as.factor(cluster_metadata)
  seu_nonepithelial$max_posterior_nonepithelial <- NA_real_
  seu_nonepithelial$max_posterior_nonepithelial[match(names(max_post), cell_ids)] <- max_post

  # Save cluster assignments
  assign_df <- data.frame(
    cell_id = names(clusters),
    celltype = clusters,
    max_posterior = max_post[names(clusters)],
    stringsAsFactors = FALSE
  )
  data.table::fwrite(assign_df, file.path(output_dir, "nonepithelial_cluster_assignments.csv"))

  # Flightpath PDF
  pdf(file.path(output_dir, "nonepithelial_flightpath.pdf"), width = 10, height = 8)
  tryCatch({
    cols <- InSituType::colorCellTypes(freqs = table(result$clust), palette = "brewers")
    flightpath <- InSituType::flightpath_layout(logliks = result$logliks, profiles = result$profiles)
    par(mar = c(0, 0, 2, 0))
    plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[result$clust],
         main = "Non-epithelial Flightpath")
    text(flightpath$clustpos[, 1], flightpath$clustpos[, 2],
         rownames(flightpath$clustpos), cex = 0.7)
    par(mar = c(5, 4, 4, 2) + 0.1)
  }, error = function(e) {
    message("  Warning: flightpath plot failed: ", e$message)
  })
  dev.off()

  message(Sys.time(), " - Non-epithelial typing complete: ", length(unique(clusters)), " clusters, ",
          length(clusters), " cells")

  list(
    seurat_obj = seu_nonepithelial,
    ist_result = result
  )
}


# =============================================================================
# Function 4: type_epithelial
# =============================================================================

type_epithelial <- function(seu_epithelial,
                             neg_probes = NULL,
                             n_clusters = 10,
                             min_posterior = 0.5,
                             cohorting = TRUE,
                             n_cohorts = 10,
                             output_dir = "results/celltyping") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message(Sys.time(), " - Typing epithelial cells (PanCK-high)...")

  # Extract counts (cells x genes)
  counts <- Matrix::t(seu_epithelial[["RNA"]]$counts)

  # Extract or use provided negative probes
  if (is.null(neg_probes)) {
    negmean <- Matrix::colMeans(seu_epithelial[["negprobes"]])
  } else {
    negmean <- neg_probes
  }

  # Cohorting (PanCK only for epithelial)
  coh <- NULL
  if (cohorting) {
    md <- seu_epithelial@meta.data
    if ("Mean.PanCK" %in% colnames(md)) {
      cohort_features <- as.matrix(md[, "Mean.PanCK", drop = FALSE])
      coh <- fastCohorting(cohort_features, gaussian_transform = TRUE, n_cohorts = n_cohorts)
      message("  Cohorting with: Mean.PanCK only (", n_cohorts, " cohorts)")
    } else {
      message("  Mean.PanCK not found, skipping cohorting")
    }
  }

  # Run unsupervised InSituType (no reference profiles)
  message("  Running unsupervised InSituType (", n_clusters, " clusters)...")
  result <- insitutype(
    x = counts,
    neg = negmean,
    cohort = coh,
    bg = NULL,
    n_clusts = n_clusters,
    reference_profiles = NULL,
    update_reference_profiles = FALSE
  )

  # Extract clusters and rename to Epi_a, Epi_b, ...
  clusters <- result$clust
  orig_labels <- sort(unique(clusters))
  epi_labels <- paste0("Epi_", letters[seq_along(orig_labels)])
  label_map <- setNames(epi_labels, orig_labels)
  clusters <- label_map[clusters]
  names(clusters) <- names(result$clust)

  message("  Cluster mapping: ", paste(paste0(orig_labels, " -> ", epi_labels), collapse = ", "))

  # Compute posteriors from logliks
  logliks <- result$logliks
  posteriors <- exp(logliks - apply(logliks, 1, max))
  posteriors <- posteriors / rowSums(posteriors)
  max_post <- apply(posteriors, 1, max)

  # Flag low-confidence cells
  n_flagged <- sum(max_post < min_posterior)
  if (n_flagged > 0) {
    message("  Flagging ", n_flagged, " cells with max posterior < ", min_posterior, " as Unassigned_Epi")
    clusters[max_post < min_posterior] <- "Unassigned_Epi"
  }

  # Add to metadata
  cell_ids <- colnames(seu_epithelial)
  cluster_metadata <- rep(NA_character_, length(cell_ids))
  names(cluster_metadata) <- cell_ids
  cluster_metadata[names(clusters)] <- clusters

  seu_epithelial$celltype_epithelial <- as.factor(cluster_metadata)
  seu_epithelial$max_posterior_epithelial <- NA_real_
  seu_epithelial$max_posterior_epithelial[match(names(max_post), cell_ids)] <- max_post

  # Save cluster assignments
  assign_df <- data.frame(
    cell_id = names(clusters),
    celltype = clusters,
    max_posterior = max_post[names(result$clust)],
    stringsAsFactors = FALSE
  )
  data.table::fwrite(assign_df, file.path(output_dir, "epithelial_cluster_assignments.csv"))

  # Flightpath PDF
  pdf(file.path(output_dir, "epithelial_flightpath.pdf"), width = 10, height = 8)
  tryCatch({
    cols <- InSituType::colorCellTypes(freqs = table(result$clust), palette = "brewers")
    flightpath <- InSituType::flightpath_layout(logliks = result$logliks, profiles = result$profiles)
    par(mar = c(0, 0, 2, 0))
    plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[result$clust],
         main = "Epithelial Flightpath")
    text(flightpath$clustpos[, 1], flightpath$clustpos[, 2],
         rownames(flightpath$clustpos), cex = 0.7)
    par(mar = c(5, 4, 4, 2) + 0.1)
  }, error = function(e) {
    message("  Warning: flightpath plot failed: ", e$message)
  })
  dev.off()

  message(Sys.time(), " - Epithelial typing complete: ", length(unique(clusters)), " clusters, ",
          length(clusters), " cells")

  list(
    seurat_obj = seu_epithelial,
    ist_result = result
  )
}


# =============================================================================
# Function 5: merge_annotations
# =============================================================================

merge_annotations <- function(seurat_obj_full,
                               seu_nonepithelial,
                               seu_epithelial,
                               output_dir = "results/celltyping") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message(Sys.time(), " - Merging annotations...")

  md_full <- seurat_obj_full@meta.data
  all_cells <- rownames(md_full)

  # Initialize columns
  md_full$celltype_final <- NA_character_
  md_full$max_posterior <- NA_real_

  # Non-epithelial annotations
  nonepi_cells <- colnames(seu_nonepithelial)
  md_nonepi <- seu_nonepithelial@meta.data
  matched_nonepi <- intersect(all_cells, nonepi_cells)
  md_full[matched_nonepi, "celltype_final"] <- as.character(md_nonepi[matched_nonepi, "celltype_nonepithelial"])
  if ("max_posterior_nonepithelial" %in% colnames(md_nonepi)) {
    md_full[matched_nonepi, "max_posterior"] <- md_nonepi[matched_nonepi, "max_posterior_nonepithelial"]
  }

  # Epithelial annotations
  epi_cells <- colnames(seu_epithelial)
  md_epi <- seu_epithelial@meta.data
  matched_epi <- intersect(all_cells, epi_cells)
  md_full[matched_epi, "celltype_final"] <- as.character(md_epi[matched_epi, "celltype_epithelial"])
  if ("max_posterior_epithelial" %in% colnames(md_epi)) {
    md_full[matched_epi, "max_posterior"] <- md_epi[matched_epi, "max_posterior_epithelial"]
  }

  md_full$celltype_final <- as.factor(md_full$celltype_final)

  # Write back to Seurat object
  seurat_obj_full@meta.data <- md_full

  message("  Merged annotations for ", sum(!is.na(md_full$celltype_final)), " / ", nrow(md_full), " cells")
  message("  Final cell types: ", paste(sort(unique(na.omit(md_full$celltype_final))), collapse = ", "))

  # Compute cluster summary
  clusters_vec <- md_full$celltype_final
  total_cells <- nrow(md_full)

  # Find fluorescence columns
  panck_col <- grep("panck", colnames(md_full), ignore.case = TRUE, value = TRUE)
  cd45_col  <- grep("cd45",  colnames(md_full), ignore.case = TRUE, value = TRUE)
  panck_col <- if (length(panck_col) > 0) panck_col[1] else NULL
  cd45_col  <- if (length(cd45_col) > 0)  cd45_col[1]  else NULL

  cluster_ids <- unique(clusters_vec[!is.na(clusters_vec)])

  cluster_summary <- do.call(rbind, lapply(cluster_ids, function(cl) {
    mask <- clusters_vec == cl & !is.na(clusters_vec)
    md_cl <- md_full[mask, , drop = FALSE]
    n <- nrow(md_cl)

    data.frame(
      cluster = cl,
      n_cells = n,
      pct_of_total = round(100 * n / total_cells, 2),
      mean_nCount_RNA = if ("nCount_RNA" %in% colnames(md_cl)) round(mean(md_cl$nCount_RNA, na.rm = TRUE), 2) else NA,
      median_nCount_RNA = if ("nCount_RNA" %in% colnames(md_cl)) round(median(md_cl$nCount_RNA, na.rm = TRUE), 2) else NA,
      mean_nFeature_RNA = if ("nFeature_RNA" %in% colnames(md_cl)) round(mean(md_cl$nFeature_RNA, na.rm = TRUE), 2) else NA,
      mean_PanCK = if (!is.null(panck_col)) round(mean(md_cl[[panck_col]], na.rm = TRUE), 2) else NA,
      mean_CD45 = if (!is.null(cd45_col)) round(mean(md_cl[[cd45_col]], na.rm = TRUE), 2) else NA,
      stringsAsFactors = FALSE
    )
  }))
  rownames(cluster_summary) <- NULL
  cluster_summary <- cluster_summary[order(cluster_summary$n_cells, decreasing = TRUE), ]

  data.table::fwrite(cluster_summary, file.path(output_dir, "all_clusters_summary.csv"))
  message("  Saved cluster summary: ", file.path(output_dir, "all_clusters_summary.csv"))

  # Print summary table
  message("\n  Cluster summary:")
  for (i in seq_len(nrow(cluster_summary))) {
    row <- cluster_summary[i, ]
    message(sprintf("    %-30s  n=%6d  (%.1f%%)", row$cluster, row$n_cells, row$pct_of_total))
  }

  seurat_obj_full
}


# =============================================================================
# Function 6: run_celltyping_pipeline (wrapper)
# =============================================================================

run_celltyping_pipeline <- function(seurat_obj,
                                     reference_profiles,
                                     neg_probes = NULL,
                                     panck_threshold = 3500,
                                     panck_col = "Mean.PanCK",
                                     n_free_clusters_nonepithelial = 3,
                                     n_clusters_epithelial = 10,
                                     min_posterior = 0.5,
                                     cohorting = TRUE,
                                     n_cohorts = 10,
                                     output_dir = "results/celltyping") {

  message(Sys.time(), " - Starting PanCK-gated cell typing pipeline")
  message("=" , strrep("=", 60))

  # Step 1: PanCK threshold
  message(Sys.time(), " - Step 1: PanCK threshold determination")
  gate_result <- determine_panck_threshold(
    seurat_obj, panck_col = panck_col,
    threshold = panck_threshold, output_dir = output_dir
  )

  # Step 2: Split by PanCK
  message(Sys.time(), " - Step 2: Splitting cells by PanCK gate")
  split_result <- split_by_panck(
    seurat_obj, panck_col = panck_col,
    threshold = panck_threshold
  )
  seurat_obj <- split_result$full

  # Step 3: Type non-epithelial cells
  message(Sys.time(), " - Step 3: Typing non-epithelial cells")
  nonepi_result <- type_nonepithelial(
    split_result$nonepithelial,
    reference_profiles = reference_profiles,
    neg_probes = neg_probes,
    n_free_clusters = n_free_clusters_nonepithelial,
    min_posterior = min_posterior,
    cohorting = cohorting, n_cohorts = n_cohorts,
    output_dir = output_dir
  )

  # Step 4: Type epithelial cells
  message(Sys.time(), " - Step 4: Typing epithelial cells")
  epi_result <- type_epithelial(
    split_result$epithelial,
    neg_probes = neg_probes,
    n_clusters = n_clusters_epithelial,
    min_posterior = min_posterior,
    cohorting = cohorting, n_cohorts = n_cohorts,
    output_dir = output_dir
  )

  # Step 5: Merge annotations
  message(Sys.time(), " - Step 5: Merging annotations")
  seurat_obj <- merge_annotations(
    seurat_obj,
    nonepi_result$seurat_obj,
    epi_result$seurat_obj,
    output_dir = output_dir
  )

  message("=" , strrep("=", 60))
  message(Sys.time(), " - Pipeline complete!")
  message("  Final cell types: ", length(unique(na.omit(seurat_obj$celltype_final))))

  seurat_obj
}


# Usage:
# source("celltyping_pipeline.R")
#
# seu <- readRDS("path/to/seurat_object.rds")
# ref <- read.csv("path/to/BreastCancer_6k.profiles.csv", row.names = 1)
#
# result <- run_celltyping_pipeline(
#   seurat_obj = seu,
#   reference_profiles = ref,
#   panck_threshold = 3500,
#   n_free_clusters_nonepithelial = 3,
#   n_clusters_epithelial = 10
# )
#
# # Inspect:
# table(result$celltype_final)
