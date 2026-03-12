# These are helper functions for CosMX data analysis
# import them to your analysis using source("helper_functions.R")
# Author: Frederik Stihler

## IO UTILS

# Function: Load all Seurat objects from a folder into a named list
# Args:
#   folder: directory containing .RDS Seurat objects
# Returns:
#   Named list of Seurat objects
load_seurat_objects <- function(folder, assays_to_remove = c("QC_Normalization.RNA.1_1")) {
  files <- list.files(folder, pattern = "^seuratObject_.*\\.RDS$", full.names = TRUE)
  print(files)
  # Extract slide names from filenames
  slide_names <- sub("^seuratObject_(.*)\\.RDS$", "\\1", basename(files))
  print(slide_names)
  
  # Load into list
  seurat_list <- setNames(lapply(files, readRDS), slide_names)
  
  # Remove specified assays
  seurat_list <- lapply(seurat_list, function(obj) {
    for (assay_name in assays_to_remove) {
      if (assay_name %in% names(obj@assays)) {
        print(paste("Removing assay:", assay_name))
        obj[[assay_name]] <- NULL
      }
    }
    return(obj)
  })

  return(seurat_list)
}

# Function: Merge all Seurat objects with slide-specific cell ID prefixes
# Args:
#   seurat_list: named list of Seurat objects
#   project_name: name for merged object
# Returns:
#   Merged Seurat object
merge_seurat_objects <- function(seurat_list, project_name = "CosMx") {
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    project = project_name
  )
  return(merged)
}

## ANNOTATIONS

make_casewhen <- function(borders) {
  borders <- sort(borders)
  starts <- c(1, borders[-length(borders)] + 1)
  ends <- borders
  regions <- seq_along(borders)
  
  lines <- paste0(
    "fov %in% ", starts, ":", ends, " ~ \"", regions, "\""
  )
  
  case_block <- paste(lines, collapse = ",\n    ")
  full_code <- paste0(
    "case_when(\n    ",
    case_block,
    ",\n    TRUE ~ \"Unknown\"\n  )"
  )
  
  cat(full_code)
  invisible(full_code)
}

process_tma_csv <- function(file_path, sep = ",", skip_rows = 1) {
  # Load CSV
  metadata_csv <- read.csv(file_path, stringsAsFactors = FALSE, sep = sep, skip = skip_rows)
  
  # Rename columns
  colnames(metadata_csv) <- c(
    "study_id",
    "accession_id",
    "accession_year",
    "part_id",
    "block_id",
    "tma_number",
    "tma_map_tumor",
    "tma_map_normal",
    "batch"
  )
  
  # Keep only relevant columns
  metadata_csv <- metadata_csv[, c("study_id", "tma_number", "tma_map_tumor", "tma_map_normal", "batch")]

  # Remove leading and trailing whitespace
  metadata_csv <- metadata_csv %>%
    mutate(across(everything(), ~trimws(.)))
  
  # Convert to long format
  metadata_long <- metadata_csv %>%
    mutate(
      tma_map_tumor = ifelse(tma_map_tumor == "No Tumor", NA, tma_map_tumor),
      tma_map_normal = ifelse(tma_map_normal == "No Normal", NA, tma_map_normal)
    ) %>%
    pivot_longer(
      cols = c(tma_map_tumor, tma_map_normal),
      names_to = "type",
      values_to = "location"
    ) %>%
    filter(!is.na(location)) %>%
    mutate(type = ifelse(type == "tma_map_tumor", "T", "N")) %>%
    separate_rows(location, sep = ",\\s*") %>%
    select(study_id, tma_number, location, batch, type)
  
  return(metadata_long)
}

add_metadata_details <- function(seu.obj, metadata_details, cols_to_keep = NULL) {
  
  if (is.null(cols_to_keep)) {
    cols_to_keep <- c(
    "study_id", "age_at_dx", "race_ethnicity", "ht_m", "wt_kg", "bmi",
    "menopausal_status", "diabetes_status", "alcohol_status_en", "smoking_status_en",
    "stage", "grade", "er_status_registry", "pr_status_registry", "her2_status_registry",
    "er_status_path", "pr_status_path", "her2_status_path",
    "received_chemo", "received_rad", "received_horm", "received_surgery",
    "tumor_size_mm", "reg_nodes_metastasis", "reg_nodes_pos", "reg_nodes_exam", "tma"
  )
  }
  
  metadata_details <- metadata_details %>%
    select(any_of(cols_to_keep))
  
  seu.obj@meta.data <- seu.obj@meta.data %>%
    mutate(temp_cellid = rownames(seu.obj@meta.data)) %>%
    left_join(
      metadata_details,
      by = "study_id"
    ) %>%
    column_to_rownames(var = "temp_cellid")

  
  if (length(unique(na.omit(seu.obj@meta.data$tma))) != 1) {
      print(unique(na.omit(seu.obj@meta.data$tma)))
      stop("tma mismatch detected")
    } else {
      message("tma check passed: only one tma detected.")
      seu.obj@meta.data <- seu.obj@meta.data %>%
        select(-tma)
    }

  return(seu.obj)
}


## VISUALIZATION

xyplot <- function(cluster_column,
                   x_column = "x_slide_mm",
                   y_column = "y_slide_mm",
                   cls = NULL,
                   clusters = NULL,
                   metadata,
                   ptsize = 0.25,
                   plotfirst = NULL,
                   plotfirst_on_top = FALSE,
                   alphasize = 1,
                   show_legend = TRUE,
                   coord_equal = TRUE,
                   continuous_palette = function(n) viridis::viridis(n, option = "plasma"),
                   aes_mappings = list(size = NULL, shape = NULL, alpha = NULL),
                   order = NULL,
                   na_color = "black",
                   theme = ggplot2::theme_bw(),
                   show_labels = FALSE,         
                   label_size = 4,               
                   label_color = "black",        
                   label_fontface = "bold",     
                   label_method = "median",      
                   label_repel = TRUE,           
                   label_max_overlaps = 50       
) {

  pd <- data.table::copy(data.table::data.table(metadata))

  if (!is.null(order)) {
    order <- as.character(order)
  } else if (is.factor(pd[[cluster_column]])) {
    order <- levels(pd[[cluster_column]])
  }

  if (is.null(clusters)) clusters <- unique(pd[[cluster_column]])
  clusters <- unique(clusters)

  mask <- pd[[cluster_column]] %in% clusters
  mask[is.na(mask)] <- FALSE
  if (any(is.na(pd[[cluster_column]]))) mask <- mask | is.na(pd[[cluster_column]])
  pd_plot <- pd[mask, , drop = FALSE]

  if (!is.null(order)) {
    present_levels <- order
    pd_plot[[cluster_column]] <- factor(as.character(pd_plot[[cluster_column]]), levels = present_levels)
    clusters_plot <- present_levels
  } else {
    if (all(!is.na(suppressWarnings(as.numeric(as.character(clusters)))))) {
      clusters_plot <- sort(as.numeric(as.character(clusters)))
      clusters_plot <- as.character(clusters_plot)
    } else {
      clusters_plot <- sort(as.character(clusters))
    }
    if (any(is.na(pd_plot[[cluster_column]]))) clusters_plot <- unique(c(clusters_plot, NA))
  }

  mapping_args <- list(x = x_column, y = y_column, color = cluster_column)
  if (!is.null(aes_mappings$size)) mapping_args$size <- aes_mappings$size
  if (!is.null(aes_mappings$shape)) mapping_args$shape <- aes_mappings$shape
  if (!is.null(aes_mappings$alpha)) mapping_args$alpha <- aes_mappings$alpha
  mapping <- do.call(ggplot2::aes_string, mapping_args)

  size_is_mapped <- !is.null(aes_mappings$size)

  if (!is.null(plotfirst)) {
    plotfirst <- intersect(clusters_plot, plotfirst)
    notplotfirst <- setdiff(clusters_plot, plotfirst)

    df_first <- pd_plot[as.character(pd_plot[[cluster_column]]) %in% as.character(plotfirst), , drop = FALSE]
    df_rest  <- pd_plot[as.character(pd_plot[[cluster_column]]) %in% as.character(notplotfirst), , drop = FALSE]

    layer_first <- ggplot2::geom_point(
      data = df_first, mapping = mapping,
      size = if (!size_is_mapped) ptsize else NULL,
      alpha = if (is.null(aes_mappings$alpha)) 1 else NULL
    )
    layer_rest <- ggplot2::geom_point(
      data = df_rest, mapping = mapping,
      size = if (!size_is_mapped) ptsize else NULL,
      alpha = if (is.null(aes_mappings$alpha)) alphasize else NULL
    )

    p <- ggplot2::ggplot() + theme
    if (plotfirst_on_top) {
      p <- p + layer_rest + layer_first
    } else {
      p <- p + layer_first + layer_rest
    }
  } else {
    p <- ggplot2::ggplot(pd_plot, mapping = mapping) +
      ggplot2::geom_point(
        size = if (!size_is_mapped) ptsize else NULL,
        alpha = if (is.null(aes_mappings$alpha)) alphasize else NULL
      ) +
      theme
  }

  if (is.numeric(pd[[cluster_column]])) {
    cols <- try(continuous_palette(256), silent = TRUE)
    if (inherits(cols, "try-error") || !is.character(cols)) {
      p <- p + ggplot2::scale_color_viridis_c(
        guide = if (show_legend) ggplot2::guide_colorbar() else "none",
        na.value = na_color
      )
    } else {
      p <- p + ggplot2::scale_color_gradientn(
        colors = cols,
        guide = if (show_legend) ggplot2::guide_colorbar() else "none",
        na.value = na_color
      )
    }
  } else {
    n_clusters <- length(clusters_plot)
    if (is.null(cls)) {
      base_pal <- unname(pals::alphabet())
      if (n_clusters <= length(base_pal)) {
        palette_used <- base_pal[seq_len(n_clusters)]
      } else {
        palette_used <- grDevices::colorRampPalette(base_pal)(n_clusters)
      }
      names(palette_used) <- as.character(clusters_plot)
      cls_final <- palette_used
    } else {
      if (is.character(cls) && length(cls) < n_clusters) {
        cls_final <- grDevices::colorRampPalette(cls)(n_clusters)
        names(cls_final) <- as.character(clusters_plot)
      } else {
        if (is.null(names(cls))) {
          cls_final <- cls[seq_len(min(length(cls), n_clusters))]
          names(cls_final) <- as.character(clusters_plot)[seq_len(length(cls_final))]
          if (length(cls) < n_clusters) {
            extra <- grDevices::colorRampPalette(cls)(n_clusters - length(cls))
            cls_final <- c(
              cls_final,
              stats::setNames(extra, as.character(clusters_plot)[(length(cls) + 1):n_clusters])
            )
          }
        } else {
          cls_final <- cls[as.character(clusters_plot)]
          missing_idx <- which(is.na(cls_final))
          if (length(missing_idx) > 0) {
            fill_colors <- grDevices::colorRampPalette(unname(cls))(length(missing_idx))
            cls_final[missing_idx] <- fill_colors
          }
        }
      }
    }

    p <- p + ggplot2::scale_color_manual(
      values = cls_final,
      guide = if (show_legend)
        ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))
      else "none",
      na.value = na_color
    )
  }

  if (isTRUE(coord_equal)) p <- p + ggplot2::coord_fixed()

  # New label layer
  if (isTRUE(show_labels) && !is.numeric(pd_plot[[cluster_column]])) {
    summary_fun <- switch(label_method, "mean" = mean, "median" = median, median)
    centroids <- pd_plot[, .(
      x = summary_fun(get(x_column), na.rm = TRUE),
      y = summary_fun(get(y_column), na.rm = TRUE)
    ), by = cluster_column]

    if (requireNamespace("ggrepel", quietly = TRUE) && isTRUE(label_repel)) {
      p <- p + ggrepel::geom_label_repel(
        data = centroids,
        aes(x = x, y = y, label = .data[[cluster_column]]),
        inherit.aes = FALSE,
        color = label_color,
        size = label_size,
        fontface = label_fontface,
        show.legend = FALSE,
        max.overlaps = label_max_overlaps
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = centroids,
        aes(x = x, y = y, label = .data[[cluster_column]]),
        inherit.aes = FALSE,
        color = label_color,
        size = label_size,
        fontface = label_fontface,
        show.legend = FALSE
      )
    }
  }

  p
}

plot_embedding <- function(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label_size = 3.5,
  alpha = 0.6,
  point_size = 0.3,
  raster = TRUE,
  shuffle = TRUE,
  seed = 123,
  palette = NULL,
  na_color = "grey80",
  palette_continuous = function(n) viridis::viridis(n, option = "plasma"),
  legend = TRUE
) {
  # Validate inputs
  if (!reduction %in% names(seu@reductions))
    stop("Reduction not found: ", reduction)
  if (!group.by %in% colnames(seu@meta.data))
    stop("Column not found in meta.data: ", group.by)

  # Extract embeddings and metadata
  emb <- as.data.frame(Embeddings(seu, reduction))
  colnames(emb)[1:2] <- c("Dim1", "Dim2")
  emb$group <- seu@meta.data[[group.by]]

  # Optional shuffle for even plotting
  if (shuffle) {
    set.seed(seed)
    emb <- emb[sample(nrow(emb)), , drop = FALSE]
  }

  # Define plotting geometry
  geom_fun <- if (raster) ggrastr::geom_point_rast else geom_point

  # Handle numeric vs categorical
  if (is.numeric(emb$group)) {
    # Continuous variable (e.g., expression, score)
    p <- ggplot(emb, aes(Dim1, Dim2, color = group)) +
      geom_fun(size = point_size, alpha = alpha) +
      scale_color_gradientn(
        colors = palette_continuous(256),
        na.value = na_color
      ) +
      labs(color = group.by)
  } else {
    emb$group <- as.factor(emb$group)
    levels_vec <- levels(emb$group)
    # palette logic
    if (is.null(palette)) {
      pal_vec <- unname(pals::alphabet())
      pal_use <- setNames(pal_vec[seq_along(levels_vec)], levels_vec)
    } else {
      # Check if the user provided a named vector
      if (!is.null(names(palette))) {
        # Filter and reorder the palette to match levels exactly
        # This ensures cluster "a" always gets the color assigned to "a"
        pal_use <- palette[levels_vec]
        
        # Handle cases where some levels might be missing from the palette
        pal_use[is.na(pal_use)] <- na_color
      } else {
        # Fallback for unnamed palettes: match by position
        pal_vec <- palette
        if (length(pal_vec) < length(levels_vec))
          pal_vec <- rep(pal_vec, length.out = length(levels_vec))
        pal_use <- setNames(pal_vec[seq_along(levels_vec)], levels_vec)
      }
    }

    p <- ggplot(emb, aes(Dim1, Dim2, color = group)) +
      geom_fun(size = point_size, alpha = alpha) +
      scale_color_manual(values = pal_use, na.value = na_color) +
      labs(color = group.by)

    # optional labeling
    if (label) {
      centers <- aggregate(cbind(Dim1, Dim2) ~ group, data = emb, FUN = median)
      p <- p +
        ggrepel::geom_text_repel(
          data = centers,
          aes(x = Dim1, y = Dim2, label = group),
          size = label_size,
          color = "black",
          fontface = "bold",
          max.overlaps = Inf,
          show.legend = FALSE
        )
    }
  }

  # styling
  p <- p +
    theme_bw() +
    coord_fixed() +
    guides(
      color = if (is.numeric(emb$group))
        guide_colorbar(barwidth = 0.8, barheight = 8)
      else
        guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = if (legend) "right" else "none"
    )

  return(p)
}

umap_plots <- function( seu,
                        slide_name,
                        outdir = ".",
                        file_prefix = "Results_Report_",
                        reduction = "umapharmony",
                        cluster_col = "clusters_unsup_harmony",
                        slide_id_col = "slide_id",
                        region_col = "region",
                        condition_col = "condition",
                        slide_cls = brewer.pal(8, "Set2"),
                        region_cls = NULL,
                        condition_cls = c("T" = "red2", "N" = "thistle"),
                        shuffle = TRUE) {

    pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
        width = 12, height = 10)
    
    md_cols <- colnames(seu@meta.data)

    # Normalize cluster_col to character vector
    if (is.list(cluster_col)) {
        cluster_cols <- unlist(cluster_col)
    } else {
        cluster_cols <- as.character(cluster_col)
    }

    # DimPlots cluster columns
    for (cluster in cluster_cols) {
        if (cluster %in% md_cols) {
            p <- plot_embedding(seu, reduction = reduction, group.by = cluster, shuffle = shuffle) + 
                ggtitle(paste0("UMAP: ", reduction, ": ", cluster))
            print(p)
        }
        else {
            message("Skipping DimPlot for missing metadata column: ", cluster)
        }
    }

    # DimPlots standard columns
    if (slide_id_col %in% md_cols) {
        p <- plot_embedding(seu, reduction = reduction, group.by = slide_id_col, shuffle = shuffle, label = FALSE, palette = slide_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", slide_id_col))
        print(p)
    }

    if (region_col %in% md_cols) {
        p <- plot_embedding(seu, reduction = reduction, group.by = region_col, shuffle = shuffle, label = FALSE, legend = FALSE, palette = region_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", region_col))
        print(p)
    }
    if (condition_col %in% md_cols) {
        seu@meta.data[[condition_col]] <- factor(seu@meta.data[[condition_col]],
                                         levels = names(condition_cls))
        p <- plot_embedding(seu, reduction = reduction, group.by = condition_col, shuffle = shuffle, label = FALSE, palette = condition_cls) + 
            ggtitle(paste0("UMAP: ", reduction, ": ", condition_col))
        print(p)
    }

    dev.off()
    message("Summary plots saved for slide: ", slide_name)
}

  # Function to get cluster colors
  get_cluster_colors <- function(clusters, palette_spec = NULL) {
    n_clusters <- length(clusters)
    base_pal <- unname(pals::alphabet())

    if (is.null(palette_spec)) {
      # default alphabet with extension if needed
      if (n_clusters <= length(base_pal)) {
        cols <- base_pal[seq_len(n_clusters)]
      } else {
        cols <- grDevices::colorRampPalette(base_pal)(n_clusters)
      }
    } else if (is.character(palette_spec)) {
      # user-supplied palette vector
      if (length(palette_spec) < n_clusters) {
        cols <- c(
          palette_spec,
          grDevices::colorRampPalette(palette_spec)(n_clusters - length(palette_spec))
        )
      } else {
        cols <- palette_spec[seq_len(n_clusters)]
      }
    } else if (is.function(palette_spec)) {
      cols <- palette_spec(n_clusters)
    } else {
      cols <- grDevices::rainbow(n_clusters)
    }

    stats::setNames(cols, as.character(clusters))
  }

xy_plots <- function(seu,
                          slide_name,
                          outdir = ".",
                          file_prefix = "Spatial_Plots_",
                          cluster_col = "clusters_unsup_harmony",
                          shuffle = TRUE,
                          palettes = NULL,
                          show_legend = TRUE) {

  # Ensure cluster_col is a character vector
  cluster_cols <- as.character(unlist(cluster_col))

  # ensure palettes align with cluster_cols
  if (!is.null(palettes)) {
    if (length(palettes) != length(cluster_cols)) {
      stop("Length of 'palettes' does not match number of cluster columns. Stopping...")
    }
  }

  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
      width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)

  # ------------------------
  # Full tissue plots (one page per clustering)
  # ------------------------
    for (i in seq_along(cluster_cols)) {
    grp <- cluster_cols[i]
    palette_spec <- if (!is.null(palettes)) palettes[[i]] else NULL

    if (all(c("x_slide_mm", "y_slide_mm") %in% md_cols) && grp %in% md_cols) {
      all_clusters <- sort(unique(seu@meta.data[[grp]]))

      # Use global get_cluster_colors function
      cols <- get_cluster_colors(all_clusters, palette_spec)

      print(
        xyplot(
          grp,
          metadata = seu@meta.data,
          ptsize = 0.01,
          cls = cols,
          show_legend = show_legend
        ) +
          coord_fixed() +
          labs(title = paste("SpatialPlot:", grp))
      )
    } else {
      message("Skipping FeatureScatter for ", grp,
              ": required columns missing (x_slide_mm, y_slide_mm, or cluster_col)")
    }
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}

xy_plots_by_region <- function(
  seu,
  slide_name,
  outdir = ".",
  file_prefix = "XY_Zoom_",
  cluster_col = "clusters_unsup_harmony",
  region_col = "region",
  regions_per_page = 2,
  shuffle = TRUE,
  palettes = NULL,
  show_legend = TRUE
) {
  # Ensure cluster_col is a character vector
  cluster_cols <- as.character(unlist(cluster_col))

  # ensure palettes align with cluster_cols
  if (!is.null(palettes) && length(palettes) != length(cluster_cols)) {
    stop("Length of 'palettes' must match number of cluster columns")
  }
  
  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")), width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)
  if (!region_col %in% md_cols)
    stop("Specified region_col '", region_col, "' not found in metadata.")
  
  region_vals <- seu@meta.data[[region_col]]
  if (is.factor(region_vals)) {
    regions <- levels(region_vals)
  } else {
    regions <- mixedsort(unique(region_vals))
  }

  # helper: plot a single region for a given clustering column
  plot_region_cluster <- function(obj_region, clust, palette_spec) {
    region_md_cols <- colnames(obj_region@meta.data)
    if (!all(c("x_slide_mm", "y_slide_mm") %in% region_md_cols) || !(clust %in% region_md_cols)) {
      return(ggplot() + theme_void() + labs(title = paste(clust, "missing for region")))
    }

    clusters_plot <- sort(unique(obj_region@meta.data[[clust]]))
    cols <- get_cluster_colors(clusters_plot, palette_spec)

    xyplot(
      clust,
      metadata = obj_region@meta.data,
      ptsize = 0.01,
      cls = cols,
      show_legend = show_legend
    ) +
      coord_fixed() +
      labs(title = paste(clust, "\nRegion", unique(obj_region@meta.data[[region_col]])))
  }

  if (length(cluster_cols) == 1) {
    # single cluster column: possibly multiple regions per page
    clust <- cluster_cols[1]
    palette_spec <- if (!is.null(palettes)) palettes[[1]] else NULL

    region_plots <- lapply(regions, function(r) {
      obj_region <- subset(seu, cells = rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ]))
      plot_region_cluster(obj_region, clust, palette_spec)
    })

    # Remove NULLs
    region_plots <- Filter(Negate(is.null), region_plots)

    # arrange multiple regions per page
    for (i in seq(1, length(region_plots), by = regions_per_page)) {
      idx <- i:min(i + regions_per_page - 1, length(region_plots))
      print(wrap_plots(region_plots[idx], ncol = 2))
    }

  } else {
    # multiple cluster columns: one region per page, all clusters on same page
    for (r in regions) {
      cells <- rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ])
      if (length(cells) == 0) next
      obj_region <- subset(seu, cells = cells)

      # if palettes is NULL, create a list of NULLs matching cluster columns
      region_palettes <- palettes
      if (is.null(region_palettes)) {
        region_palettes <- vector("list", length(cluster_cols))
      }

      region_plots <- mapply(function(clust, pal) {
        plot_region_cluster(obj_region, clust, pal)
      }, cluster_cols, region_palettes, SIMPLIFY = FALSE)

      # Remove NULLs if any
      region_plots <- Filter(Negate(is.null), region_plots)
      if (length(region_plots) == 0) next

      print(wrap_plots(region_plots, ncol = 2))
    }
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}

spatial_plots <- function(seu,
                          slide_name,
                          outdir = ".",
                          file_prefix = "XY_Zoom_",
                          cluster_col = "clusters_unsup_harmony",
                          shuffle = TRUE,
                          show_legend = TRUE) {

  # Normalize cluster_col to character vector
  if (is.list(cluster_col)) {
    cluster_cols <- unlist(cluster_col)
  } else {
    cluster_cols <- as.character(cluster_col)
  }

  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")),
      width = 12, height = 10)

  md_cols <- colnames(seu@meta.data)
  regions <- mixedsort(unique(seu$region))

  # ------------------------
  # 2) Per-region plots (one page per region, all clusterings side by side)
  # ------------------------
  for (r in regions) {
    obj_region <- subset(seu, subset = region == r)
    region_md_cols <- colnames(obj_region@meta.data)

    region_plots <- lapply(cluster_cols, function(clust) {

      if (!all(c("x_slide_mm", "y_slide_mm") %in% region_md_cols) || !(clust %in% region_md_cols)) {
        return(ggplot() + theme_void() +
                 labs(title = paste(clust, "missing for region", r)))
      }

      all_clusters <- sort(unique(obj_region@meta.data[[clust]]))
      cols <- setNames(alphabet()[1:length(all_clusters)], all_clusters)

      xyplot(clust,
             metadata = obj_region@meta.data,
             ptsize = 0.01,
             cls = cols,
             show_legend = show_legend) +
        coord_fixed() +
        labs(title = paste(clust, "\n- Region", r))
    })

    # Combine all clusterings side by side for this region
    combined_region <- wrap_plots(region_plots, ncol = 2)
    print(combined_region)
  }

  dev.off()
  message("Summary plots saved for slide: ", slide_name)
}


## QC

summarize_qc_by <- function(metadata, group_col,
                            qc_vars = c("nCount_RNA", "nFeature_RNA", "Area.um2", "nCell", "nCount_negprobes", "nCount_falsecode"),
                                        digits = 0) {
  stopifnot(group_col %in% colnames(metadata))
  
  metadata %>%
    group_by(.data[[group_col]]) %>%
    summarize(across(all_of(qc_vars),
                     list(
                       mean = ~mean(.x, na.rm = TRUE),
                       median = ~median(.x, na.rm = TRUE),
                       min = ~min(.x, na.rm = TRUE),
                       max = ~max(.x, na.rm = TRUE),
                       range = ~diff(range(.x, na.rm = TRUE)),
                       sd = ~sd(.x, na.rm = TRUE),
                       variance = ~var(.x, na.rm = TRUE),
                       IQR = ~IQR(.x, na.rm = TRUE),
                       count = ~sum(!is.na(.x)),
                       total = ~sum(.x, na.rm = TRUE)
                     ),
                     .names = "{.col}_{.fn}"
    ),
    .groups = "drop") %>%
    pivot_longer(
      cols = -all_of(group_col),
      names_to = c("Variable", "Statistic"),
      names_sep = "_(?=[^_]+$)",  # split only at the LAST underscore
      values_to = "Value"
    ) %>%
    dplyr::mutate(Value = round(Value, digits))
}


calculate_qc_summaries <- function(metadata,
                                   slide_col = "slide_name",
                                   region_col = "region",
                                   fov_col = "fov",
                                  #  qc_vars = c("nCount_RNA", "nFeature_RNA", "Area.um2", "nCell", "nCount_negprobes", "nCount_falsecode"),
                                   qc_vars = c("nCount_RNA", "nFeature_RNA", "Area.um2"),
                                   digits = 2) {
  stats_list <- list(
    slide_level  = summarize_qc_by(metadata, slide_col, qc_vars, digits),
    region_level = summarize_qc_by(metadata, region_col, qc_vars, digits),
    fov_level    = summarize_qc_by(metadata, fov_col, qc_vars, digits)
  )

  # # compute SNR per group
  # if ("nCount_negprobes" %in% qc_vars && "nCount_RNA" %in% qc_vars) {
    
  #   stats_list <- lapply(stats_list, function(df) {
  #     df %>%
  #       group_by(.data[[names(df)[1]]]) %>%
  #       mutate(SNR = Value[Variable == "nCount_RNA" & Statistic == "mean"] /
  #                   Value[Variable == "nCount_negprobes" & Statistic == "mean"])
  #   })  
  # }
  return(stats_list)
}


create_qc_plots <- function(seu, batch_col = "slide_id", fov_col = "fov", 
                             count_col = "nCount_RNA", feature_col = "nFeature_RNA",
                             neg_col = "nCount_negprobes", fc_col = "nCount_falsecode",
                             pal = NULL) {

    md <- seu@meta.data
  
    # Determine number of negative and system control probes
    nNeg <- length(unique(seu[["negprobes"]]$counts@Dimnames[[1]]))  # total negative probes
    nFC  <- length(unique(seu[["falsecode"]]$counts@Dimnames[[1]]))  # total system control probes
    
    # Compute FOV-level QC summary
    qc_stats <- md %>%
        group_by(.data[[batch_col]], .data[[fov_col]]) %>%
        summarise(
        Number_Cells_Per_FOV = mean(nCell),
        Mean_Transcripts_Per_Cell_Per_FOV = mean(.data[[count_col]]),
        Mean_Unique_Transcripts_Per_Cell_Per_FOV = mean(.data[[feature_col]]),
        Total_Transcripts_Per_FOV = sum(.data[[count_col]]),
        Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV = sum(.data[[neg_col]]) / (unique(nCell) * nNeg) ,
        Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV = sum(.data[[fc_col]]) / (unique(nCell) * nFC) ,
        .groups = "drop"
        )
    
    # Default palette if not provided
    if (is.null(pal)) {
        n_batches <- length(unique(md[[batch_col]]))
        pal <- scales::hue_pal()(n_batches)
    }
    
    # Metrics to plot
    plot_metrics <- c(
        "Mean_Transcripts_Per_Cell_Per_FOV",
        "Mean_Unique_Transcripts_Per_Cell_Per_FOV",
        "Total_Transcripts_Per_FOV",
        "Number_Cells_Per_FOV",
        "Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV",
        "Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV"
    )
    
    # Loop over metrics and plot
    for (metric in plot_metrics) {
        # skip if metric not present
        if (!metric %in% colnames(qc_stats)) next
        
        p <- ggplot(qc_stats, aes_string(x = batch_col, y = metric, fill = batch_col)) +
        geom_violin(alpha = 0.3) +
        geom_jitter(width = 0.25, height = 0, size = 0.5) +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal) +
        ylab(metric) +
        xlab("Batch / Slide") +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5)
        )
        
        print(p)
    }
    }

# Function: QC analysis and PDF report per Seurat object
# Args:
#   seu: Seurat object
#   slide_name: name of the slide for labeling outputs
#   outdir: directory to save PDF reports
# Returns:
#   Nothing, saves PDF report
qc_report <- function(seu, slide_name, outdir = ".", file_prefix = "QC_Report_", filter_log = NULL, count_column = "nCount_RNA", region_col = "region", flag_col = "qcCellsFlagged", boxplots = TRUE, fov_plots = TRUE) {
  pdf(file.path(outdir, paste0(file_prefix, slide_name, ".pdf")), width = 14, height = 10)

  # Tissue plot from metadata
  md <- seu@meta.data

  # If filter_log is provided, print the row for this slide
  if (!is.null(filter_log) && slide_name %in% filter_log$Slide) {
    slide_log <- filter_log[filter_log$Slide == slide_name, ]
    log_text <- paste(
      paste(names(slide_log), slide_log, sep = ": "), collapse = "\n"
    )
    grid::grid.newpage()
    grid::grid.text(
      log_text,
      x = 0.5, y = 0.5,
      gp = grid::gpar(fontsize = 12), just = "center"
    )
  }

  if (all(c("x_slide_mm", "y_slide_mm") %in% colnames(md))) {
    tissue_plot <- xyplot(region_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                      metadata = md, 
                      alphasize = 1, show_legend = FALSE, show_labels = TRUE,
                      label_color = "black",
                      label_repel = TRUE) +
                      ggtitle(paste("Tissue Layout -", slide_name))
    print(tissue_plot)
    region_colors <- ggplot2::ggplot_build(tissue_plot)$plot$scales$get_scales("colour")$palette.cache
  }

  if ("nCount_RNA" %in% colnames(md) & "nFeature_RNA" %in% colnames(md)) {
    md$log10_nCount_RNA <- log10(md$nCount_RNA + 1)
    md$log10_nFeature_RNA <- log10(md$nFeature_RNA + 1)
    
    tissue_plot_nCount <- xyplot("log10_nCount_RNA", x_column = "x_slide_mm", y_column = "y_slide_mm", 
                      metadata = md, 
                      alphasize = 1, show_legend = TRUE, show_labels = FALSE,
                      continuous_palette = function(n) viridis::viridis(n, option = "inferno")) +
                      ggtitle(paste("log10(nCountRNA + 1) -", slide_name))
    print(tissue_plot_nCount)

    # tissue_plot_nFeature <- xyplot("log10_nFeature_RNA", x_column = "x_slide_mm", y_column = "y_slide_mm", 
    #                   metadata = md, 
    #                   alphasize = 1, show_legend = TRUE, show_labels = FALSE,
    #                   continuous_palette = function(n) viridis::viridis(n, option = "inferno")) +
    #                   ggtitle(paste("log10(FeatureRNA + 1) -", slide_name))
    # print(tissue_plot_nFeature)
  }
  
  if (flag_col %in% colnames(md)) {
    flag_cls <- c("TRUE" = "red4", "FALSE" = "gray")

    tissue_plot_qc <- xyplot(flag_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                   cls = flag_cls, metadata = md, ptsize = 0.25,
                   alphasize = 1, show_legend = TRUE, show_labels = FALSE) +
                   ggtitle(paste(flag_cls,": AtoMX QC Flags"))
    print(tissue_plot_qc)
  }

  # Tissue plot with condition
  condition_col <- "condition"

  if (condition_col %in% colnames(md)) {
    condition_cls <- c("T" = "red2", "N" = "thistle")

    tissue_plot_cond <- xyplot(condition_col, x_column = "x_slide_mm", y_column = "y_slide_mm", 
                   cls = condition_cls,metadata = md, ptsize = 0.25,
                   alphasize = 1, show_legend = TRUE, show_labels = FALSE) +
                   ggtitle(paste("Tissue plot ", condition_col))
    print(tissue_plot_cond)
  }
  
  # Violin plots with median bars
  vp1 <- VlnPlot(seu, features = count_column, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
    ggtitle(paste(count_column," per Cell"))

  count_column_2 <- "nFeature_RNA"
  vp2 <- VlnPlot(seu, features = count_column_2, pt.size = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
    ggtitle(paste(count_column_2," per Cell"))
  
  print(vp1 + vp2 + plot_layout(ncol = 2) + plot_annotation(title = paste("QC Metrics Violin Plots -", slide_name)))

  # Summary stats 
  qc_stats <- calculate_qc_summaries(
      metadata = md,
      slide_col = "slide_id",
      region_col = region_col,
      fov_col = "fov",
      digits = 2
    )

  ## Slide Level
  slide_table <- qc_stats$slide_level %>%
    pivot_wider(names_from = Statistic, values_from = Value)

  grid::grid.newpage()
  grid::grid.text("Slide-Level QC Summary", y = 0.95, gp = grid::gpar(fontsize = 16, fontface = "bold"))

  qc_table_grob <- gridExtra::tableGrob(
      slide_table,
      rows = NULL,
      theme = gridExtra::ttheme_default(
        core = list(fg_params = list(cex = 0.7)),
        colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
      )
    )

  grid::grid.draw(qc_table_grob)

  # Signal-to-noise ratio
  if ("negprobes" %in% names(seu)) {
    snr <- mean(Matrix::colMeans(seu[["RNA"]]$counts)) / mean(Matrix::colMeans(seu[["negprobes"]]$counts))
  } else {
    snr <- NA
  }

  # Flagged cells at the bottom
  if (flag_col %in% colnames(md)) {
    flagged_count <- sum(md[[flag_col]], na.rm = TRUE)
    total_cells <- ncol(seu)
    flagged_text <- paste0("SNR: ", round(snr, 2),
                          "\nTotal cells: ", total_cells,
                          "\nFlagged cells: ", flagged_count,
                          " (", round(100 * flagged_count / total_cells, 2), "%)")
    grid::grid.text(
      flagged_text,
      x = 0.5, y = 0.3,  # lower vertical position
      gp = grid::gpar(fontsize = 14)
    )
  }

  # Violin plots split by condition
  if (condition_col %in% colnames(seu@meta.data)) {
    p_condition <- VlnPlot(seu, features = count_column, pt.size = 0, group.by = condition_col, cols = condition_cls) +
      ggtitle(paste(count_column," per Cell by Condition")) +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_condition)
  }

  # Violin plots split by region
  if (region_col %in% colnames(seu@meta.data)) {
    p_region <- VlnPlot(seu, features = count_column, pt.size = 0, group.by = region_col, cols = region_colors) +
      ggtitle(paste(count_column," per Cell by Region")) +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_region)
  }

  if (boxplots) {
    if (region_col %in% colnames(seu@meta.data)) {
      p_region <- ggplot(seu@meta.data, aes_string(x = region_col, y = count_column, fill = region_col)) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = region_colors) +
        labs(title = paste(count_column," per Region"), x = "Region", y = count_column) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_region)
    }
  }

  if (region_col %in% colnames(seu@meta.data)) {
    p_region <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0, group.by = region_col, cols = region_colors) +
      ggtitle("nFeature_RNA per Cell by Region") +
      stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_region)
  }
  
  if (boxplots) {
    if (region_col %in% colnames(seu@meta.data)) {
      p_region <- ggplot(seu@meta.data, aes_string(x = region_col, y = "nFeature_RNA", fill = region_col)) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = region_colors) +
        labs(title = "nFeature_RNA per Region", x = "Region", y = "nFeature_RNA") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_region)
    }
  }
    
    # Bar plot of cell counts per region
    cell_counts_df <- as.data.frame(table(seu@meta.data[[region_col]]))
    colnames(cell_counts_df) <- c("Region", "Cell_Count")
    
    p_counts <- ggplot(cell_counts_df, aes(x = Region, y = Cell_Count, fill = Region)) +
      geom_bar(stat = "identity", color = "black") +
      labs(title = "Number of Cells per Region", x = "Region", y = "Cell Count") +
      scale_fill_manual(values = region_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(legend.position = "none")
    print(p_counts)
  
  # Cell size distribution
  if ("Area.um2" %in% colnames(md)) {
    area_plot <- ggplot(md, aes(x = Area.um2)) +
      geom_histogram(bins = 100, fill = "grey70", color = "black") +
      geom_vline(xintercept = c(40, 600), color = "red", linetype = "dashed", linewidth = 0.5) +
      theme_minimal() +
      ggtitle(paste("Distribution of Cell Sizes (µm²) -", slide_name)) +
      xlab("Area (µm²)") + ylab("Cell count")
    print(area_plot)
  }

if (boxplots) {
  if (region_col %in% colnames(seu@meta.data)) {
  p_area <- ggplot(seu@meta.data, aes_string(x = region_col, y = "Area.um2", fill = region_col)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = region_colors) +
    labs(title = "Cell Area per Region", x = "Region", y = "Area (µm²)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none"
    )
  print(p_area)
}
}

## FOV Border proximity
if ("fov_border_proximity" %in% colnames(seu@meta.data)) {
    tab <- table(seu@meta.data$fov_border_proximity, useNA = "no")
    required <- c("border", "interior")
    tab_req <- tab[required]
    tab_req[is.na(tab_req)] <- 0
    proximity_cls <- c("border" = "#06b8aa", "interior" = "#6e6edf")
    if (all(tab_req > 0)) {
      p_area_fov_border <- ggplot(seu@meta.data, aes_string(x = "fov_border_proximity", y = "Area.um2", fill = "fov_border_proximity")) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = proximity_cls) +
        labs(title = "Cell Area per FOV Border Proximity", x = "FOV Border Proximity", y = "Area (µm²)") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_area_fov_border)

      p_count_fov_border <- ggplot(seu@meta.data, aes_string(x = "fov_border_proximity", y = "nCount_RNA", fill = "fov_border_proximity")) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = proximity_cls) +
        labs(title = "nCount_RNA per FOV Border Proximity", x = "FOV Border Proximity", y = "nCount_RNA") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none"
        )
      print(p_count_fov_border)

      n_border   <- as.integer(tab["border"])
      n_interior <- as.integer(tab["interior"])
      n_total    <- n_border + n_interior

      pct_border   <- 100 * n_border / n_total
      pct_interior <- 100 * n_interior / n_total

      grid::grid.newpage()
      txt <- paste0(
        "FOV Border Proximity Summary\n\n",
        "Border cells:   ", n_border,   " (", sprintf("%.1f", pct_border),   "%)\n",
        "Interior cells: ", n_interior, " (", sprintf("%.1f", pct_interior), "%)\n",
        "Total cells:    ", n_total
      )

      grid::grid.text(
        txt,
        x = 0.05,
        y = 0.95,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 10)
      )

      df <- seu@meta.data[, c("fov", "fov_border_proximity")]
      df_sum <- as.data.frame(
        table(df$fov, df$fov_border_proximity),
        stringsAsFactors = FALSE
      )
      colnames(df_sum) <- c("fov", "proximity", "n")

      df_sum <- df_sum |>
        dplyr::group_by(fov) |>
        dplyr::mutate(pct = 100 * n / sum(n)) |>
        dplyr::ungroup()

      p_fov_pct <- ggplot(
        df_sum,
        aes(x = fov, y = pct, fill = proximity)
      ) +
        geom_col(width = 0.9) +
        scale_fill_manual(values = proximity_cls) +
        labs(
          title = "Border vs Interior Composition per FOV",
          x = "FOV",
          y = "Percentage"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1)
        )

      print(p_fov_pct)

      df_border <- df_sum[df_sum$proximity == "border", ]

      p_fov_abs <- ggplot(
        df_border,
        aes(x = fov, y = n)
      ) +
        geom_col(fill = proximity_cls["border"], width = 0.9) +
        labs(
          title = "Absolute Number of Border Cells per FOV",
          x = "FOV",
          y = "Number of Border Cells"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1)
        )

      print(p_fov_abs)
  
}
}


## FOV Level
if (fov_plots) {
create_qc_plots(seu, batch_col = "slide_id", fov_col = "fov", count_col = count_column, feature_col = "nFeature_RNA", 
                  neg_col = "nCount_negprobes", fc_col = "nCount_falsecode")
}

  
  dev.off()
}

# Function: Run QC for list of Seurat objects
qc_all_slides <- function(seurat_list, outdir = ".", file_prefix = "QC_Report_", filter_log = NULL, count_column = "nCount_RNA", region_col = "region", flag_col = "qcCellsFlagged", boxplots = TRUE, fov_plots = TRUE) {
  for (nm in names(seurat_list)) {
    message("Running QC for: ", nm)
    qc_report(seurat_list[[nm]], nm, outdir, file_prefix = file_prefix, filter_log = filter_log, count_column = count_column, region_col = region_col, flag_col = flag_col, boxplots = boxplots, fov_plots = fov_plots)
  }
}


filter_seurat_objects <- function(seurat_list, thresholds_list = NULL, low_thres = 0.025, high_thres = 0.975, area_upper = 600, area_lower = 40, flag_col = "remove_flagged_cells", region_col = "region", min_cells_per_region = 500, filter_for_fov_border_proximity = TRUE, fov_border_col = "fov_border_proximity") {
  filtered_list <- list()
  log_list <- list()
  
  for (slide_name in names(seurat_list)) {
    message("Filtering: ", slide_name)
    seu <- seurat_list[[slide_name]]
    # seu@images <- list()
    # seu <- suppressWarnings(suppressMessages(UpdateSeuratObject(seu)))
    # seu <- UpdateSeuratObject(seu)
    total_cells <- ncol(seu)

    # Count flagged cells before any filtering
    message("Counting flagged cells before filtering ", slide_name)
    flagged_before_total <- if (flag_col %in% colnames(seu@meta.data)) {
      sum(seu@meta.data[[flag_col]], na.rm = TRUE)
    } else 0
    
    # Determine thresholds
    message("Determining thresholds ", slide_name)
    counts_vector <- seu$nCount_RNA
    if (!is.null(thresholds_list) && slide_name %in% names(thresholds_list)) {
      lower <- thresholds_list[[slide_name]]$lower
      upper <- thresholds_list[[slide_name]]$upper
    } else {
      lower <- quantile(counts_vector, low_thres)
      upper <- quantile(counts_vector, high_thres)
    }
    
    # Quantile / threshold filtering
    bottom_removed <- sum(counts_vector <= lower)
    top_removed <- sum(counts_vector >= upper)
    seu_filtered <- subset(seu, subset = nCount_RNA > lower & nCount_RNA < upper)
    
    remaining_after_counts <- ncol(seu_filtered)
    
    bottom_pct <- bottom_removed / total_cells * 100
    top_pct <- top_removed / total_cells * 100

    # Count flagged cells after threshold filtering
    message("Filtering flagged cells for ", slide_name)
    flagged_after_threshold <- if (flag_col %in% colnames(seu_filtered@meta.data)) {
      sum(seu_filtered@meta.data[[flag_col]], na.rm = TRUE)
    } else 0
    
    # Filter flagged cells
    # seu_filtered <- subset(seu_filtered, subset = remove_flagged_cells == FALSE)
    keep_cells <- rownames(seu_filtered@meta.data)[!seu_filtered@meta.data[[flag_col]]]
    seu_filtered <- subset(seu_filtered, cells = keep_cells)
    flagged_removed <- flagged_after_threshold
    flagged_pct <- flagged_removed / total_cells * 100

    ## ---- FOV border proximity filtering
    border_removed <- 0
    if (
      filter_for_fov_border_proximity &&
      fov_border_col %in% colnames(seu_filtered@meta.data)
    ) {
      border_removed <- sum(
        seu_filtered@meta.data[[fov_border_col]] == "border",
        na.rm = TRUE
      )

      if (border_removed > 0) {
        keep_cells <- rownames(seu_filtered@meta.data)[
          seu_filtered@meta.data[[fov_border_col]] != "border"
        ]
        seu_filtered <- subset(seu_filtered, cells = keep_cells)
      }
    }

    # Filter by Area.um2
    area_removed <- 0
    if ("Area.um2" %in% colnames(seu_filtered@meta.data)) {
      area_vals <- seu_filtered@meta.data$Area.um2
      keep_cells <- rownames(seu_filtered@meta.data)[
        is.na(area_vals) |
        (area_vals >= area_lower & area_vals <= area_upper)
      ]
      area_removed <- ncol(seu_filtered) - length(keep_cells)
      if (area_removed > 0) {
        seu_filtered <- subset(seu_filtered, cells = keep_cells)
      }
    }

    # Region-based Filtering
    regions_removed_count <- 0
    if (region_col %in% colnames(seu_filtered@meta.data)) {
      # Count cells currently in each region
      region_counts <- table(seu_filtered@meta.data[[region_col]])
      # Identify which regions meet the minimum threshold
      keep_regions <- names(region_counts[region_counts >= min_cells_per_region])
      # Calculate how many cells are in the 'dropped' regions for the log
      regions_removed_count <- sum(region_counts[!(names(region_counts) %in% keep_regions)])
      # Perform the subsetting
      seu_filtered <- subset(seu_filtered, 
                             cells = rownames(seu_filtered@meta.data)[seu_filtered@meta.data[[region_col]] %in% keep_regions])
      seu_filtered@meta.data[[region_col]] <- droplevels(as.factor(seu_filtered@meta.data[[region_col]]))

      message(paste0("  Regions dropped: ", length(region_counts) - length(keep_regions), 
                     " (Total cells removed from these regions: ", regions_removed_count, ")"))
    }

    remaining_cells <- ncol(seu_filtered)
    
    # Store filtered object
    filtered_list[[slide_name]] <- seu_filtered

    total_removed <- bottom_removed +
                     top_removed +
                     flagged_removed +
                     border_removed +
                     area_removed +
                     regions_removed_count
    
    #LOG
    log_list[[slide_name]] <- data.frame(
      Slide = slide_name,
      Total_Cells = total_cells,
      Bottom_Threshold = lower,
      Bottom_Removed = bottom_removed,
      Bottom_Percent = round(bottom_pct, 2),
      Top_Threshold = upper,
      Top_Removed = top_removed,
      Top_Percent = round(top_pct, 2),
      Flagged_Before_Threshold = flagged_before_total,
      Flagged_After_Threshold = flagged_after_threshold,
      Flagged_Removed = flagged_removed,
      Flagged_Percent = round(flagged_pct, 2),
      Border_Removed = border_removed,
      Area_Upper_Threshold = area_upper,
      Area_Lower_Threshold = area_lower,
      Area_Removed = area_removed,
      Remaining_Cells = remaining_cells,
      Remaining_Percent = round(remaining_cells / total_cells * 100, 2),
      Cells_Removed = total_removed,
      Cells_Removed_Percent = round((total_removed) / total_cells * 100, 2),
      stringsAsFactors = FALSE
    )
    
    message(paste0("Filtered slide: ", slide_name,
                   "  | Total cells: ", total_cells,
                   "  | Cells removed: ", total_removed,
                   "  | Cells removed (calc): ", total_cells - remaining_cells,
                   " | Remaining cells: ", remaining_cells))
  }
  
  filter_log <- do.call(rbind, log_list)
  rownames(filter_log) <- NULL
  
  return(list(filtered_objects = filtered_list, filter_log = filter_log))
}


## BATCH CORRECTION
apply_scPearsonPCA <- function(seu, 
                              nfeatures = 3000,
                              resolution = 0.8,
                              slot_names = list(
                                                pca = "pearsonpca",
                                                umap = "pearsonumap",
                                                graph = "pearsongraph",
                                                clusters = "pearson_clusters"
                                              )
  ) {
    message("Calculating total counts and gene frequencies...")
    counts_temp <- Seurat::GetAssayData(seu, assay = "RNA", layer = "counts")

    tc <- Matrix::colSums(counts_temp) ## total counts per cell (across all genes)
    genefreq <- scPearsonPCA::gene_frequency(counts_temp) ## gene frequency (across all cells)

    if (!sum(genefreq)==1) {
      message("Gene frequencies do not sum to 1.")
      return(NULL)
        } else {
        message("Gene frequencies sum to 1. Check passed.")
        }

    message("Finding highly variable genes...")    
    seu <- Seurat::FindVariableFeatures(seu, nfeatures = nfeatures)
    hvgs <- Seurat::VariableFeatures(seu)

    message("Performing scPearsonPCA...")   
    pcaobj <- sparse_quasipoisson_pca_seurat(seu[["RNA"]]$counts[hvgs,]
                               ,totalcounts = tc
                               ,grate = genefreq[hvgs]
                               ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                               ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                               ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
                               )

    message("Making umap...")
    umapobj <- scPearsonPCA::make_umap(pcaobj)
    seu[[slot_names$pca]] <- pcaobj$reduction.data
    seu[[slot_names$umap]] <- umapobj$ump  ## umap
    seu[[slot_names$graph]] <- Seurat::as.Graph(umapobj$grph) ## nearest neighbors / adjacency matrix used for unsupervised clustering
    
    message("Finding clusters...")
    seu <- Seurat::FindClusters(seu, graph = slot_names$graph, resolution = resolution)
    cluster_colname <- paste0(slot_names$clusters, "_res", resolution)
    seu@meta.data[[cluster_colname]] <- seu@meta.data$seurat_clusters

    list(
    tc = tc,
    hvgs = hvgs,
    seu = seu
  )
}

## (minor side-note), this helper function `make_umap` calls the same package 
## and function as Seurat::RunUMAP, but also returns the 
## 'nearest-neighbor graph' used in the UMAP algorithm.
## Using the same nearest neighbor graph for clustering and UMAP, 
## can give better concordance between UMAP and unsupervised clusters.

make_umap <- function(pcaobj, min_dist=0.01, n_neighbors=30, metric="cosine",key ="UMAP_" ){
  ump <- 
    uwot::umap(pcaobj$reduction.data@cell.embeddings
               ,n_neighbors = n_neighbors
               ,nn_method = "annoy"
               ,metric = metric
               ,min_dist = min_dist
               ,ret_extra = c("fgraph","nn")
               ,verbose = TRUE)
  
  umpgraph <- ump$fgraph
  dimnames(umpgraph) <- list(rownames(ump$nn[[1]]$idx), rownames(ump$nn[[1]]$idx))
  colnames(ump$embedding) <- paste0(key, c(1,2)) 
  ump <- Seurat::CreateDimReducObject(embeddings = ump$embedding, key = key)
  return(list(grph = umpgraph
              ,ump = ump))
}


# From: https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/blob/Main/_code/HieraType/R/utils.R
#' Normalize cell expression in a raw counts matrix by their totalcounts
#' param counts_matrix a cells x genes matrix of raw counts to be normalized
#' param tc optional vector of totalcounts.  Useful if providing a counts_matrix based on a subset of genes.  
#' If not provided, totalcounts per cells is taken to be the rowSums of the counts matrix.
totalcount_norm <- function(counts_matrix, tc = NULL){
  if(is.null(tc)) tc <- Matrix::rowSums(counts_matrix)
  scale.factor <- mean(tc)
  tc[tc==0] <- 1
  return(Matrix::Diagonal(x = scale.factor/tc, names = TRUE) %*% counts_matrix)
}

## Functions relevant for Cell Typing with InSituType

plot_composition <- function(seu, cluster_col, split_by, type = "relative", palette = NULL) {
  
  # 1. Extract and aggregate data
  df <- seu@meta.data %>%
    dplyr::group_by(!!rlang::sym(cluster_col), !!rlang::sym(split_by)) %>%
    dplyr::tally() %>%
    dplyr::ungroup()
  
  colnames(df) <- c("Cluster", "SplitVar", "Count")
  
  # 2. Determine scaling and labels
  # 'fill' scales bars to 1 (100%), 'stack' keeps raw counts
  pos <- if (type == "relative") "fill" else "stack"
  label_y <- if (type == "relative") "Proportion of Cells" else "Number of Cells"

  # 3. Setup Palette (Default to alphabet if NULL)
  if (is.null(palette)) {
    unique_items <- as.character(unique(df$SplitVar))
    n_items <- length(unique_items)
    
    if (n_items <= 26) {
      # Use alphabet directly if small enough
      pal_vec <- as.character(pals::alphabet(n = n_items))
    } else {
      # Interpolate alphabet colors to handle 26+ items
      pal_vec <- colorRampPalette(as.character(pals::alphabet()))(n_items)
    }
    palette <- setNames(pal_vec, unique_items)
  }

  p <- ggplot(df, aes(x = Cluster, y = Count, fill = SplitVar)) +
    geom_bar(stat = "identity", position = pos) +
    scale_fill_manual(values = palette) + # Use the alphabet palette
    theme_bw() +
    labs(
      title = paste("Composition of", cluster_col, "by", split_by),
      subtitle = paste("Mode:", type),
      x = "Cluster",
      y = label_y,
      fill = split_by
    ) +
    theme(
      # Updated to 90 degrees and adjusted justification to align with ticks
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank() # Clean up the X-axis grid for bar charts
    )
  
  return(p)
}


# CT_QC_plot:
# This function generates a comprehensive PDF report of cluster-level quality control (QC) plots for cell typing validation.
#
# The PDF includes:
#   - Barplot of cell counts per cluster.
#   - Violin plots for QC features across clusters (e.g., marker expression, gene/cell counts).
#   - Optional heatmap and flightpath/trajectory plots if InSituType results provided.
#   - Embedding plots of clusters, annotations, and study IDs for visual inspection.

CT_QC_plot <- function(seu, cluster_col, cluster_pal=NULL, annotation_col = NULL, annotation_pal = NULL, IST_obj = NULL, out_dir = ".", reduction = NULL, split_by_col = "study_id") {
    
    pdf(file.path(out_dir, paste0("Cluster_QC_", cluster_col, ".pdf")), width = 10, height = 10)

    print(count_per_ct(seu, cluster_col=cluster_col, cluster_pal = cluster_pal))
    sp_plot <- xyplot(cluster_col,
                   x_column = "x_slide_mm",
                   y_column = "y_slide_mm",
                   cls = cluster_cls,
                   metadata = seu@meta.data)
    print(sp_plot)
    print(feature_vlns(seu, cluster_col=cluster_col, features = c("Mean.PanCK", "Mean.CD45", "nCount_RNA", "nFeature_RNA"), cluster_pal = cluster_pal))


    if (!is.null(IST_obj)) {
        heatmap(sweep(IST_obj$profiles, 1, pmax(apply(IST_obj$profiles, 1, max), .2), "/"), scale = "none",
            main = "Cluster mean expression profiles")
        fp_layout(IST_obj, cluster_pal = cluster_pal)
        print(flightpath_plot(flightpath_result = NULL, insitutype_result = IST_obj, col = cluster_pal[IST_obj$clust]))
    }

    if (reduction %in% names(seu@reductions)) {
        p1 <- plot_embedding(
            seu,
            reduction = reduction,
            group.by = cluster_col,
            label = TRUE,
            palette = cluster_pal,
            legend = TRUE
            ) 
        print(p1)

        if (!is.null(annotation_col)) {
            p2 <- plot_embedding(
                seu,
                reduction = reduction,
                group.by = annotation_col,
                label = TRUE,
                palette = annotation_pal,
                legend = TRUE
                ) 
            print(p2)
        }

        p3 <- plot_embedding(
            seu,
            reduction = reduction,
            group.by = split_by_col,
            label = TRUE,
            palette = NULL,
            legend = TRUE
            ) 
        print(p3)
    }

    print(plot_composition(seu, cluster_col = cluster_col, split_by = split_by_col, type = "relative"))
    print(plot_composition(seu, cluster_col = cluster_col, split_by = split_by_col, type = "absolute"))
    print(plot_composition(seu, cluster_col = split_by_col, split_by = cluster_col, type = "relative", palette = cluster_pal))
    print(plot_composition(seu, cluster_col = split_by_col, split_by = cluster_col, type = "absolute", palette = cluster_pal))

    dev.off()
    }

count_per_ct <- function(seu, cluster_col, xlab = "Cell Type",
                         ylab = "Cell Count", title = "Cell Count per Cell Type",
                         cluster_pal = NULL) {

    df <- seu@meta.data %>%
        dplyr::count(celltype = .data[[cluster_col]]) %>%
        dplyr::arrange(desc(n))

    p <- ggplot(df, aes(x = reorder(celltype, -n), y = n, fill = celltype)) +
        geom_bar(stat = "identity") +
        labs(x = xlab, y = ylab, title = title) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")

    if (!is.null(cluster_pal)) {
        p <- p + scale_fill_manual(values = cluster_pal)}
    p
}

feature_vlns <- function(seu, cluster_col, features = c("Mean.PanCK", "Mean.CD45"), cluster_pal = NULL) {
    levels_ct <- levels(seu[[cluster_col]][,1])
    cols_vec <- cluster_pal[levels_ct]

    plots <- lapply(features, function(feature) {
        VlnPlot(
            seu,
            features = feature,
            group.by = cluster_col,
            pt.size = 0,
            cols = cols_vec
        ) +
            ggtitle(paste0(feature, " by ", cluster_col)) +
            theme_bw() +
            theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)
            ) +
            coord_cartesian(clip = "off")
        })

    plots
}

plot_celltype_highlights <- function(seu, 
                                     cluster_col, 
                                     cluster_cls, 
                                     out_dir,
                                     ptsize = 0.01,
                                     background_alpha = 0.3,  # Alpha for background celltypes
                                     highlight_alpha = 1,      # Alpha for highlighted celltype
                                     x_column = "x_slide_mm",
                                     y_column = "y_slide_mm",
                                     show_legend = TRUE) {
  
  # Get unique celltypes from the cluster column
  celltypes <- unique(seu@meta.data[[cluster_col]])
  celltypes <- celltypes[!is.na(celltypes)]
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Loop through each celltype
  for (ct in celltypes) {
    message(paste0("Plotting highlighted celltype: ", ct))
    
    # Use xyplot's plotfirst functionality to control layering
    # First plot background with low alpha, then overlay the celltype of interest
    p <- xyplot(
      cluster_col,
      x_column = x_column,
      y_column = y_column,
      cls = cluster_cls,  # Use original colors
      clusters = NULL,
      metadata = seu@meta.data,
      ptsize = ptsize,
      plotfirst = ct,           # Plot this celltype separately
      plotfirst_on_top = TRUE,  # Plot it on top
      alphasize = background_alpha,  # Background gets low alpha
      show_legend = show_legend,
      coord_equal = TRUE,
      continuous_palette = function(n) viridis::viridis(n, option = "plasma"),
      aes_mappings = list(size = NULL, shape = NULL, alpha = NULL),
      order = NULL,
      na_color = "black",
      theme = ggplot2::theme_bw(),
      show_labels = FALSE,
      label_size = 4,
      label_color = "black",
      label_fontface = "bold",
      label_method = "median",
      label_repel = TRUE,
      label_max_overlaps = 50
    )
    
    # The plotfirst layer is plotted with alpha=1, but we want to make it black
    # So we need to add another layer on top
    highlight_data <- seu@meta.data[seu@meta.data[[cluster_col]] == ct & !is.na(seu@meta.data[[cluster_col]]), ]
    
    if (nrow(highlight_data) > 0) {
      p <- p + 
        ggplot2::geom_point(
          data = highlight_data,
          ggplot2::aes(x = !!rlang::sym(x_column), y = !!rlang::sym(y_column)),
          color = "black",  # Highlighted celltype in black
          size = ptsize,
          alpha = highlight_alpha
        ) +
        ggplot2::ggtitle(paste0("Highlighted: ", ct))
    } else {
      p <- p + ggplot2::ggtitle(paste0("Highlighted: ", ct, " (no cells found)"))
    }
    
    # Save the plot
    safe_filename <- gsub("[^a-zA-Z0-9_-]", "_", ct)
    output_file <- file.path(out_dir, paste0("highlight_", safe_filename, ".png"))
    
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = 12,
      height = 10,
      dpi = 300
    )
    
    message(paste0("  Saved: ", output_file))
  }
  
  message(paste0("Completed all ", length(celltypes), " celltype highlight plots"))
  
  return(invisible(NULL))
}

plot_celltypes_by_region <- function(
  seu,
  cluster_col,
  cluster_cls,
  region_col = "region",
  out_dir,
  mode = c("highlight", "all_celltypes"),
  celltypes = NULL,
  regions = NULL,
  regions_per_page = 4,
  ptsize = 0.01,
  background_alpha = 0.3,
  x_column = "x_slide_mm",
  y_column = "y_slide_mm",
  show_legend = TRUE,
  include_summary = TRUE
) {
  
  require(patchwork)
  require(gtools)
  
  mode <- match.arg(mode)
  
  # Validate region_col exists
  if (!region_col %in% colnames(seu@meta.data)) {
    stop("Specified region_col '", region_col, "' not found in metadata.")
  }
  
  # Get all celltypes
  all_celltypes <- unique(seu@meta.data[[cluster_col]])
  all_celltypes <- all_celltypes[!is.na(all_celltypes)]
  
  if (is.null(celltypes)) {
    celltypes <- all_celltypes
  } else {
    # Validate requested celltypes exist
    missing <- setdiff(celltypes, all_celltypes)
    if (length(missing) > 0) {
      warning("Celltypes not found in data: ", paste(missing, collapse = ", "))
      celltypes <- intersect(celltypes, all_celltypes)
    }
  }
  
  # Get all regions using mixedsort
  region_vals <- seu@meta.data[[region_col]]
  if (is.null(regions)) {
    if (is.factor(region_vals)) {
      regions <- levels(region_vals)
    } else {
      regions <- gtools::mixedsort(unique(region_vals))
    }
  } else {
    regions <- gtools::mixedsort(regions)
  }
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Determine grid layout based on regions_per_page
  if (regions_per_page == 2) {
    ncol_grid <- 1
    nrow_grid <- 2
  } else if (regions_per_page == 4) {
    ncol_grid <- 2
    nrow_grid <- 2
  } else {
    stop("regions_per_page must be either 2 or 4")
  }
  
  # ========================================
  # MODE: all_celltypes - One PDF for all regions
  # ========================================
  if (mode == "all_celltypes") {
    message("Mode: all_celltypes - creating one PDF with all regions")
    
    pdf_path <- file.path(out_dir, "all_celltypes_regions.pdf")
    pdf(pdf_path, width = 12, height = 10)
    
    # Optional: Create summary page showing all regions
    if (include_summary) {
      message("  Creating summary page...")
      
      summary_plot <- xyplot(
        cluster_col,
        x_column = x_column,
        y_column = y_column,
        cls = cluster_cls,
        clusters = NULL,
        metadata = seu@meta.data,
        ptsize = ptsize,
        alphasize = 1,
        show_legend = show_legend,
        coord_equal = TRUE,
        na_color = "black",
        theme = ggplot2::theme_bw(),
        show_labels = FALSE
      ) +
        ggplot2::ggtitle("Summary: All celltypes (all regions)") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
      
      print(summary_plot)
    }
    
    # Create plots for each region
    region_plots <- list()
    
    for (r in regions) {
      message(paste0("  Processing region: ", r))
      
      # Subset to this region
      region_cells <- rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ])
      
      if (length(region_cells) == 0) {
        # Create empty plot for missing region
        p <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste0("Region ", r, "\n(No cells)")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        region_plots[[as.character(r)]] <- p
        next
      }
      
      obj_region <- subset(seu, cells = region_cells)
      region_metadata <- obj_region@meta.data
      
      # Check if coordinates exist
      if (!all(c(x_column, y_column) %in% colnames(region_metadata))) {
        p <- ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste0("Region ", r, "\n(Missing coordinates)")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        region_plots[[as.character(r)]] <- p
        next
      }
      
      # Create plot with all celltypes
      p <- xyplot(
        cluster_col,
        x_column = x_column,
        y_column = y_column,
        cls = cluster_cls,
        clusters = NULL,
        metadata = region_metadata,
        ptsize = ptsize,
        alphasize = 1,
        show_legend = show_legend,
        coord_equal = TRUE,
        na_color = "black",
        theme = ggplot2::theme_bw(),
        show_labels = FALSE
      ) +
        ggplot2::ggtitle(paste0("Region ", r)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
      
      region_plots[[as.character(r)]] <- p
    }
    
    # Arrange plots into pages
    n_plots <- length(region_plots)
    
    for (i in seq(1, n_plots, by = regions_per_page)) {
      idx <- i:min(i + regions_per_page - 1, n_plots)
      plots_subset <- region_plots[idx]
      
      # Use patchwork to arrange
      combined_plot <- patchwork::wrap_plots(
        plots_subset,
        ncol = ncol_grid,
        nrow = nrow_grid
      )
      
      print(combined_plot)
    }
    
    dev.off()
    message(paste0("Saved: ", pdf_path))
    
  # ========================================
  # MODE: highlight - One PDF per celltype
  # ========================================
  } else {  # mode == "highlight"
    message("Mode: highlight - creating one PDF per celltype")
    
    # Loop through each celltype
    for (ct in celltypes) {
      message(paste0("Processing celltype: ", ct))
      
      # Create PDF for this celltype
      pdf_name <- gsub("[^a-zA-Z0-9_-]", "_", ct)
      pdf_path <- file.path(out_dir, paste0(pdf_name, "_regions.pdf"))
      pdf(pdf_path, width = 12, height = 10)
      
      # Optional: Create summary page showing all regions
      if (include_summary) {
        message("  Creating summary page...")
        
        summary_plot <- xyplot(
          cluster_col,
          x_column = x_column,
          y_column = y_column,
          cls = cluster_cls,
          clusters = NULL,
          metadata = seu@meta.data,
          ptsize = ptsize,
          plotfirst = ct,
          plotfirst_on_top = TRUE,
          alphasize = background_alpha,
          show_legend = show_legend,
          coord_equal = TRUE,
          na_color = "black",
          theme = ggplot2::theme_bw(),
          show_labels = FALSE
        )
        
        # Add black overlay for highlighted celltype
        highlight_data <- seu@meta.data[
          seu@meta.data[[cluster_col]] == ct & !is.na(seu@meta.data[[cluster_col]]), 
        ]
        
        if (nrow(highlight_data) > 0) {
          summary_plot <- summary_plot +
            ggplot2::geom_point(
              data = highlight_data,
              ggplot2::aes(x = !!rlang::sym(x_column), y = !!rlang::sym(y_column)),
              color = "black",
              size = ptsize,
              alpha = 1
            )
        }
        
        summary_plot <- summary_plot +
          ggplot2::ggtitle(paste0("Summary: ", ct, " (all regions)")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
        
        print(summary_plot)
      }
      
      # Create plots for each region
      region_plots <- list()
      
      for (r in regions) {
        # Subset to this region
        region_cells <- rownames(seu@meta.data[seu@meta.data[[region_col]] == r, ])
        
        if (length(region_cells) == 0) {
          # Create empty plot for missing region
          p <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = paste0("Region ", r, "\n(No cells)")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
          region_plots[[as.character(r)]] <- p
          next
        }
        
        obj_region <- subset(seu, cells = region_cells)
        region_metadata <- obj_region@meta.data
        
        # Check if coordinates exist
        if (!all(c(x_column, y_column) %in% colnames(region_metadata))) {
          p <- ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = paste0("Region ", r, "\n(Missing coordinates)")) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
          region_plots[[as.character(r)]] <- p
          next
        }
        
        # Count cells of this celltype in this region
        n_cells_ct <- sum(region_metadata[[cluster_col]] == ct, na.rm = TRUE)
        
        # Create base plot with faded colors
        p <- xyplot(
          cluster_col,
          x_column = x_column,
          y_column = y_column,
          cls = cluster_cls,
          clusters = NULL,
          metadata = region_metadata,
          ptsize = ptsize,
          plotfirst = ct,
          plotfirst_on_top = TRUE,
          alphasize = background_alpha,
          show_legend = show_legend,
          coord_equal = TRUE,
          na_color = "black",
          theme = ggplot2::theme_bw(),
          show_labels = FALSE
        )
        
        # Add black overlay for highlighted celltype
        highlight_data <- region_metadata[
          region_metadata[[cluster_col]] == ct & !is.na(region_metadata[[cluster_col]]), 
        ]
        
        if (nrow(highlight_data) > 0) {
          p <- p +
            ggplot2::geom_point(
              data = highlight_data,
              ggplot2::aes(x = !!rlang::sym(x_column), y = !!rlang::sym(y_column)),
              color = "black",
              size = ptsize,
              alpha = 1
            )
        }
        
        # Add title with cell count
        title_text <- if (n_cells_ct == 0) {
          paste0("Region ", r, "\n", ct, " (absent)")
        } else {
          paste0("Region ", r, "\n", ct, " (n=", n_cells_ct, ")")
        }
        
        p <- p +
          ggplot2::ggtitle(title_text) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
        
        region_plots[[as.character(r)]] <- p
      }
      
      # Arrange plots into pages
      n_plots <- length(region_plots)
      
      for (i in seq(1, n_plots, by = regions_per_page)) {
        idx <- i:min(i + regions_per_page - 1, n_plots)
        plots_subset <- region_plots[idx]
        
        # Use patchwork to arrange
        combined_plot <- patchwork::wrap_plots(
          plots_subset,
          ncol = ncol_grid,
          nrow = nrow_grid
        )
        
        print(combined_plot)
      }
      
      dev.off()
      message(paste0("  Saved: ", pdf_path))
    }
    
    message(paste0("Completed all ", length(celltypes), " celltype region plots"))
  }
  
  return(invisible(NULL))
}

plot_umap_celltype_highlights <- function(
  seu,
  cluster_col,
  cluster_cls,
  out_dir,
  reduction = "umap",
  celltypes = NULL,
  background_alpha = 0.2,
  highlight_alpha = 1,
  point_size = 0.3,
  raster = TRUE,
  shuffle = TRUE,
  seed = 123,
  show_labels = FALSE,
  label_size = 3.5,
  show_legend = TRUE,
  output_format = c("png", "pdf"),
  include_summary = FALSE,
  na_color = "grey80"
) {
  
  output_format <- match.arg(output_format)
  
  # Validate reduction exists
  if (!reduction %in% names(seu@reductions)) {
    stop("Reduction not found: ", reduction)
  }
  
  # Validate cluster_col exists
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Column not found in meta.data: ", cluster_col)
  }
  
  # Get all celltypes
  all_celltypes <- unique(seu@meta.data[[cluster_col]])
  all_celltypes <- all_celltypes[!is.na(all_celltypes)]
  
  if (is.null(celltypes)) {
    celltypes <- all_celltypes
  } else {
    # Validate requested celltypes exist
    missing <- setdiff(celltypes, all_celltypes)
    if (length(missing) > 0) {
      warning("Celltypes not found in data: ", paste(missing, collapse = ", "))
      celltypes <- intersect(celltypes, all_celltypes)
    }
  }
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # If PDF mode, open PDF device
  if (output_format == "pdf") {
    pdf_path <- file.path(out_dir, paste0("celltype_", reduction, "_highlights.pdf"))
    pdf(pdf_path, width = 10, height = 8)
    message("Creating multi-page PDF: ", pdf_path)
    
    # Optional summary page
    if (include_summary) {
      message("  Creating summary page...")
      summary_plot <- plot_embedding(
        seu = seu,
        reduction = reduction,
        group.by = cluster_col,
        label = TRUE,
        label_size = label_size,
        alpha = 1,
        point_size = point_size,
        raster = raster,
        shuffle = shuffle,
        seed = seed,
        palette = cluster_cls,
        na_color = na_color,
        legend = show_legend
      ) +
        ggplot2::ggtitle(paste0("Summary: All celltypes (", reduction, ")")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
      
      print(summary_plot)
    }
  }
  
  # Loop through each celltype
  for (ct in celltypes) {
    message(paste0("Processing celltype: ", ct))
    
    # Count cells of this celltype
    n_cells <- sum(seu@meta.data[[cluster_col]] == ct, na.rm = TRUE)
    
    # Create base plot with faded colors using plot_embedding
    p <- plot_embedding(
      seu = seu,
      reduction = reduction,
      group.by = cluster_col,
      label = show_labels,
      label_size = label_size,
      alpha = background_alpha,  # Faded background
      point_size = point_size,
      raster = raster,
      shuffle = shuffle,
      seed = seed,
      palette = cluster_cls,
      na_color = na_color,
      legend = show_legend
    )
    
    # Extract embedding coordinates for overlay
    emb <- as.data.frame(Seurat::Embeddings(seu, reduction))
    colnames(emb)[1:2] <- c("Dim1", "Dim2")
    emb$group <- seu@meta.data[[cluster_col]]
    
    # Get highlight data for this celltype
    highlight_data <- emb[emb$group == ct & !is.na(emb$group), ]
    
    # Add black overlay for highlighted celltype
    if (nrow(highlight_data) > 0) {
      geom_fun <- if (raster) ggrastr::geom_point_rast else ggplot2::geom_point
      
      p <- p +
        geom_fun(
          data = highlight_data,
          ggplot2::aes(x = Dim1, y = Dim2),
          color = "black",
          size = point_size,
          alpha = highlight_alpha,
          inherit.aes = FALSE
        )
    }
    
    # Add title
    title_text <- paste0("Highlighted: ", ct, " (n=", n_cells, ")")
    p <- p +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))
    
    # Save based on output format
    if (output_format == "png") {
      safe_filename <- gsub("[^a-zA-Z0-9_-]", "_", ct)
      output_file <- file.path(out_dir, paste0("highlight_", safe_filename, "_", reduction, ".png"))
      
      ggplot2::ggsave(
        filename = output_file,
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      message(paste0("  Saved: ", output_file))
    } else {
      # PDF mode - just print to open device
      print(p)
    }
  }
  
  # Close PDF if in PDF mode
  if (output_format == "pdf") {
    dev.off()
    message(paste0("Completed PDF with ", length(celltypes), " celltype pages"))
  } else {
    message(paste0("Completed all ", length(celltypes), " PNG files"))
  }
  
  return(invisible(NULL))
}

fp_layout <- function(IST_obj, cluster_pal = NULL) {
    if (is.null(cluster_pal)) {
        cols <- InSituType::colorCellTypes(freqs = table(IST_obj$clust), palette = "brewers")
    } else {
        cols <- cluster_pal
    }

    flightpath <- InSituType::flightpath_layout(logliks = IST_obj$logliks, profiles = IST_obj$profiles)
    par(mar = c(0,0,0,0))
    plot(flightpath$cellpos, pch = 16, cex = 0.2, col = cols[IST_obj$clust])
    text(flightpath$clustpos[, 1], flightpath$clustpos[, 2], rownames(flightpath$clustpos), cex = 0.7)
    par(mar = c(5, 4, 4, 2) + 0.1)
}


plot_all_cts <- function(seu, cluster_col, cluster_pal = NULL, out_dir = ".", save_pngs = FALSE, folder_name = "CT_PLOTS") {
    # Extract the coordinates from metadata
    xy <- as.matrix(seu@meta.data[, c("x_slide_mm", "y_slide_mm")])

    if (save_pngs) {
        # Create the output directory if it doesn't exist
        out_dir <- file.path(out_dir, folder_name)
        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }
    }

    for (ct in unique(seu@meta.data[[cluster_col]])) {
        if (save_pngs) {
            png(file.path(out_dir, paste0("celltype_", ct, "_spread.png")), 
            width = diff(range(xy[,1]))*.7, height = diff(range(xy[,2]))*.7, units = "in", 
            res = 400)  # res of 400 is pretty good; 600 is publication-quality
            par(mar = c(0,0,0,0))
        }
        plot(xy, pch = 16, col = scales::alpha(cluster_pal[seu@meta.data[[cluster_col]]], 0.3), cex = 0.1,
            xlab = "", ylab = "", xaxt = "n", yaxt = "n")
        points(xy[semisup$clust == ct, ], pch = 16, cex = 0.1, col = "black")
        legend("top", legend = ct)
        if (save_pngs) {
            dev.off()
        }
    }

    if (save_pngs) {
        par(mar = c(5, 4, 4, 2) + 0.1)
    }
}

plot_hm_pdf <- function(IST_obj, out_dir = ".", file_name = "Cluster_mean_expression_profiles.pdf") {
    pdf(file.path(out_dir, paste0(file_name)), width = 6, height = 20)
    
    heatmap(sweep(IST_obj$profiles, 1, pmax(apply(IST_obj$profiles, 1, max), .2), "/"), scale = "none",
        main = "Cluster mean expression profiles")

    dev.off()
}



## DEG specific functions
plot_volcano <- function(df, 
                         x_col, 
                         y_col, 
                         x_trans = NULL, 
                         y_trans = "-log10", 
                         x_cut = 0.5, 
                         p_cut = 0.05, 
                         genes_of_interest = NULL,
                         sig_color = "red", 
                         insig_color = "grey50", 
                         poi_color = "blue") {
  
  # 1. Prepare data and apply transformations
  plt_df <- df
  plt_df$X <- plt_df[[x_col]]
  plt_df$Y <- plt_df[[y_col]]
  
  if (!is.null(x_trans) && x_trans == "log2") plt_df$X <- log2(plt_df$X)
  if (!is.null(y_trans) && y_trans == "-log10") plt_df$Y_plot <- -log10(plt_df$Y) else plt_df$Y_plot <- plt_df$Y
  
  # 2. Define Significance categories
  # Note: logic uses the transformed X (Effect) but the RAW p-value threshold (p_cut)
  plt_df$status <- "Not Significant"
  plt_df$status[plt_df$Y < p_cut & abs(plt_df$X) > x_cut] <- "Significant"
  
  # 3. Handle specific genes of interest (POI)
  plt_df$is_poi <- FALSE
  if (!is.null(genes_of_interest)) {
    plt_df$is_poi <- rownames(plt_df) %in% genes_of_interest
  }
  
  # 4. Determine plotting order (so blue points are on top)
  plt_df <- plt_df[order(plt_df$status, plt_df$is_poi), ]

  # 5. Build Plot
  g <- ggplot(plt_df, aes(x = X, y = Y_plot)) +
    # Background points
    geom_point(aes(color = status), alpha = 0.6, size = 1.5) +
    # Highlight Genes of Interest
    geom_point(data = subset(plt_df, is_poi), color = poi_color, size = 2.5) +
    # Vertical lines for effect size
    geom_vline(xintercept = c(-x_cut, x_cut), linetype = "dashed", color = "darkgrey") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    # Horizontal line for p-value
    geom_hline(yintercept = ifelse(y_trans == "-log10", -log10(p_cut), p_cut), 
               linetype = "dashed", color = "darkgrey") +
    # Labels for significant OR points of interest
    geom_text_repel(
      data = subset(plt_df, status == "Significant" | is_poi),
      aes(label = rownames(subset(plt_df, status == "Significant" | is_poi))),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5,
      # Color labels blue if they are POIs, otherwise black
      color = ifelse(subset(plt_df, status == "Significant" | is_poi)$is_poi, poi_color, "black")
    ) +
    scale_color_manual(values = c("Significant" = sig_color, "Not Significant" = insig_color)) +
    labs(
      x = paste("Effect Size (", x_col, ")"),
      y = paste(y_trans, "(", y_col, ")"),
      color = "Status"
    ) +
    theme_classic() +
    theme(legend.position = "top")

  return(g)
}

get_gene_bin_plot_list <- function(seu_obj, 
                                   genes, 
                                   bin_col, 
                                   bins_to_include = NULL, 
                                   bin_order = NULL,
                                   bin_colors = NULL) {
  
  # 1. Verify genes
  genes <- intersect(genes, rownames(seu_obj))
  if(length(genes) == 0) stop("None of the provided genes found.")
  
  # 2. Extract and Prepare Data
  # We use data.table for efficiency
  meta <- as.data.table(seu_obj@meta.data)
  expr_data <- Seurat::GetAssayData(seu_obj, slot = "data")[genes, , drop = FALSE]
  
  plot_list <- list()
  
  # 3. Loop through genes and create individual plots
  for (goi in genes) {
    
    # Construct DT for this specific gene
    dt <- data.table(
      expr = as.numeric(expr_data[goi, ]),
      bin = meta[[bin_col]]
    )
    
    # Filter and order bins
    if (!is.null(bins_to_include)) dt <- dt[bin %in% bins_to_include]
    if (!is.null(bin_order)) dt[, bin := factor(bin, levels = bin_order)]
    
    # Calculate Mean and Standard Error
    summary_dt <- dt[, .(
      mean_expr = mean(expr, na.rm = TRUE),
      se = sd(expr, na.rm = TRUE) / sqrt(.N)
    ), by = bin]
    
    # 4. Create the Plot
    p <- ggplot(summary_dt, aes(x = bin, y = mean_expr, fill = bin)) +
      geom_col(color = "black", width = 0.7) +
      geom_errorbar(aes(ymin = mean_expr - se, ymax = mean_expr + se), 
                    width = 0.2, color = "black") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Adds 10% space at top
      theme_classic() +
      labs(
        title = paste("Gene:", goi),
        subtitle = paste("Cell type:", celltype_oi),
        x = "Spatial Zone",
        y = "Mean Normalized Expression (+/- SE)",
        fill = "Zone"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none" # Hide legend to maximize space on the PDF page
      )
    
    if (!is.null(bin_colors)) p <- p + scale_fill_manual(values = bin_colors)
    
    plot_list[[goi]] <- p
  }
  
  return(plot_list)
}

generate_spatial_dotplot <- function(seu_obj, 
                                     genes, 
                                     bin_col,
                                     scale = TRUE,
                                     bins_to_include = NULL, 
                                     bin_order = NULL) {
  
  # 1. Filter genes to ensure they exist
  genes <- intersect(genes, rownames(seu_obj))
  
  # 2. Handle subsetting and ordering
  plot_obj <- seu_obj
  
  # Filter bins if requested
  if (!is.null(bins_to_include)) {
    plot_obj <- subset(plot_obj, cells = colnames(plot_obj)[plot_obj[[bin_col, drop=TRUE]] %in% bins_to_include])
  }
  
  # Enforce factor levels for the bin column to control plot order
  if (!is.null(bin_order)) {
    plot_obj@meta.data[[bin_col]] <- factor(plot_obj@meta.data[[bin_col]], levels = bin_order)
  }

  # 3. Generate DotPlot using group.by
  p <- DotPlot(
    object = plot_obj, 
    features = genes, 
    group.by = bin_col,      # Use group.by instead of changing Idents
    cols = c("lightgrey", "firebrick"), 
    scale = scale, 
    col.min = -2.5, 
    col.max = 2.5
  ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Spatial Gene Expression",
      x = "Genes",
      y = "Spatial Zone",
      size = "Percent Expressed",
      color = "Average Expression (Z-score)"
    )
  
  return(p)
}




# FastReSeg

estimate_MeanProfile_sparse <- function(counts, clust, bg, s) {
  stopifnot(
    inherits(counts, "dgCMatrix"),
    length(clust) == nrow(counts),
    length(bg) == nrow(counts),
    length(s) == nrow(counts)
  )

  clust <- as.factor(clust)
  clusters <- levels(clust)
  n_genes <- ncol(counts)

  # Preallocate result: genes x clusters
  means <- matrix(0, nrow = n_genes, ncol = length(clusters))
  colnames(means) <- clusters
  rownames(means) <- colnames(counts)

  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    idx <- which(clust == cl)
    if (length(idx) == 0) next

    # Sparse-safe per-cluster mean
    cluster_counts <- counts[idx, , drop = FALSE]        # subset cells in cluster
    scaled_counts <- sweep(pmax(cluster_counts - bg[idx], 0), 1, s[idx], "/")
    gene_means <- Matrix::colMeans(scaled_counts)

    means[, i] <- gene_means
  }

  means
}


# EDA
plot_module_details <- function(res, xy, output_file = "insitucor_module_details.pdf") {
  
  require(ComplexHeatmap)
  require(viridis)
  require(grid)
  require(circlize)
  
  if (file.exists(output_file)) file.remove(output_file)
  
  col_fun <- colorRamp2(seq(min(res$celltypeinvolvement), max(res$celltypeinvolvement), length = 100),
                        colorRampPalette(c("white", "darkblue"))(100))
  
  module_names <- colnames(res$scores_env)
  n_modules <- length(module_names)
  
  pdf(output_file, width = 11, height = 8.5, onefile = TRUE)
  
  draw(Heatmap(res$celltypeinvolvement, 
               col = col_fun,
               name = "Score",
               column_title = "Contribution of cell types to module scores",
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10)))
  
  for (i in seq_along(module_names)) {
    module_name <- module_names[i]
    cat("Processing module", i, "of", n_modules, ":", module_name, "\n")
    
    grid.newpage()
    grid.text(module_name, gp = gpar(fontsize = 24, fontface = "bold"))
    
    if (!is.null(res$attributionmats) && i <= length(res$attributionmats)) {
      mat <- res$attributionmats[[i]]
      col_fun_attr <- colorRamp2(seq(min(mat), max(mat), length = 100),
                                 colorRampPalette(c("white", "darkblue"))(100))
      
      draw(Heatmap(mat,
                   col = col_fun_attr,
                   name = "Score",
                   column_title = paste0(module_name, "\nContribution of cell types to module genes"),
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10)))
    } else {
      grid.newpage()
      grid.text("Attribution matrix not available", gp = gpar(fontsize = 14, fontface = "italic"))
    }
    
    tempscore <- res$scores_sc[, i]
    colors <- viridis_pal(option = "B")(101)[1 + pmin(round(100 * (tempscore / quantile(tempscore, 0.995))), 100)]
    plot(xy, asp = 1, pch = 16, cex = 0.5, main = paste0(module_name, "\nSingle cell score"),
         xlab = "x mm", ylab = "y mm", col = colors, cex.main = 1.2)
    
    tempscore <- res$scores_env[, i]
    colors <- viridis_pal(option = "B")(101)[1 + pmin(round(100 * (tempscore / quantile(tempscore, 0.995))), 100)]
    plot(xy, asp = 1, pch = 16, cex = 0.5, main = paste0(module_name, "\nNeighborhood score"),
         xlab = "x mm", ylab = "y mm", col = colors, cex.main = 1.2)
  }
  
  dev.off()
  
  cat("\n✓ Created:", output_file, "\n")
  cat("✓ Pages:", 1 + (n_modules * 4), "\n")
  
  invisible(list(n_modules = n_modules, module_names = module_names, output_file = output_file))
}

summarize_insitucor <- function(res, corthresh = 0.25, output_file = "insitucor_summary.pdf") {
  
  tryCatch({
    # Start PDF device
    pdf(output_file, width = 11, height = 8.5, onefile = TRUE)
    
    module_dims <- dim(res$modules)
    n_unique_modules <- length(unique(res$modules$module))
    
    stats_text <- textGrob(
      paste0(
        "InSituCor Results Summary\n\n\n",
        "Module Dimensions:\n",
        module_dims[1], " rows × ", module_dims[2], " columns\n\n\n",
        "Number of Unique Modules:\n",
        n_unique_modules
      ),
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    grid.draw(stats_text)

    
    genes_per_module <- table(res$modules$module)
    genes_df <- data.frame(
      module = names(genes_per_module),
      n_genes = as.numeric(genes_per_module)
    )
    
    p_barplot <- ggplot(genes_df, aes(x = reorder(module, n_genes), y = n_genes)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Number of Genes per Module",
           x = "Module",
           y = "Number of Genes") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.y = element_text(size = 10)
      )

    print(p_barplot)

    
    # Determine how many rows to show (max 50 for readability)
    n_rows_to_show <- min(25, nrow(res$modules))
    table_data <- head(res$modules, n_rows_to_show)
    
    # Create table grob
    table_grob <- tableGrob(
      table_data,
      rows = NULL,
      theme = ttheme_default(
        core = list(fg_params = list(cex = 0.6)),
        colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
      )
    )
    
    # Add title and footer
    title_grob <- textGrob(
      "Module Table", 
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    footer_grob <- textGrob(
      paste0("Showing first ", n_rows_to_show, " of ", nrow(res$modules), " total rows"),
      gp = gpar(fontsize = 10, fontface = "italic")
    )
    
    grid.arrange(
      title_grob,
      table_grob,
      footer_grob,
      heights = c(0.1, 0.85, 0.05),
      ncol = 1
    )

    plotCorrelationNetwork(
      res$condcor, 
      modules = res$modules, 
      show_gene_names = FALSE, 
      corthresh = corthresh
    )
    grid.text(
      paste0("Correlation Network (corthresh = ", corthresh, ")"),
      x = 0.5, y = 0.97, 
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    
    plotCorrelationNetwork(
      res$condcor, 
      modules = res$modules, 
      show_gene_names = TRUE, 
      corthresh = corthresh
    )
    grid.text(
      paste0("Correlation Network with Gene Names (corthresh = ", corthresh, ")"),
      x = 0.5, y = 0.97,
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    
    dev.off()
    
    # Verify file was created
    if (file.exists(output_file) && file.size(output_file) > 0) {
      cat("✓ Summary report created successfully:", output_file, "\n")
      cat("✓ File size:", round(file.size(output_file) / 1024, 2), "KB\n")
    } else {
      warning("PDF file may be corrupted or empty")
    }
    
    # Return summary statistics
    invisible(list(
      dimensions = module_dims,
      n_modules = n_unique_modules,
      genes_per_module = genes_df,
      output_file = output_file
    ))
    
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    stop("Error creating PDF: ", e$message)
  })
}

analyze_module_correlations <- function(res, output_file = "module_correlation_stats.pdf") {
  
  require(ggplot2)
  require(gridExtra)
  
  cormat <- res$condcor
  module_list <- split(res$modules$gene, res$modules$module)
  
  # Calculate mean correlation for each module
  calculate_mean_cor <- function(genes, cormat) {
    genes <- genes[genes %in% rownames(cormat)]
    if (length(genes) < 2) return(NA)
    
    submat <- cormat[genes, genes]
    n <- length(genes)
    meancor <- (sum(submat) - n) / (n^2 - n)
    return(meancor)
  }
  
  module_stats <- data.frame(
    module = names(module_list),
    n_genes = sapply(module_list, length),
    mean_cor = sapply(module_list, calculate_mean_cor, cormat = cormat),
    stringsAsFactors = FALSE
  )
  
  module_stats <- module_stats[order(module_stats$n_genes, decreasing = TRUE), ]
  
  pdf(output_file, width = 11, height = 8.5)
  
  # Plot 1: Histogram of mean correlations
  p1 <- ggplot(module_stats, aes(x = mean_cor)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    geom_vline(xintercept = median(module_stats$mean_cor, na.rm = TRUE), 
               linetype = "dashed", color = "red", size = 1) +
    labs(title = "Distribution of Mean Conditional Correlations per Module",
         subtitle = paste0("Median: ", round(median(module_stats$mean_cor, na.rm = TRUE), 3)),
         x = "Mean Conditional Correlation",
         y = "Number of Modules") +
    theme_minimal(base_size = 14)
  
  # Plot 2: Scatter plot - module size vs mean correlation
  p2 <- ggplot(module_stats, aes(x = n_genes, y = mean_cor)) +
    geom_point(size = 3, alpha = 0.6, color = "darkblue") +
    geom_smooth(method = "loess", se = TRUE, color = "red") +
    labs(title = "Module Size vs Mean Correlation",
         x = "Number of Genes in Module",
         y = "Mean Conditional Correlation") +
    theme_minimal(base_size = 14)
  
  # Plot 3: Bar plot of correlations ordered by module
  module_stats$module <- factor(module_stats$module, levels = module_stats$module)
  p3 <- ggplot(module_stats, aes(x = module, y = mean_cor, fill = n_genes)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "# Genes") +
    labs(title = "Mean Correlation by Module (ordered by size)",
         x = "Module",
         y = "Mean Conditional Correlation") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  print(p1)
  print(p2)
  print(p3)
  
  # Summary statistics page
  grid.newpage()
  summary_text <- paste0(
    "Module Correlation Summary Statistics\n\n",
    "Number of Modules: ", nrow(module_stats), "\n\n",
    "Mean Correlation Statistics:\n",
    "  Min:    ", round(min(module_stats$mean_cor, na.rm = TRUE), 3), "\n",
    "  Q1:     ", round(quantile(module_stats$mean_cor, 0.25, na.rm = TRUE), 3), "\n",
    "  Median: ", round(median(module_stats$mean_cor, na.rm = TRUE), 3), "\n",
    "  Mean:   ", round(mean(module_stats$mean_cor, na.rm = TRUE), 3), "\n",
    "  Q3:     ", round(quantile(module_stats$mean_cor, 0.75, na.rm = TRUE), 3), "\n",
    "  Max:    ", round(max(module_stats$mean_cor, na.rm = TRUE), 3), "\n\n",
    "Module Size Statistics:\n",
    "  Min:    ", min(module_stats$n_genes), "\n",
    "  Median: ", median(module_stats$n_genes), "\n",
    "  Max:    ", max(module_stats$n_genes)
  )
  
  grid.text(summary_text, x = 0.5, y = 0.5, 
            gp = gpar(fontsize = 14, fontfamily = "mono"))
  
  dev.off()
  
  cat("\n✓ Created:", output_file, "\n")
  
  invisible(module_stats)
}



prepare_tma_data <- function(counts, norm, metadata, xy, cores = NULL, region_col = "region", min_cells_per_core = 3000, study_id_col = "study_id") {
  
  # If no cores specified, use all unique cores
  if (is.null(cores)) {
    cores <- unique(metadata[[region_col]])
  }
  
  # Initialize output list
  prepared_data <- list()
  
  cat("Preparing data for", length(cores), "cores...\n")
  
  for (core_id in cores) {
    # Subset to current core
    core_cells <- metadata[[region_col]] == core_id
    
    # Skip if too few cells
    if (sum(core_cells) < min_cells_per_core) {
      cat("  Skipping", core_id, "- too few cells (", sum(core_cells), ")\n")
      next
    }
    
    # Extract core-specific data
    core_data <- list(
      core_id = core_id,
      counts = counts[core_cells, ],
      norm = norm[core_cells, ],
      metadata = metadata[core_cells, ],
      xy = xy[core_cells, ],
      n_cells = sum(core_cells),
      patient_id = unique(metadata[[study_id_col]][core_cells])[1]  # Assuming one patient per core
    )
    
    # Verify data integrity
    if (core_data$n_cells != nrow(core_data$counts)) {
      stop("Dimension mismatch for core ", core_id)
    }
    
    prepared_data[[core_id]] <- core_data
    cat("  ", core_id, ":", core_data$n_cells, "cells from patient", core_data$patient_id, "\n")
  }
  
  return(prepared_data)
}


# ==============================================================================
# Generalized InsituCor pipeline functions (per-unit: core OR patient)
# ==============================================================================

#' Prepare data for InsituCor analysis, grouped by an arbitrary unit
#'
#' @param unit_col Column to group by ("region" for per-core, "study_id" for per-patient)
#' @param tissue_col Column to use as `tissue` in InsituCor (prevents cross-tissue neighbors).
#'   For per-patient analysis set to "region" so cores within a patient don't share neighbors.
#'   For per-core analysis set to NULL.
prepare_insitucor_data <- function(counts, norm, metadata, xy,
                                    units = NULL,
                                    unit_col = "region",
                                    tissue_col = NULL,
                                    min_cells_per_unit = 3000) {

  if (is.null(units)) {
    units <- unique(metadata[[unit_col]])
  }

  prepared_data <- list()
  cat("Preparing data for", length(units), "units (grouped by", unit_col, ")...\n")

  for (uid in units) {
    unit_cells <- metadata[[unit_col]] == uid

    if (sum(unit_cells) < min_cells_per_unit) {
      cat("  Skipping", uid, "- too few cells (", sum(unit_cells), ")\n")
      next
    }

    unit_data <- list(
      unit_id    = uid,
      counts     = counts[unit_cells, ],
      norm       = norm[unit_cells, ],
      metadata   = metadata[unit_cells, ],
      xy         = xy[unit_cells, ],
      n_cells    = sum(unit_cells),
      tissue_col = tissue_col
    )

    if (unit_data$n_cells != nrow(unit_data$counts)) {
      stop("Dimension mismatch for unit ", uid)
    }

    prepared_data[[as.character(uid)]] <- unit_data
    cat("  ", uid, ":", unit_data$n_cells, "cells\n")
  }

  return(prepared_data)
}


#' Run InsituCor on a single analysis unit
#'
#' Automatically derives the `tissue` argument from `unit_data$tissue_col` so that
#' per-patient runs correctly prevent cross-core neighbor links.
run_insitucor_single_unit <- function(unit_data,
                                      conditionon_cols = c("fov", "nCount_RNA", "negmean", "final_annotation"),
                                      k = 50,
                                      radius = NULL,
                                      annotation_col = "final_annotation",
                                      min_module_size = 3,
                                      min_module_cor = 0.25,
                                      max_cells = 1e5,
                                      roundcortozero = 0.1,
                                      attribution_subset_size = 1000,
                                      verbose = TRUE,
                                      output_dir = NULL,
                                      study_name = NULL) {

  cat("\n=== Running InsituCor on unit:", unit_data$unit_id, "===\n")
  cat("Cells:", unit_data$n_cells, "\n")

  if (unit_data$n_cells > max_cells) {
    cat("Subsampling from", unit_data$n_cells, "to", max_cells, "cells will be performed in the InsituCor run\n")
  }

  counts_use   <- unit_data$counts
  norm_use     <- unit_data$norm
  metadata_use <- unit_data$metadata
  xy_use       <- unit_data$xy

  # Prepare conditioning variables
  conditionon <- as.data.frame(metadata_use[, conditionon_cols, drop = FALSE])

  if (any(is.na(conditionon))) {
    cat("WARNING: Missing values in conditioning variables. Removing affected cells.\n")
    cc <- complete.cases(conditionon)
    counts_use   <- counts_use[cc, ]
    norm_use     <- norm_use[cc, ]
    metadata_use <- metadata_use[cc, ]
    xy_use       <- xy_use[cc, ]
    conditionon  <- conditionon[cc, , drop = FALSE]
  }

  # Derive tissue from stored tissue_col
  if (!is.null(unit_data$tissue_col)) {
    tissue <- metadata_use[[unit_data$tissue_col]]
    cat("Using tissue column:", unit_data$tissue_col,
        "(", length(unique(tissue)), "unique values)\n")
  } else {
    tissue <- NULL
  }

  res <- tryCatch({
    InSituCor::insitucor(
      counts    = norm_use,
      conditionon = conditionon,
      celltype  = metadata_use[[annotation_col]],
      neighbors = NULL,
      xy        = xy_use,
      k         = k,
      radius    = radius,
      tissue    = tissue,
      min_module_size    = min_module_size,
      max_module_size    = 25,
      min_module_cor     = min_module_cor,
      gene_weighting_rule = "inverse_sqrt",
      roundcortozero     = roundcortozero,
      max_cells          = max_cells,
      attribution_subset_size = attribution_subset_size,
      verbose   = verbose
    )
  }, error = function(e) {
    cat("ERROR in InsituCor for unit", unit_data$unit_id, ":", e$message, "\n")
    return(NULL)
  })

  if (!is.null(output_dir) && !is.null(res)) {
    dir.create(file.path(output_dir, "InsituCor_Results_perUnit"),
               showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir, "InsituCor_Results_perUnit",
                             paste0("InsituCor_Results_", study_name,
                                    "_unit", unit_data$unit_id,
                                    "_minModuleSize", min_module_size,
                                    "_minModuleCor", min_module_cor,
                                    "_k", k, ".rds"))
    saveRDS(res, file = output_file)
    cat("Saved results to:", output_file, "\n")
  }

  return(res)
}


#' Run InsituCor sequentially across multiple units
run_insitucor_multiple_units <- function(prepared_data,
                                         fun = run_insitucor_single_unit,
                                         conditionon_cols = c("fov", "nCount_RNA", "negmean", "final_annotation"),
                                         k = 50,
                                         annotation_col = "final_annotation",
                                         min_module_size = 3,
                                         min_module_cor = 0.25,
                                         output_dir = NULL,
                                         study_name = "",
                                         verbose = FALSE) {

  cat("\nRunning InsituCor sequentially on", length(prepared_data), "units...\n")
  results <- lapply(prepared_data, fun,
      conditionon_cols = conditionon_cols,
      k = k,
      annotation_col = annotation_col,
      min_module_size = min_module_size,
      min_module_cor = min_module_cor,
      output_dir = output_dir,
      study_name = study_name,
      verbose = verbose
  )

  results <- results[!sapply(results, is.null)]
  cat("\nSuccessfully analyzed", length(results), "units\n")
  return(results)
}


#' Extract unit ID from result filename
extract_unit_id_from_filename <- function(filename) {
  pattern <- ".*_unit(.+?)_minModuleSize.*"
  uid <- sub(pattern, "\\1", filename)
  if (uid == filename) {
    pattern2 <- ".*unit(.+?)\\.rds$"
    uid <- sub(pattern2, "\\1", filename)
  }
  uid
}


#' Load per-unit InsituCor results from disk
load_insitucor_unit_results <- function(results_dir, unit_ids = NULL) {

  files <- list.files(path = results_dir, pattern = "\\.rds$",
                      full.names = TRUE, recursive = TRUE)

  if (length(files) == 0) {
    stop("No RDS files found in results_dir: ", results_dir)
  }

  if (is.null(unit_ids)) {
    unit_ids <- unique(sapply(files, extract_unit_id_from_filename))
  }

  results <- vector("list", length(files))
  names(results) <- unit_ids

  for (i in seq_along(files)) {
    results[[i]] <- readRDS(files[i])
  }

  results
}


run_insitucor_single_core <- function(core_data, 
                                      conditionon_cols = c("fov", "nCount_RNA", "negmean", "final_annotation"),
                                      k = 50,
                                      radius = NULL,
                                      tissue = NULL,
                                      annotation_col = "final_annotation",
                                      min_module_size = 3,
                                      min_module_cor = 0.25,
                                      max_cells = 1e5,
                                      roundcortozero = 0.1,
                                      attribution_subset_size = 1000,
                                      verbose = TRUE,
                                      output_dir = NULL,
                                      study_name = NULL) {
  
  cat("\n=== Running InsituCor on core:", core_data$core_id, "===\n")
  cat("Cells:", core_data$n_cells, "\n")
  
  # Check if subsampling needed
  if (core_data$n_cells > max_cells) {
    cat("Subsampling from", core_data$n_cells, "to", max_cells, "cells will be performed in the InsituCor run\n")
  }
  counts_use <- core_data$counts
  norm_use <- core_data$norm
  metadata_use <- core_data$metadata
  xy_use <- core_data$xy
  
  # Prepare conditioning variables
  conditionon <- as.data.frame(metadata_use[, conditionon_cols, drop = FALSE])
  
  # Check for missing values in conditioning variables
  if (any(is.na(conditionon))) {
    cat("WARNING: Missing values in conditioning variables. Removing affected cells.\n")
    complete_cases <- complete.cases(conditionon)
    counts_use <- counts_use[complete_cases,]
    norm_use <- norm_use[complete_cases,]
    metadata_use <- metadata_use[complete_cases, ]
    xy_use <- xy_use[complete_cases, ]
    conditionon <- conditionon[complete_cases, , drop = FALSE]
  }
  
  # Run InsituCor
  res <- tryCatch({
    InSituCor::insitucor(
      # Fundamental input data
      counts = norm_use, 
      conditionon = conditionon, 
      celltype = metadata_use[[annotation_col]],
      
      # Args for neighbor definition
      neighbors = NULL, 
      xy = xy_use, 
      k = k, 
      radius = radius, 
      tissue = tissue, 
      
      # Args for module definition
      min_module_size = min_module_size, 
      max_module_size = 25, 
      min_module_cor = min_module_cor,
      gene_weighting_rule = "inverse_sqrt",   
      
      # Args for controlling memory and compute
      roundcortozero = roundcortozero, 
      max_cells = max_cells,
      
      # Args for cell type attribution scoring
      attribution_subset_size = attribution_subset_size,
      
      verbose = verbose
    )
  }, error = function(e) {
    cat("ERROR in InsituCor for core", core_data$core_id, ":", e$message, "\n")
    return(NULL)
  })
  
  # Save results if output directory specified
  if (!is.null(output_dir) && !is.null(res)) {
    dir.create(file.path(output_dir, "InsituCor_Results_perCore"), showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir, "InsituCor_Results_perCore" ,paste0("InsituCor_Results_", study_name, "_core", core_data$core_id,"_minModuleSize", min_module_size, "_minModuleCor", min_module_cor, "_k", k,".rds"))
    saveRDS(res, file = output_file)
    cat("Saved results to:", output_file, "\n")
  }
  
  return(res)
}

run_insitucor_multiple_cores <- function(prepared_data, 
    fun, 
    conditionon_cols = c("fov", "nCount_RNA", "negmean", "final_annotation"), 
    k = 50, 
    annotation_col = "final_annotation", 
    min_module_size = 3,
    min_module_cor = 0.25,
    output_dir = NULL,
    study_name = "",
    verbose = FALSE
) {
    cat("\nRunning InsituCor sequentially...\n")
    results <- lapply(prepared_data, fun, 
        conditionon_cols = conditionon_cols, 
        k = k, 
        annotation_col = annotation_col, 
        min_module_size = min_module_size,
        min_module_cor = min_module_cor,
        output_dir = output_dir,
        study_name = study_name,
        verbose = verbose
    )
  
  # Remove NULL results (failed cores)
  results <- results[!sapply(results, is.null)]
  
  cat("\nSuccessfully analyzed", length(results), "cores\n")
  return(results)
}

inspect_results_directory <- function(results_dir) {
  
  # Find all RDS files
  rds_files <- list.files(
    path = results_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(rds_files) == 0) {
    cat("No RDS files found in", results_dir, "\n")
    return(NULL)
  }
  
  # Extract information from each file
  file_info <- data.frame(
    file = basename(rds_files),
    size_mb = file.size(rds_files) / 1024^2,
    modified = file.mtime(rds_files),
    stringsAsFactors = FALSE
  )
  
  cat("Found", nrow(file_info), "RDS files\n")
  cat("\nFiles found:", paste(file_info$file, collapse = ", "), "\n")
  cat("Total size:", round(sum(file_info$size_mb), 2), "MB\n")
  
  return(file_info)
}

extract_core_id_from_filename <- function(filename) {
  pattern <- ".*_core(.+?)_minModuleSize.*"
  core_id <- sub(pattern, "\\1", filename)

  if (core_id == filename) {
    pattern2 <- ".*core(.+?)\\.rds$"
    core_id <- sub(pattern2, "\\1", filename)
  }

  core_id
}


load_insitucor_results <- function(results_dir, core_ids = NULL) {

  files <- list.files(
    path = results_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files) == 0) {
    stop("No RDS files found in results_dir")
  }

  if (is.null(core_ids)) {
    core_ids <- unique(sapply(files, extract_core_id_from_filename))
  }

  results <- vector("list", length(files))
  names(results) <- core_ids

  for (i in seq_along(files)) {
    results[[i]] <- readRDS(files[i])
  }

  results
}

extract_condcor_array <- function(results) {
  
  # Get dimensions from first result
  first_result <- results[[1]]
  genes <- rownames(first_result$condcor)
  n_genes <- length(genes)
  n_cores <- length(results)

  cat("Extracting condcor matrices for", n_cores, "units and", n_genes, "genes...\n")
  
  # Initialize array
  condcors <- array(NA, 
                    dim = c(n_cores, n_genes, n_genes),
                    dimnames = list(names(results), genes, genes))
  
  # Fill array
  for (i in seq_along(results)) {
    core_name <- names(results)[i]
    condcors[core_name, , ] <- as.matrix(results[[i]]$condcor[genes, genes])
  }
  
  return(condcors)
}

# build_consensus_network <- function(condcors, 
#                                     method = "lowq", 
#                                     threshold = 0.2,
#                                     min_fraction = NULL) {
  
#   n_cores <- dim(condcors)[1]
#   n_genes <- dim(condcors)[2]
  
#   # Determine quantile from min_fraction
#   if (!is.null(min_fraction)) {
#     if (min_fraction < 0 || min_fraction > 1) {
#       stop("min_fraction must be between 0 and 1")
#     }
#     # Convert fraction to quantile
#     # min_fraction = 0.5 means "present in at least 50% of cores"
#     quantile <- 1 - min_fraction
#     cat("Using min_fraction =", min_fraction, 
#         "→ quantile =", round(quantile, 3),
#         "(require edge in ≥", ceiling(n_cores * min_fraction), "cores)\n")
#   } else {
#     min_fraction <- 1
#     quantile <- 1 - min_fraction
#     cat("No min_fraction provided, all cores required.\n",
#         "Using min_fraction =", min_fraction, 
#         "→ quantile =", round(quantile, 3),
#         "(require edge in ≥", ceiling(n_cores * min_fraction), "cores)\n")
#   }
  
#   cat("Building consensus network using method:", method, "\n")
#   cat("Threshold:", threshold, ", Quantile:", round(quantile, 3), "\n")

#   # # Check symmetry of individual cores
#   # for (i in seq_len(dim(condcors)[1])) {
#   #   cat("Symmetry Check (core", i, "):", 
#   #       isSymmetric(condcors[i, , ]), "\n")
#   # }

#   # # Get indices of upper triangle (excluding diagonal)
#   # tri_idx <- which(upper.tri(matrix(0, n_genes, n_genes)), arr.ind = TRUE)
#   # n_pairs <- nrow(tri_idx)
#   # cat("Number of pairs:", n_pairs, "\n")
#   # # Extract upper triangle for all cores: [n_cores, n_pairs]
#   # condcors_upper <- matrix(NA, nrow = n_cores, ncol = nrow(tri_idx))
#   # for (i in 1:n_cores) {
#   #   # Extract upper triangle values for core i
#   #   condcors_upper[i, ] <- condcors[i, tri_idx[,1], tri_idx[,2]]
#   # }

#   # # Compute summary statistic per gene pair
#   # if (method == "lowq") {
#   #   summary_vec <- apply(condcors_upper, 2, quantile, probs = quantile, na.rm = TRUE)
#   # } else if (method == "maxq") {
#   #   summary_vec <- apply(condcors_upper, 2, max, na.rm = TRUE)
#   # } else if (method == "mean") {
#   #   summary_vec <- colMeans(condcors_upper, na.rm = TRUE)
#   # } else {
#   #   stop("Unknown method")
#   # }

#   # # Fill full symmetric matrix
#   # summary_cor <- matrix(0, n_genes, n_genes)
#   # summary_cor[tri_idx] <- summary_vec
#   # summary_cor[lower.tri(summary_cor)] <- t(summary_cor)[lower.tri(summary_cor)]
#   # diag(summary_cor) <- 0

#   # # After building consensus network
#   # cat("Symmetry Check (final consensus):", 
#   #     isSymmetric(summary_cor), "\n")

#   # # Restore dimnames
#   # rownames(summary_cor) <- dimnames(condcors)[[2]]
#   # colnames(summary_cor) <- dimnames(condcors)[[3]]


#   # Flatten array: [cores, genes*genes]
#   condcors_mat <- matrix(condcors, nrow = n_cores, ncol = n_genes^2)

#   # Compute summary statistic
#   if (method == "lowq") {
#     summary_vec <- apply(condcors_mat, 2, quantile, probs = quantile, na.rm = TRUE)
#   } else if (method == "maxq") {
#     summary_vec <- apply(condcors_mat, 2, max, na.rm = TRUE)
#   } else if (method == "mean") {
#     summary_vec <- colMeans(condcors_mat, na.rm = TRUE)
#   } else {
#     stop("Unknown method")
#   }

#   # Reshape back to [genes, genes]
#   summary_cor <- matrix(summary_vec, nrow = n_genes, ncol = n_genes)
#   rownames(summary_cor) <- dimnames(condcors)[[2]]
#   colnames(summary_cor) <- dimnames(condcors)[[3]]
  
#   # # Compute summary statistics
#   # if (method == "lowq") {
#   #   summary_cor <- apply(condcors, 2:3, quantile, probs = quantile, na.rm = TRUE)
#   # } else if (method == "maxq") {
#   #   summary_cor <- apply(condcors, 2:3, max, na.rm = TRUE)
#   # } else if (method == "mean") {
#   #   summary_cor <- apply(condcors, 2:3, mean, na.rm = TRUE)
#   # }
  
#   diag(summary_cor) <- 0
#   consensus <- (summary_cor > threshold) * 1
#   keep_genes <- (rowSums(consensus) > 0) & (colSums(consensus) > 0)
#   consensus <- consensus[keep_genes, keep_genes]
  
#   cat("Consensus network has", nrow(consensus), "genes and", 
#       sum(consensus)/2, "edges\n")
  
#   return(list(
#     adjacency = consensus,
#     correlations = summary_cor[keep_genes, keep_genes]
#   ))
# }

name_module <- function(mat) {
  meancors <- Matrix::colMeans(mat)
  n <- ncol(mat)
  top3 <- colnames(mat)[order(meancors, decreasing = TRUE)[1:min(n, 3)]]
  if (n <= 3) {
    name <- paste0(c(top3, n), collapse = "_")
  } else {
    name <- paste0(c(top3[1:2], n), collapse = "_")
  }
  # replace problematic special characters
  name <- make.names(name)
  return(name)
}

build_consensus_network <- function(condcors, 
                                    method = "lowq", 
                                    threshold = 0.2,
                                    min_fraction = NULL) {
  
  n_cores <- dim(condcors)[1]
  n_genes <- dim(condcors)[2]
  
  # Determine quantile from min_fraction
  if (!is.null(min_fraction)) {
    if (min_fraction < 0 || min_fraction > 1) {
      stop("min_fraction must be between 0 and 1")
    }
    quantile <- 1 - min_fraction
    cat("Using min_fraction =", min_fraction, 
        "→ quantile =", round(quantile, 3),
        "(require edge in ≥", ceiling(n_cores * min_fraction), "cores)\n")
  } else {
    min_fraction <- 1
    quantile <- 0
    cat("No min_fraction provided, all cores required.\n")
  }
  
  cat("Building consensus network using method:", method, "\n")
  cat("Threshold:", threshold, "\n")
  
  # Initialize result matrix
  summary_cor <- matrix(0, n_genes, n_genes)
  
  # Determine target rank for lowq method
  if (method == "lowq") {
    target_rank <- ceiling(n_cores * quantile)
    if (target_rank < 1) target_rank <- 1
  }
  
  # Process row by row (memory efficient, still leverages symmetry)
  cat("Processing", n_genes, "rows...\n")
  
  for (row in 1:n_genes) {
    if (row %% 500 == 0) cat("  Row", row, "of", n_genes, "\n")
    
    # Only process upper triangle: columns from row+1 to n_genes
    if (row < n_genes) {
      cols <- (row + 1):n_genes
      n_cols <- length(cols)
      
      # Extract values for this row across all cores: [n_cores, n_cols]
      row_vals <- matrix(NA, nrow = n_cores, ncol = n_cols)
      for (i in 1:n_cores) {
        row_vals[i, ] <- condcors[i, row, cols]
      }
      
      # Compute summary statistic
      if (method == "mean") {
        row_summary <- colMeans(row_vals, na.rm = TRUE)
        
      } else if (method == "maxq") {
        row_summary <- apply(row_vals, 2, max, na.rm = TRUE)
        
      } else if (method == "lowq") {
        row_summary <- apply(row_vals, 2, quantile, probs = quantile, na.rm = TRUE)
      }
      
      # Store in upper triangle
      summary_cor[row, cols] <- row_summary
    }
  }
  
  # Mirror to lower triangle
  summary_cor <- summary_cor + t(summary_cor)
  diag(summary_cor) <- 0
  
  rownames(summary_cor) <- dimnames(condcors)[[2]]
  colnames(summary_cor) <- dimnames(condcors)[[3]]
  
  # Create consensus network
  cat("Creating consensus network...\n")
  consensus <- (summary_cor > threshold) * 1
  keep_genes <- (rowSums(consensus) > 0) & (colSums(consensus) > 0)
  consensus <- consensus[keep_genes, keep_genes]
  
  cat("Consensus network has", nrow(consensus), "genes and", 
      sum(consensus)/2, "edges\n")
  
  return(list(
    adjacency = consensus,
    correlations = summary_cor[keep_genes, keep_genes]
  ))
}


define_modules_consensus<- function(consensus, nh_expression=NULL ,resolution = 0.02, min_module_size = 3, max_module_size = 20, weighting = "inverse_sqrt") {
  require(igraph)
   adj_matrix <- as.matrix(consensus$adjacency)
   if (!isSymmetric(adj_matrix)) {
     stop("Consensus adjacency matrix is not symmetric.")
   }

  message("Converting adjacency matrix to igraph object...")
  gr <- igraph::graph_from_adjacency_matrix(adj_matrix, 
                                weighted = TRUE,
                                mode = "undirected")

  message("Performing leiden clustering...")
  leid <- igraph::cluster_leiden(gr, weights = igraph::E(gr)$weight, resolution = resolution)
  clust <- leid$membership
  names(clust) <- leid$names
  
  message("Splitting excessively large clusters...")
  # split excessively large clusters - just once:
  for (cid in unique(clust)) {
  genes <- names(clust)[clust == cid]
  if (length(genes) > max_module_size) {
      # subcluster:
      subgraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix[genes, genes], 
                                                      mode = 'undirected', 
                                                      weighted = TRUE)
      subleid <- igraph::cluster_leiden(subgraph, weights = igraph::E(subgraph)$weight, resolution_parameter = resolution * 2)
      subclust <- subleid$membership
      names(subclust) <- subleid$names
      # replace original cluster id:
      clust[genes] = paste0(clust[genes], "subclust", subclust[genes])
  }
  }

  # throw out tiny clusters:
  message("Identified ", length(unique(clust)), " clusters, throwing out clusters with less than ", min_module_size, " genes...")
  clustersizes <- table(clust)
  clusternames <- names(clustersizes)[clustersizes >= min_module_size]
  
  message("Defining modules...")
  modules = list()
  for (cid in clusternames) {
      genes = colnames(adj_matrix)[clust == cid]
      newname <- name_module(adj_matrix[genes, genes])
      modules[[newname]] <- genes
      }

  weights <- define_gene_weights_consensus(modules, weighting = weighting, nh_expression = nh_expression)
  # re-format as data frame:
  weightsdf <- rbindlist(lapply(weights, function(x) data.frame(gene = names(x), weight = x)), idcol="module")


  modules = modules[order(sapply(modules, length), decreasing = TRUE)]
  message("Defined ", length(modules), " modules.")

  # convert to data frame
  modulesdf <- do.call(rbind, lapply(names(modules), function(m) {
    data.frame(module = m, gene = modules[[m]], stringsAsFactors = FALSE)
  }))

  # return:
  out = list(modules = modules,
             modulesdf = modulesdf,
             weights = weights,
             weightsdf = weightsdf)                                            
  return(out)
}

# nh_expression The neighborhood expression matrix. Colnames must contain the genes given in clusts
define_gene_weights_consensus <- function(modules, nh_expression , weighting = "inverse_sqrt") {
  if (!is.element(weighting, c("inverse_sqrt", "inverse", "identity"))) {
    stop("weighting argument must be one of \'inverse_sqrt\', \'inverse\', or \'identity\'")
  }
  if (weighting != "identity") {
    meanneighborhoodexpr <- Matrix::colMeans(nh_expression)
  }

  wts <- list()
  for (i in 1:length(modules)) {
    if (weighting == "inverse_sqrt") {
      wts[[i]] <- meanneighborhoodexpr[modules[[i]]]^-0.5
    }
    if (weighting == "inverse") {
      wts[[i]] <- meanneighborhoodexpr[modules[[i]]]^-1
    }
    if (weighting == "identity") {
      wts[[i]] <- rep(1, length(modules[[i]]))
      names(wts[[i]]) <- modules[[i]]
    }
    wts[[i]] <- wts[[i]] / sum(wts[[i]])
  }
  names(wts) <- names(modules)
  return(wts)
}


# ============================================================================
# Variable Gene Pairs Analysis (Inter-Unit Heterogeneity)
# ============================================================================

#' Find gene pairs with heterogeneous conditional correlations across units
#'
#' @param condcors 3D array [units x genes x genes] of conditional correlations
#' @param high_threshold Correlation above which a pair is "strong" in a unit
#' @param low_threshold Correlation below which a pair is "weak" in a unit
#' @param min_high_units Minimum number of units where pair must be strong
#' @param min_low_units Minimum number of units where pair must be weak
#' @return List with pairs_df, adjacency matrix, and n_pairs
find_variable_gene_pairs <- function(condcors,
                                      high_threshold = 0.4,
                                      low_threshold = 0.05,
                                      min_high_units = 3,
                                      min_low_units = 3) {

  n_units <- dim(condcors)[1]
  n_genes <- dim(condcors)[2]
  gene_names <- dimnames(condcors)[[2]]

  cat("Searching for variable gene pairs across", n_units, "units...\n")
  cat("  High threshold:", high_threshold, "(min", min_high_units, "units)\n")
  cat("  Low threshold:", low_threshold, "(min", min_low_units, "units)\n")

  # Collect variable pairs
  pairs_list <- list()
  pair_count <- 0

  for (i in 1:(n_genes - 1)) {
    if (i %% 500 == 0) cat("  Processing gene", i, "of", n_genes, "\n")
    for (j in (i + 1):n_genes) {
      cors <- condcors[, i, j]
      n_high <- sum(cors > high_threshold, na.rm = TRUE)
      n_low <- sum(cors < low_threshold, na.rm = TRUE)

      if (n_high >= min_high_units && n_low >= min_low_units) {
        pair_count <- pair_count + 1
        pairs_list[[pair_count]] <- data.frame(
          gene1 = gene_names[i],
          gene2 = gene_names[j],
          n_high = n_high,
          n_low = n_low,
          max_cor = max(cors, na.rm = TRUE),
          min_cor = min(cors, na.rm = TRUE),
          range_cor = max(cors, na.rm = TRUE) - min(cors, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (pair_count == 0) {
    cat("No variable gene pairs found with current thresholds.\n")
    cat("Consider lowering min_high_units/min_low_units or adjusting thresholds.\n")
    return(list(
      pairs_df = data.frame(gene1 = character(0), gene2 = character(0),
                            n_high = integer(0), n_low = integer(0),
                            max_cor = numeric(0), min_cor = numeric(0),
                            range_cor = numeric(0)),
      adjacency = matrix(0, 0, 0),
      n_pairs = 0
    ))
  }

  pairs_df <- do.call(rbind, pairs_list)
  pairs_df <- pairs_df[order(-pairs_df$range_cor), ]
  rownames(pairs_df) <- NULL

  # Build adjacency matrix of variable genes
  variable_genes <- unique(c(pairs_df$gene1, pairs_df$gene2))
  n_var_genes <- length(variable_genes)
  adj <- matrix(0, n_var_genes, n_var_genes,
                dimnames = list(variable_genes, variable_genes))
  for (k in 1:nrow(pairs_df)) {
    adj[pairs_df$gene1[k], pairs_df$gene2[k]] <- 1
    adj[pairs_df$gene2[k], pairs_df$gene1[k]] <- 1
  }

  cat("Found", pair_count, "variable gene pairs involving",
      n_var_genes, "unique genes.\n")

  return(list(
    pairs_df = pairs_df,
    adjacency = adj,
    n_pairs = pair_count
  ))
}


#' Extract per-unit correlations for variable gene pairs into a matrix
#'
#' @param condcors 3D array [units x genes x genes] of conditional correlations
#' @param variable_pairs Output from find_variable_gene_pairs()
#' @return Matrix with rows = gene pairs, columns = units
extract_variable_pair_matrix <- function(condcors, variable_pairs) {

  pairs_df <- variable_pairs$pairs_df
  if (nrow(pairs_df) == 0) {
    cat("No variable pairs to extract.\n")
    return(matrix(nrow = 0, ncol = 0))
  }

  unit_names <- dimnames(condcors)[[1]]
  n_pairs <- nrow(pairs_df)
  n_units <- length(unit_names)

  pair_matrix <- matrix(NA, nrow = n_pairs, ncol = n_units,
                        dimnames = list(
                          paste0(pairs_df$gene1, "_", pairs_df$gene2),
                          unit_names
                        ))

  for (k in 1:n_pairs) {
    g1 <- pairs_df$gene1[k]
    g2 <- pairs_df$gene2[k]
    pair_matrix[k, ] <- condcors[, g1, g2]
  }

  return(pair_matrix)
}


#' Cluster variable gene pairs into modules using Leiden algorithm
#'
#' @param variable_pairs Output from find_variable_gene_pairs()
#' @param resolution Leiden clustering resolution
#' @param min_module_size Minimum genes per module
#' @return Named list of gene sets (character vectors)
define_variable_modules <- function(variable_pairs,
                                     resolution = 0.5,
                                     min_module_size = 3) {
  require(igraph)

  adj <- variable_pairs$adjacency
  if (nrow(adj) == 0) {
    cat("No variable pairs — cannot define modules.\n")
    return(list())
  }

  cat("Building variable gene network from", nrow(adj), "genes...\n")
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = NULL)

  cat("Running Leiden clustering (resolution =", resolution, ")...\n")
  leid <- igraph::cluster_leiden(gr, resolution = resolution)
  clust <- leid$membership
  names(clust) <- igraph::V(gr)$name

  # Filter by module size
  cluster_sizes <- table(clust)
  keep_clusters <- names(cluster_sizes)[cluster_sizes >= min_module_size]

  modules <- list()
  for (cid in keep_clusters) {
    genes <- names(clust)[clust == as.integer(cid)]
    # Name the module using the adjacency submatrix
    sub_adj <- adj[genes, genes, drop = FALSE]
    module_name <- name_module(sub_adj)
    modules[[module_name]] <- genes
  }

  # Sort by size (descending)
  modules <- modules[order(sapply(modules, length), decreasing = TRUE)]

  cat("Defined", length(modules), "variable modules",
      "(discarded clusters with <", min_module_size, "genes).\n")
  for (nm in names(modules)) {
    cat("  ", nm, ":", length(modules[[nm]]), "genes\n")
  }

  return(modules)
}


#' Score gene sets using simple fold-enrichment (paper method)
#'
#' Per cell: mean expression of module genes / mean expression of all genes.
#' This is the simple approach from the paper for ad-hoc variable gene sets.
#'
#' @param counts Normalized expression matrix (cells x genes, or genes x cells sparse)
#' @param gene_sets Named list of character vectors (gene sets)
#' @return Matrix: cells x modules
score_gene_sets <- function(counts, gene_sets) {

  if (length(gene_sets) == 0) {
    cat("No gene sets to score.\n")
    return(matrix(nrow = 0, ncol = 0))
  }

  # Handle orientation: ensure genes in rows, cells in columns
  # Check if gene names match rownames or colnames
  test_genes <- gene_sets[[1]]
  if (all(test_genes %in% rownames(counts))) {
    # genes in rows — standard orientation for colMeans approach
    counts_gr <- counts  # genes x cells
  } else if (all(test_genes %in% colnames(counts))) {
    # genes in columns — transpose
    counts_gr <- Matrix::t(counts)
  } else {
    # Find which orientation has more matches
    row_match <- sum(test_genes %in% rownames(counts))
    col_match <- sum(test_genes %in% colnames(counts))
    if (row_match >= col_match) {
      counts_gr <- counts
    } else {
      counts_gr <- Matrix::t(counts)
    }
    warning("Not all genes found in counts matrix. Using best-matching orientation.")
  }

  # Background: mean expression per cell across all genes
  bg <- Matrix::colMeans(counts_gr)
  bg[bg == 0] <- 1e-10  # avoid division by zero

  n_cells <- ncol(counts_gr)
  score_mat <- matrix(NA, nrow = n_cells, ncol = length(gene_sets),
                      dimnames = list(colnames(counts_gr), names(gene_sets)))

  for (i in seq_along(gene_sets)) {
    genes <- gene_sets[[i]]
    genes_present <- genes[genes %in% rownames(counts_gr)]
    if (length(genes_present) == 0) {
      warning("No genes from module '", names(gene_sets)[i],
              "' found in counts matrix.")
      score_mat[, i] <- NA
      next
    }
    if (length(genes_present) < length(genes)) {
      cat("  Module", names(gene_sets)[i], ": using",
          length(genes_present), "of", length(genes), "genes\n")
    }
    module_mean <- Matrix::colMeans(counts_gr[genes_present, , drop = FALSE])
    score_mat[, i] <- module_mean / bg
  }

  return(score_mat)
}


#' Plot variable gene pairs analysis results
#'
#' Multi-page PDF with heatmap of per-unit correlations and spatial plots.
#'
#' @param pair_matrix Matrix from extract_variable_pair_matrix()
#' @param variable_modules Named list of gene sets from define_variable_modules()
#' @param scores Score matrix from score_gene_sets()
#' @param xy Spatial coordinates matrix (cells x 2)
#' @param metadata Cell metadata data.frame
#' @param region_col Column name for tissue/region
#' @param annotation_col Column name for cell type annotation
#' @param output_file Path to output PDF
plot_variable_analysis <- function(pair_matrix, variable_modules,
                                    scores, xy, metadata,
                                    region_col, annotation_col,
                                    output_file) {

  require(ComplexHeatmap)
  require(viridis)
  require(circlize)
  require(grid)

  if (file.exists(output_file)) file.remove(output_file)

  pdf(output_file, width = 11, height = 8.5, onefile = TRUE)

  tryCatch({

    # --- Page 1: Summary text ---
    grid.newpage()
    summary_text <- paste0(
      "Variable Gene Pairs Analysis\n\n",
      "Gene pairs with heterogeneous conditional correlations across units.\n\n",
      "Number of variable gene pairs: ", nrow(pair_matrix), "\n",
      "Number of units: ", ncol(pair_matrix), "\n",
      "Number of variable modules: ", length(variable_modules), "\n\n",
      "Module sizes:\n",
      paste(sapply(names(variable_modules), function(nm) {
        paste0("  ", nm, ": ", length(variable_modules[[nm]]), " genes")
      }), collapse = "\n")
    )
    grid.text(summary_text, x = 0.5, y = 0.5,
              gp = gpar(fontsize = 14, fontfamily = "mono"),
              just = "centre")

    # --- Page 2: Heatmap of per-unit correlations ---
    if (nrow(pair_matrix) > 0 && ncol(pair_matrix) > 0) {
      col_fun <- colorRamp2(c(-0.2, 0, 0.2, 0.5, 0.8),
                            c("blue", "white", "lightyellow", "orange", "red"))

      # Limit label display for large pair counts
      show_row_names <- nrow(pair_matrix) <= 80

      ht <- Heatmap(pair_matrix,
                    name = "Cond. Cor.",
                    col = col_fun,
                    column_title = "Per-Unit Conditional Correlations of Variable Gene Pairs",
                    row_names_gp = gpar(fontsize = if (show_row_names) max(4, 8 - nrow(pair_matrix) / 20) else 6),
                    column_names_gp = gpar(fontsize = 8),
                    show_row_names = show_row_names,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE,
                    row_title = paste0(nrow(pair_matrix), " gene pairs"),
                    column_title_gp = gpar(fontsize = 14, fontface = "bold"))
      draw(ht)
    }

    # --- Pages 3+: Spatial plots per variable module ---
    if (length(variable_modules) > 0 && nrow(scores) > 0) {
      units <- unique(metadata[[region_col]])

      for (mod_idx in seq_along(variable_modules)) {
        mod_name <- names(variable_modules)[mod_idx]
        mod_genes <- variable_modules[[mod_idx]]

        if (!mod_name %in% colnames(scores)) next

        mod_scores <- scores[, mod_name]
        # Cap at 99.5th percentile for color scaling
        cap <- quantile(mod_scores, 0.995, na.rm = TRUE)
        if (cap == 0) cap <- 1

        grid.newpage()

        # Title page for module
        grid.text(
          paste0("Variable Module: ", mod_name, "\n",
                 "Genes: ", paste(mod_genes, collapse = ", ")),
          x = 0.5, y = 0.5,
          gp = gpar(fontsize = 16, fontface = "bold"),
          just = "centre"
        )

        # Spatial plot: all units together
        grid.newpage()
        scaled <- pmin(mod_scores / cap, 1)
        cols <- viridis::viridis(101, option = "B")[1 + round(100 * scaled)]

        pushViewport(viewport(x = 0.5, y = 0.45, width = 0.9, height = 0.8))

        # Normalize coordinates to viewport
        x_range <- range(xy[, 1], na.rm = TRUE)
        y_range <- range(xy[, 2], na.rm = TRUE)
        x_norm <- (xy[, 1] - x_range[1]) / diff(x_range)
        y_norm <- (xy[, 2] - y_range[1]) / diff(y_range)

        grid.points(x = unit(x_norm, "npc"),
                    y = unit(y_norm, "npc"),
                    pch = 16, size = unit(0.3, "mm"),
                    gp = gpar(col = cols))
        popViewport()

        grid.text(paste0(mod_name, " — Spatial Scores"),
                  x = 0.5, y = 0.95,
                  gp = gpar(fontsize = 14, fontface = "bold"))

        # Per-unit spatial plots (grid layout)
        n_units <- length(units)
        if (n_units > 0) {
          ncols <- min(4, n_units)
          nrows <- ceiling(n_units / ncols)

          # Limit to 2 pages of per-unit plots
          max_per_page <- ncols * min(nrows, 4)
          pages_needed <- ceiling(n_units / max_per_page)

          for (page in 1:pages_needed) {
            grid.newpage()
            grid.text(paste0(mod_name, " — Per-Unit Scores (page ", page, ")"),
                      x = 0.5, y = 0.97,
                      gp = gpar(fontsize = 12, fontface = "bold"))

            start_idx <- (page - 1) * max_per_page + 1
            end_idx <- min(page * max_per_page, n_units)
            page_units <- units[start_idx:end_idx]
            n_page <- length(page_units)
            page_nrows <- ceiling(n_page / ncols)

            for (u_idx in seq_along(page_units)) {
              u <- page_units[u_idx]
              row_pos <- ceiling(u_idx / ncols)
              col_pos <- ((u_idx - 1) %% ncols) + 1

              cell_mask <- metadata[[region_col]] == u
              if (sum(cell_mask) == 0) next

              x_vp <- (col_pos - 0.5) / ncols
              y_vp <- 1 - (row_pos - 0.5) / page_nrows * 0.9 - 0.05

              pushViewport(viewport(x = x_vp, y = y_vp,
                                    width = 0.9 / ncols,
                                    height = 0.85 / page_nrows))

              xy_u <- xy[cell_mask, , drop = FALSE]
              scores_u <- mod_scores[cell_mask]
              scaled_u <- pmin(scores_u / cap, 1)
              cols_u <- viridis::viridis(101, option = "B")[1 + round(100 * scaled_u)]

              x_r <- range(xy_u[, 1], na.rm = TRUE)
              y_r <- range(xy_u[, 2], na.rm = TRUE)
              if (diff(x_r) == 0) x_r <- x_r + c(-1, 1)
              if (diff(y_r) == 0) y_r <- y_r + c(-1, 1)
              xn <- (xy_u[, 1] - x_r[1]) / diff(x_r)
              yn <- (xy_u[, 2] - y_r[1]) / diff(y_r)

              grid.points(x = unit(xn, "npc"), y = unit(yn, "npc"),
                          pch = 16, size = unit(0.2, "mm"),
                          gp = gpar(col = cols_u))
              grid.text(as.character(u), x = 0.5, y = 1.05,
                        gp = gpar(fontsize = 7, fontface = "bold"))
              grid.rect(gp = gpar(col = "grey80", fill = NA, lwd = 0.5))

              popViewport()
            }
          }
        }
      }
    }

    dev.off()

    if (file.exists(output_file) && file.size(output_file) > 0) {
      cat("Created:", output_file, "\n")
      cat("File size:", round(file.size(output_file) / 1024, 2), "KB\n")
    }

  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    stop("Error creating variable analysis PDF: ", e$message)
  })
}


summarize_insitucor_consensus <- function(consensus, modules_df, output_file = "insitucor_consensus_summary.pdf", vertex_size = 2, condcors = NULL, method = NULL, threshold = NULL, min_fraction = NULL, min_module_size = NULL, max_module_size = NULL, resolution = NULL) {
  require(igraph)
  require(ggrepel)
  require(ggraph)

  # Print cores included, print params used!
  
  adj_matrix <- as.matrix(consensus$adjacency)

  tryCatch({
    # Start PDF device
    pdf(output_file, width = 11, height = 8.5, onefile = TRUE)
    
    text_lines <- c(
    "InSituCor Consensus Results Summary",
    ""
  )

  ## condcors block
  if (!is.null(condcors)) {
    n_cores <- dim(condcors)[1]
    core_names <- dimnames(condcors)[[1]]

    text_lines <- c(
      text_lines,
      paste0("Consensus Network built for ", n_cores, " cores:"),
      paste(core_names, collapse = ", "),
      ""
    )
  }

  ## method block
  if (!is.null(method)) {
    text_lines <- c(
      text_lines,
      paste0("Method used to build consensus network: ", method),
      ""
    )
  }

  ## threshold block
  if (!is.null(threshold)) {
    text_lines <- c(
      text_lines,
      paste0("Conditional Correlation Threshold: ", threshold),
      ""
    )
  }

  ## min_fraction block
  if (!is.null(min_fraction) && !is.null(condcors)) {
    quantile <- 1 - min_fraction
    required_cores <- ceiling(n_cores * min_fraction)
    min_fraction <- round(min_fraction, 3)
    quantile <- round(quantile, 3)

    text_lines <- c(
      text_lines,
      paste0(
        "Minimum fraction requirement: ",
        min_fraction,
        " → quantile : ",
        quantile,
        " (edge must be present in ≥ ",
        required_cores,
        " cores)"
      ),
      ""
    )
  }
  if (!is.null(resolution) && !is.null(min_module_size) && !is.null(max_module_size)) {
    text_lines <- c(
    text_lines,
    paste0(
      "Module definition parameters: ",
      "resolution = ", resolution, ", ",
      "minimum module size = ", min_module_size, ", ",
      "maximum module size = ", max_module_size
    ),
    ""
  )
  }

  ## module summary
  module_dims <- dim(modules_df)
  n_unique_modules <- length(unique(modules_df$module))

  text_lines <- c(
    text_lines,
    "Module Summary:",
    paste0("Module dataframe dimensions: ",
           module_dims[1], " rows × ", module_dims[2], " columns"),
    paste0("Number of unique modules: ", n_unique_modules)
  )

  stats_text <- textGrob(
    paste(text_lines, collapse = "\n"),
    x = 0.01,
    y = 0.99,
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )

  grid.newpage()
  grid.draw(stats_text)

    
    genes_per_module <- table(modules_df$module)
    genes_df <- data.frame(
      module = names(genes_per_module),
      n_genes = as.numeric(genes_per_module)
    )
    
    p_barplot <- ggplot(genes_df, aes(x = reorder(module, n_genes), y = n_genes)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Number of Genes per Module",
           x = "Module",
           y = "Number of Genes") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.y = element_text(size = 10)
      )

    print(p_barplot)

    
    # Determine how many rows to show (max 50 for readability)
    n_rows_to_show <- min(25, nrow(modules_df))
    table_data <- head(modules_df, n_rows_to_show)
    
    # Create table grob
    table_grob <- tableGrob(
      table_data,
      rows = NULL,
      theme = ttheme_default(
        core = list(fg_params = list(cex = 0.6)),
        colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
      )
    )
    
    # Add title and footer
    title_grob <- textGrob(
      "Module Table", 
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    footer_grob <- textGrob(
      paste0("Showing first ", n_rows_to_show, " of ", nrow(modules_df), " total rows"),
      gp = gpar(fontsize = 10, fontface = "italic")
    )
    
    grid.arrange(
      title_grob,
      table_grob,
      footer_grob,
      heights = c(0.1, 0.85, 0.05),
      ncol = 1
    )

    message("Creating Network Plots..")
    genes <- modules_df$gene
    modules_vec <- modules_df$module

    # set module colors:
    distinctcolors <-c('#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5', 
                        '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
                        '#e377c2','#17becf','#7f7f7f','#E41A1C','#377EB8',
                        '#BEBADA','#FB8072','#80B1D3','#FDB462','#984EA3',
                        '#A65628','#F781BF','#999999','#FF7F00','#BC80BD',
                        '#8DD3C7','#BEBADA','#FB8072','#80B1D3','#B3DE69',
                        '#BC80BD','#66C2A5','#FC8D62','#8DA0CB','#E78AC3',
                        '#A6D854','#FFD92F','#E5C494','#CCEBC5','#FCCDE5',
                        '#8c564b','#4DAF4A',
                        sample(colors()[!grepl("grey", colors())], 200, replace = FALSE))
    distinctcolors <- distinctcolors[seq_len(length(unique(modules_vec)) + 1)]
    names(distinctcolors) <- c("", unique(modules_vec))

    # get a good layout:
    xum <- uwot::umap(adj_matrix[genes, genes], 
                    spread = 25, 
                    min_dist = 0.1,
                    n_neighbors = max(min(length(genes)-2, 15), 1))
    plot(xum, pch = 16, col = distinctcolors[modules_vec])

    gr0 <- igraph::graph_from_adjacency_matrix(adj_matrix[genes, genes], 
                                    weighted = TRUE,
                                    mode = "undirected")
    #  color by module:
    igraph::V(gr0)$module <- modules_vec[match(genes, names(igraph::V(gr0)))]
    V(gr0)$color <- distinctcolors[igraph::V(gr0)$module]
    V(gr0)$size <- vertex_size                            
    
    # Umap based, without labels
    igraph::plot.igraph(gr0, 
                  layout = xum, 
                  vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + FALSE]), 
                  vertex.size = 0, #edge.size = 1, 
                  vertex.color = distinctcolors[igraph::V(gr0)$module], 
                  vertex.label.color = distinctcolors[igraph::V(gr0)$module]
                  )

    # Umap based, with labels
    igraph::plot.igraph(gr0, 
                  layout = xum, 
                  vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), 
                  vertex.size = vertex_size, #edge.size = 1, 
                  vertex.color = distinctcolors[igraph::V(gr0)$module], 
                  vertex.label.color = distinctcolors[igraph::V(gr0)$module]
                  )

    # Base graph plot FR
    plot(gr0, 
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       edge.width = 0.5,
       layout = layout_with_fr(gr0)) 

    # Advanced graph plot FR
    p1 <- ggraph(gr0, layout = "fr") +
      geom_edge_link(width = 0.5, alpha = 0.4) +
      geom_node_point(
        aes(color = I(color), size = I(size))
      ) +
      geom_node_text(
        # aes(label = name, color = I(color)),
        aes(label = name, color = I(color)),
        size = 2,
        repel = TRUE
      ) +
      scale_size_identity() +
      theme_void()
    
    print(p1)
    
    dev.off()
    
    # Verify file was created
    if (file.exists(output_file) && file.size(output_file) > 0) {
      cat("✓ Summary report created successfully:", output_file, "\n")
      cat("✓ File size:", round(file.size(output_file) / 1024, 2), "KB\n")
    } else {
      warning("PDF file may be corrupted or empty")
    }
    
    # Return summary statistics
    invisible(list(
      dimensions = module_dims,
      n_modules = n_unique_modules,
      genes_per_module = genes_df,
      output_file = output_file
    ))
    
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    stop("Error creating PDF: ", e$message)
  })
}


pathway_summary_plots <- function(pw_cols, seu, reduction, cluster_cols = NULL, output_file) {

  pdf(output_file, width = 10, height = 10)

  for (pw in pw_cols) {
    sp_plot_full <- xyplot(pw, metadata = seu@meta.data, ptsize = 0.01) + coord_fixed() +
      labs(title = pw, color = "AUCell_score")
    print(sp_plot_full)
  }

  for (pw in pw_cols) {
    reduction_plot_full <- plot_embedding(
              seu,
              reduction = reduction,
              group.by = pw,
              label = TRUE,
              palette = NULL,
              legend = TRUE
              ) +
        labs(title = pw, color = "AUCell_score")
    print(reduction_plot_full)
  }
  for (cluster_col in cluster_cols) {
    if (!is.null(cluster_cols)) {
      for (pw in pw_cols) {
      p <- plot_violin_box(seu@meta.data, pw, cluster_col, plot_type = "box", color_pal = NULL)
      p <- plot_violin_box(seu@meta.data, pw, cluster_col, plot_type = "violin", color_pal = NULL)
      print(p)
    }
    }
  }

  # Create named vector: 
  pw_labels <- setNames(paste0("pw", seq_along(pw_cols)), pw_cols)

  for (cluster_col in cluster_cols) {
    if (!is.null(cluster_cols)) {
      
      dp <- DotPlot(object = seu, features = pw_cols, group.by = cluster_col, cols = c("white", "red"), scale = FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
        scale_x_discrete(labels = pw_labels)
      print(dp)

    }
  }

  dev.off()
}

plot_violin_box <- function(auc_df, pathway, meta_col, 
                                test_pairs = NULL, plot_type = "violin",
                                color_pal = NULL) {
  
  # Calculate median per group
  group_medians <- auc_df %>%
    group_by(!!sym(meta_col)) %>%
    summarise(median_score = median(.data[[pathway]], na.rm = TRUE), .groups = "drop")
  
  # Set up colors
  levels_meta <- if(is.factor(auc_df[[meta_col]])) {
    levels(auc_df[[meta_col]])
  } else {
    unique(auc_df[[meta_col]])
  }
  
  cols_vec <- if(!is.null(color_pal)) {
    color_pal[levels_meta]
  } else {
    NULL  # ggplot default colors
  }
  
  # Base plot
  p <- ggplot(auc_df, aes(x = .data[[meta_col]], y = .data[[pathway]], 
                          fill = .data[[meta_col]])) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = pathway,
      x = meta_col,
      y = "AUCell Score"
    )
  
  # Add violin or boxplot
  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7, trim = FALSE)
  } else {
    p <- p + geom_boxplot(alpha = 0.7, outlier.size = 0.5)
  }
  
  
  # Apply color palette if provided
  if (!is.null(cols_vec)) {
    p <- p + scale_fill_manual(values = cols_vec)
  }
  
  return(p)
}



# Pathways

#' Calculate summary statistics for all pathways
#' 
#' @param auc_df Data frame with cells x pathways (+ metadata columns)
#' @param pw_cols Vector of pathway column names
#' @return Data frame with pathway statistics
calculate_pathway_stats <- function(auc_df, pw_cols) {
  
  pathway_stats <- data.frame(
    pathway = pw_cols,
    mean = colMeans(auc_df[, pw_cols], na.rm = TRUE),
    median = apply(auc_df[, pw_cols], 2, median, na.rm = TRUE),
    sd = apply(auc_df[, pw_cols], 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Add derived metrics
  pathway_stats$cv <- pathway_stats$sd / ifelse(pathway_stats$mean == 0, NA, pathway_stats$mean )
  pathway_stats$min <- apply(auc_df[, pw_cols], 2, min, na.rm = TRUE)
  pathway_stats$max <- apply(auc_df[, pw_cols], 2, max, na.rm = TRUE)
  pathway_stats$range <- pathway_stats$max - pathway_stats$min
  
  # Add percentile ranks
  pathway_stats$mean_percentile <- rank(pathway_stats$mean) / nrow(pathway_stats)
  pathway_stats$sd_percentile <- rank(pathway_stats$sd) / nrow(pathway_stats)
  pathway_stats$cv_percentile <- rank(pathway_stats$cv) / nrow(pathway_stats)
  
  return(pathway_stats)
}

visualize_pathway_distributions <- function(pathway_stats, min_mean, min_sd, output_file, min_cv = NULL) {

    # Mean Histogram
    p1 <- ggplot(pathway_stats, aes(x = mean)) +
      geom_histogram(bins = 50, fill = "cornflowerblue", alpha = 0.7) +
      geom_vline(xintercept = min_mean, 
                color = "red", linetype = "dashed") +
      labs(title = "Distribution of Pathway Mean Scores",
          x = "Mean AUCell Score", y = "Count") +
      theme_bw()

    # SD Histogram
    p2 <- ggplot(pathway_stats, aes(x = sd)) +
      geom_histogram(bins = 50, fill = "chartreuse4", alpha = 0.7) +
      geom_vline(xintercept = min_sd, 
                color = "red", linetype = "dashed") +
      labs(title = "Distribution of Pathway Standard Deviations",
          x = "SD", y = "Count") +
      theme_bw()

    # # CV Histogram
    # p3 <- ggplot(pathway_stats, aes(x = cv)) +
    #   geom_histogram(bins = 50, fill = "orchid3", alpha = 0.7) +
    #   geom_vline(xintercept = min_cv, 
    #             color = "red", linetype = "dashed") +
    #   labs(title = "Distribution of Coefficient of Variation",
    #       x = "CV (SD/Mean)", y = "Count") +
    #   theme_bw() +
    #   xlim(0, quantile(pathway_stats$cv, 0.99, na.rm = TRUE))  # Remove extreme outliers

    # Mean-SD scatterplot
    p4 <- ggplot(pathway_stats, aes(x = mean, y = sd)) +
      geom_point(alpha = 0.4, size = 2, color = "darkslategray") +
      geom_hline(yintercept = min_sd, 
                color = "red", linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = min_mean, 
                color = "red", linetype = "dashed", alpha = 0.7) +
      labs(title = "Pathway Mean vs SD",
          subtitle = "Red lines: min SD, min mean",
          x = "Mean Score", y = "Standard Deviation") +
      theme_bw()

    # Ranked Variance Plot
    pathway_stats_sorted <- pathway_stats[order(pathway_stats$mean, decreasing = TRUE), ]
    pathway_stats_sorted$rank <- 1:nrow(pathway_stats_sorted)
    
    p5 <- ggplot(pathway_stats_sorted, aes(x = rank, y = mean)) +
      geom_line(color = "cornflowerblue", size = 1) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_hline(yintercept = min_mean, 
                color = "red", linetype = "dashed") +
      labs(title = "Pathway Mean Ranking",
          subtitle = "Red line: Min Mean",
          x = "Rank", y = "Mean AUCell Score") +
      theme_bw()
    
    # Ranked Variance Plot
    pathway_stats_sorted <- pathway_stats[order(pathway_stats$sd, decreasing = TRUE), ]
    pathway_stats_sorted$rank <- 1:nrow(pathway_stats_sorted)
    
    p6 <- ggplot(pathway_stats_sorted, aes(x = rank, y = sd)) +
      geom_line(color = "chartreuse4", size = 1) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_hline(yintercept = min_sd, 
                color = "red", linetype = "dashed") +
      labs(title = "Pathway Variance Ranking",
          subtitle = "Red line: Min SD",
          x = "Rank", y = "Standard Deviation") +
      theme_bw()

    # Ranked mean bar plot (top 50 for readability)
    pathway_stats_mean_sorted <- pathway_stats[order(pathway_stats$mean, decreasing = TRUE), ]
    pathway_stats_mean_sorted$rank <- 1:nrow(pathway_stats_mean_sorted)
    pathway_stats_mean_sorted$pathway_short <- gsub("^REACTOME_|^HALLMARK_", "", pathway_stats_mean_sorted$pathway)

    p7 <- ggplot(pathway_stats_mean_sorted, aes(x = reorder(pathway_short, mean), y = mean)) +
      geom_bar(stat = "identity", fill = "cornflowerblue", alpha = 0.8) +
      geom_hline(yintercept = min_mean, color = "red", linetype = "dashed") +
      coord_flip() +
      labs(title = "Pathways by Mean Score",
          x = NULL, y = "Mean AUCell Score") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 7))

    pdf(output_file, width = 10, height = 10)
    print(p1)
    print(p2)
    print(p4)
    print(p5)
    print(p6)
    print(p7)
    dev.off()
    
}

filter_pathways_basic <- function(pathway_stats, 
                                   min_mean = NULL,
                                   min_sd = NULL,
                                   min_cv = NULL) {
  
  if (is.null(min_mean)) min_mean <- quantile(pathway_stats$mean, 0.25)
  if (is.null(min_sd)) min_sd <- median(pathway_stats$sd)
  message(sprintf("Using absolute thresholds: mean > %.3f, sd > %.3f",
                  min_mean, min_sd))
  
  # Apply filters
  pass_mean <- pathway_stats$mean >= min_mean
  pass_sd <- pathway_stats$sd >= min_sd
  
  # Optional CV filter
  if (!is.null(min_cv)) {
    pass_cv <- pathway_stats$cv >= min_cv
    pass_filter <- pass_mean & pass_sd & pass_cv
    message(sprintf("Pathways passing filters: %d/%d (mean & sd & cv)",
                    sum(pass_filter), nrow(pathway_stats)))
  } else {
    pass_filter <- pass_mean & pass_sd
    message(sprintf("Pathways passing filters: %d/%d (mean & sd)",
                    sum(pass_filter), nrow(pathway_stats)))
  }
  
  filtered_pathways <- pathway_stats$pathway[pass_filter]
  
  return(filtered_pathways)
}

#' Merge filtered pathways with manually curated list
#' 
#' @param filtered_pathways Vector of pathways from automated filtering
#' @param manual_pathways Vector of manually selected pathways
#' @param all_pathways Vector of all available pathways (for validation)
#' @return List with final pathways and metadata
merge_pathway_lists <- function(filtered_pathways, manual_pathways, all_pathways) {
  
  # Validate manual pathways exist
  valid_manual <- manual_pathways[manual_pathways %in% all_pathways]
  invalid_manual <- setdiff(manual_pathways, all_pathways)
  
  if (length(invalid_manual) > 0) {
    warning(sprintf("Manual pathways not found in data: %s", 
                    paste(invalid_manual, collapse = ", ")))
  }
  
  # Combine lists (union)
  final_pathways <- unique(c(filtered_pathways, valid_manual))
  
  # Categorize pathways
  only_filtered <- setdiff(filtered_pathways, valid_manual)
  only_manual <- setdiff(valid_manual, filtered_pathways)
  both <- intersect(filtered_pathways, valid_manual)
  
  message(sprintf("Final pathway count: %d", length(final_pathways)))
  message(sprintf("  - From filtering only: %d", length(only_filtered)))
  message(sprintf("  - From manual list only: %d", length(only_manual)))
  message(sprintf("  - In both lists: %d", length(both)))
  
  return(list(
    final_pathways = final_pathways,
    only_filtered = only_filtered,
    only_manual = only_manual,
    both = both,
    invalid_manual = invalid_manual
  ))
}


pathway_summary_plots_spatial <- function(pw_cols, seu, reduction, cluster_cols = NULL, output_file) {

  pdf(output_file, width = 10, height = 10)

  for (pw in pw_cols) {
    sp_plot_full <- xyplot(pw, metadata = seu@meta.data, ptsize = 0.01) + coord_fixed() +
      labs(title = pw, color = "AUCell_score")
    print(sp_plot_full)
  }

  # Create named vector: 
  pw_labels <- setNames(paste0("pw", seq_along(pw_cols)), pw_cols)

  for (cluster_col in cluster_cols) {
    if (!is.null(cluster_cols)) {
      
      dp <- DotPlot(object = seu, features = pw_cols, group.by = cluster_col, cols = c("white", "red"), scale = FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
        scale_x_discrete(labels = pw_labels)
      print(dp)

    }
  }

  dev.off()
}

pathway_summary_plots_aggr <- function(pw_cols, auc_df, cluster_col = NULL, output_file) {

  plot_title <- "PW Enrichment"
  
  pdf(output_file, width = 10, height = 10)

  # Aggregate mean pathway scores per cluster
  agg_data <- aggregate(auc_df[, pw_cols], 
                        by = list(cluster = auc_df[[cluster_col]]), 
                        FUN = mean)
  rownames(agg_data) <- agg_data$cluster
  agg_data <- agg_data[, -1]

  # --- Plot 1: Raw AUCell scores ---
  # Determine color scale limits from 95th percentile
  abs_vals <- abs(as.matrix(agg_data))
  limit_raw <- quantile(abs_vals, 0.95, na.rm = TRUE)
  
  col_fun_raw <- colorRamp2(c(0, limit_raw/2, limit_raw), 
                            c("blue", "white", "red"))
  
  h1 <- Heatmap(as.matrix(agg_data),
                name = "AUCell\nScore",
                col = col_fun_raw,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_dend_side = "right",
                show_row_dend = TRUE,
                column_names_rot = 90,
                row_names_side = "left",
                column_names_side = "bottom",
                column_names_max_height = unit(16, "cm"),
                column_title = paste0(plot_title, " - Raw AUCell Scores"),
                heatmap_legend_param = list(
                  at = c(0, limit_raw/2, limit_raw),
                  labels = c("0", sprintf("%.2f", limit_raw/2), sprintf("%.2f", limit_raw))
                ))
  
  # --- Plot 2: Z-scaled scores ---
  # Z-scale per pathway (across cell types)
  agg_scaled <- scale(agg_data)  # scale() operates on columns by default
  
  col_fun_zscore <- colorRamp2(c(-2, 0, 2), 
                               c("blue", "white", "red"))
  
  h2 <- Heatmap(agg_scaled,
                name = "Z-score",
                col = col_fun_zscore,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                row_dend_side = "right",
                show_row_dend = TRUE,
                column_names_rot = 90,
                row_names_side = "left",
                column_names_side = "bottom",
                column_names_max_height = unit(16, "cm"),
                column_title = paste0(plot_title, " - Z-scaled (Relative Enrichment)"),
                heatmap_legend_param = list(
                  at = c(-2, -1, 0, 1, 2)
                ))
  
  # Draw both heatmaps
  print(h1)
  grid.newpage()
  print(h2)

  dev.off()
  
  # Return aggregated data for further use
  return(list(raw = agg_data, scaled = agg_scaled))
}


combine_llm_runs <- function(runs_parent_dir, anonmap_filepath) {

    # Load anonymization map
    anon_map_dt <- NULL
    if (!is.null(anonmap_filepath) && file.exists(anonmap_filepath)) {
    anon_map_dt <- fread(anonmap_filepath)
    # expect columns: anonym_clust, semisup_clust
    stopifnot(all(c("anonym_clust", "semisup_clust") %in% names(anon_map_dt)))
    message(sprintf("Loaded anon map: %d clusters", nrow(anon_map_dt)))
    }

    # Discover run directories
    run_dirs <- list.dirs(runs_parent_dir, full.names = TRUE, recursive = FALSE)
    message(sprintf("Found %d run directories", length(run_dirs)))

    # Build master list of all cluster names
    if (!is.null(anon_map_dt)) {
    all_IST_clusters <- anon_map_dt$semisup_clust
    } else {
    # fallback: collect from all label files
    all_IST_clusters <- unique(unlist(lapply(run_dirs, function(rd) {
        lfs <- list.files(rd, pattern = "_labels\\.rds$", full.names = TRUE)
        unique(unlist(lapply(lfs, function(lf) names(readRDS(lf)))))
    })))
    }

    # Extract the labels from each run

    all_labels <- rbindlist(lapply(run_dirs, function(rd) {

    run_id   <- basename(rd)
    is_anon  <- grepl("_anonTRUE", run_id, fixed = TRUE)

    label_files <- list.files(rd, pattern = "_labels\\.rds$", full.names = TRUE)

    if (length(label_files) == 0) {
        message(sprintf("  [%s] No label files found, skipping.", run_id))
        return(NULL)
    }

    rbindlist(lapply(label_files, function(lf) {

        llm_name <- sub("_labels\\.rds$", "", basename(lf))
        labs <- readRDS(lf)

        if (is.null(names(labs)) || length(labs) == 0) {
        message(sprintf("  [%s/%s] Empty or unnamed vector, skipping.", run_id, llm_name))
        return(NULL)
        }

        dt <- data.table(
        cluster_raw = names(labs),
        label       = unname(labs),
        llm_name    = llm_name,
        run_id      = run_id,
        is_anon     = is_anon
        )

        # back-translate anonymized cluster names to IST names
        if (is_anon && !is.null(anon_map_dt)) {
        dt <- merge(dt, anon_map_dt,
                    by.x = "cluster_raw", by.y = "anonym_clust",
                    all.x = TRUE)
        if (any(is.na(dt$semisup_clust))) {
            unmatched <- dt$cluster_raw[is.na(dt$semisup_clust)]
            warning(sprintf("  [%s/%s] Could not map anonymized clusters: %s",
                            run_id, llm_name, paste(unmatched, collapse = ", ")))
        }
        dt[, IST_clust := semisup_clust]
        dt[, semisup_clust := NULL]
        } else {
        dt[, IST_clust := cluster_raw]
        }

        # column name = full run_id + llm_name
        dt[, source := paste0(run_id, "__", llm_name)]

        dt[, .(IST_clust, cluster_raw, label, source, is_anon, run_id)]
    }))
    }))

    message(sprintf("Collected %d total label entries across %d sources",
                    nrow(all_labels), length(unique(all_labels$source))))

    # == 5. Expand to full cluster × source grid (fill missing with NA) ===========
    #       This handles LLM runs that only labelled a subset of clusters.

    all_sources <- unique(all_labels$source)

    full_grid <- CJ(IST_clust = all_IST_clusters, source = all_sources)

    all_labels <- merge(full_grid, all_labels[, .(IST_clust, source, label)],
                        by = c("IST_clust", "source"),
                        all.x = TRUE)
    # missing entries now have label = NA

    # == 6. Add anonym_clust column ================================================

    if (!is.null(anon_map_dt)) {
    all_labels <- merge(all_labels, anon_map_dt,
                        by.x = "IST_clust", by.y = "semisup_clust",
                        all.x = TRUE)
    } else {
    all_labels[, anonym_clust := NA_character_]
    }

    # == 7. Pivot to wide format ===================================================

    comparison_wide <- dcast(all_labels, IST_clust + anonym_clust ~ source,
                            value.var = "label",
                            fill = NA_character_)

    # sort: known clusters first (multi-char names), then unknown (single lowercase)
    comparison_wide[, is_unknown := nchar(IST_clust) == 1 & grepl("^[a-z]$", IST_clust)]
    setorder(comparison_wide, is_unknown, IST_clust)
    comparison_wide[, is_unknown := NULL]

    return(comparison_wide[])

    }


# =============================================================================
# smiDE Helper Functions
# =============================================================================

#' Extract and format smiDE pairwise results into a clean data.frame with FDR
#'
#' @param de_obj     smiDE result object from smi_de()
#' @param de_var     character, the DE variable name
#' @param contrast   character, specific contrast to extract (NULL = all)
#' @param fdr_method character, p.adjust method (default "BH")
#' @return data.frame with columns: gene, contrast, FC, log2FC, p_raw, fdr
extract_smide_results <- function(de_obj, de_var, contrast = NULL,
                                  fdr_method = "BH") {
    res_pw <- smiDE::results(de_obj,
                             comparisons = "pairwise",
                             variable = de_var)$pairwise

    if (!is.null(contrast)) {
        res_pw <- res_pw[res_pw$contrast == contrast, ]
    }

    res_df <- data.frame(
        gene     = res_pw$target,
        contrast = as.character(res_pw$contrast),
        FC       = res_pw$fold_change,
        log2FC   = log2(res_pw$fold_change),
        p_raw    = res_pw$p.value,
        stringsAsFactors = FALSE
    )

    # BH correction per contrast
    res_df <- do.call(rbind, lapply(split(res_df, res_df$contrast), function(x) {
        x$fdr <- p.adjust(x$p_raw, method = fdr_method)
        x
    }))
    rownames(res_df) <- NULL
    return(res_df)
}


#' Select top N genes for spatial DE based on effect size and p-value
#'
#' @param results_df data.frame with gene, effect/log2FC, and p-value columns
#' @param effect_col character, column name for effect size
#' @param p_col      character, column name for p-value
#' @param top_n      integer, number of top genes (default 10)
#' @param p_quantile numeric, p-value quantile threshold for pre-filtering (default 0.01)
#' @return character vector of top gene names
select_top_genes_for_spatial <- function(results_df, effect_col = "effect",
                                         p_col = "p_raw", top_n = 10,
                                         p_quantile = 0.01) {
    p_threshold <- quantile(results_df[[p_col]], p_quantile, na.rm = TRUE)
    ranking_score <- ifelse(results_df[[p_col]] <= p_threshold,
                            abs(results_df[[effect_col]]), 0)
    top_idx <- order(ranking_score, decreasing = TRUE)[
        seq_len(min(top_n, nrow(results_df)))
    ]
    return(results_df$gene[top_idx])
}


#' Plot a single gene's expression on spatial tissue coordinates
#'
#' @param seu_obj    Seurat object (subsetted to cell type of interest)
#' @param gene       character, gene name
#' @param x_col      character, x coordinate column (default "x_slide_mm")
#' @param y_col      character, y coordinate column (default "y_slide_mm")
#' @param data_layer character, Seurat data layer (default "data")
#' @param pal_option character, viridis palette option (default "plasma")
#' @param ptsize     numeric, point size (default 0.6)
#' @return ggplot object
plot_gene_spatial_expression <- function(seu_obj, gene,
                                         x_col = "x_slide_mm",
                                         y_col = "y_slide_mm",
                                         data_layer = "data",
                                         pal_option = "plasma",
                                         ptsize = 0.6) {
    expr_vals <- Seurat::GetAssayData(seu_obj, layer = data_layer)[gene, ]
    plot_md <- data.table::as.data.table(seu_obj@meta.data)
    plot_md[[gene]] <- expr_vals

    p <- xyplot(
        cluster_column = gene,
        metadata = plot_md,
        x_column = x_col,
        y_column = y_col,
        ptsize = ptsize,
        alphasize = 0.8,
        continuous_palette = function(n) viridis::viridis(n, option = pal_option),
        theme = ggplot2::theme_minimal() +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    ) + ggplot2::labs(title = paste0(gene, " (", data_layer, ")"))

    return(p)
}


#' Plot a heatmap of mean expression per bin for selected genes
#'
#' @param seu_obj   Seurat object (subsetted to cell type of interest)
#' @param bin_col   character, metadata column defining bins/groups
#' @param genes     character vector of genes (NULL = all genes)
#' @param bin_order character vector specifying bin order (NULL = default)
#' @param title     character, plot title
#' @param scale     character, pheatmap scale param (default "row")
#' @return pheatmap object
plot_mean_expression_heatmap <- function(seu_obj, bin_col, genes = NULL,
                                         bin_order = NULL,
                                         title = "Expression Heatmap",
                                         scale = "row") {
    expr_mat <- Seurat::GetAssayData(seu_obj, layer = "data")
    bins <- seu_obj@meta.data[[bin_col]]

    if (!is.null(bin_order)) {
        bins <- factor(bins, levels = bin_order)
    } else {
        bins <- factor(bins)
    }
    bin_levels <- levels(bins)

    if (!is.null(genes)) {
        genes <- intersect(genes, rownames(expr_mat))
        expr_mat <- expr_mat[genes, , drop = FALSE]
    }

    mean_mat <- sapply(bin_levels, function(b) {
        cells_in_bin <- bins == b
        if (sum(cells_in_bin) == 0) return(numeric(nrow(expr_mat)))
        Matrix::rowMeans(expr_mat[, cells_in_bin, drop = FALSE])
    })
    mean_mat <- as.matrix(mean_mat)
    colnames(mean_mat) <- bin_levels

    # Remove zero-variance genes
    keep_genes <- apply(mean_mat, 1, sd) > 0
    mean_mat <- mean_mat[keep_genes, , drop = FALSE]

    hm <- pheatmap::pheatmap(
        mean_mat, scale = scale,
        cluster_rows = TRUE, cluster_cols = FALSE,
        main = title,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
    )
    return(hm)
}


#' Generate a compiled multi-page PDF report of smiDE results
#'
#' @param plots       named list of ggplot/pheatmap objects
#' @param run_tag     character, run identifier for title page
#' @param config      list of config parameters for title page
#' @param output_file character, path to output PDF
#' @param width       numeric, page width (default 10)
#' @param height      numeric, page height (default 8)
#' @return invisible(output_file)
generate_smide_report <- function(plots, run_tag, config, output_file,
                                   width = 10, height = 8) {
    pdf(output_file, width = width, height = height)

    # Title page
    plot.new()
    text(0.5, 0.75, "smiDE Report", cex = 2.5, font = 2)
    text(0.5, 0.60, run_tag, cex = 1.4, font = 3)
    text(0.5, 0.45, paste("Date:", Sys.Date()), cex = 1.1)
    info_lines <- c(
        paste("Cell type:", config$celltype_oi),
        paste("DE variable:", config$de_var),
        paste("Study:", config$study_name),
        paste("Family:", config$family),
        paste("FDR cutoff:", config$fdr_cut)
    )
    for (i in seq_along(info_lines)) {
        text(0.5, 0.35 - (i - 1) * 0.06, info_lines[i], cex = 1.0)
    }

    # Print each plot
    for (nm in names(plots)) {
        p <- plots[[nm]]
        if (inherits(p, "gg") || inherits(p, "ggplot") ||
            inherits(p, "patchwork")) {
            print(p)
        } else if (inherits(p, "pheatmap")) {
            print(p)
        } else if (is.list(p)) {
            for (sub_p in p) {
                if (inherits(sub_p, "gg")) print(sub_p)
            }
        }
    }

    dev.off()
    message("Report saved to: ", output_file)
    invisible(output_file)
}
