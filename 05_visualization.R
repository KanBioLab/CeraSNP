# ==============================================================================
# Pipeline: CeraSNP - Plastome & rDNA Core SNP Discovery for Cerasus Cultivars
# Module:   05_visualization.R (Genome Summary & Publication Visualizations)
# Version:  1.0.0 (GitHub Release Edition)
# License:  MIT License
# ==============================================================================

cat("\n============================================================\n")
cat(" [MODULE 5] Generating Genome Summary & Publication Visualizations...\n")
cat("============================================================\n")
flush.console()

# ------------------------------------------------------------------------------
# 0. Environment Setup & Dependency Verification
# ------------------------------------------------------------------------------
if (!exists("outdir")) {
  outdir <- file.path(getwd(), "results")
}
if (!dir.exists(outdir)) {
  stop("Fatal Error: Results directory does not exist. Please run Modules 1-4 first.")
}

suppressPackageStartupMessages(library(pheatmap))

# ------------------------------------------------------------------------------
# 1. Full-Genome Polymorphism Density & Allele Spectrum (Table 11)
# ------------------------------------------------------------------------------
cat(" -> Synthesizing genomic compartment polymorphism density (Table 11)...\n")

if (exists("mat_p") && exists("mat_r")) {
  p_cand <- if (exists("poly_p")) poly_p else find_polymorphic_sites(mat_p, complete_rate = 1.00)
  r_cand <- if (exists("poly_r")) poly_r else find_polymorphic_sites(mat_r, complete_rate = 1.00)

  p_metrics <- lapply(p_cand, function(pos) get_site_genetic_metrics(mat_p[, pos]))
  r_metrics <- lapply(r_cand, function(pos) get_site_genetic_metrics(mat_r[, pos]))

  p_types <- vapply(p_metrics, function(x) x$Poly_Type, character(1))
  r_types <- vapply(r_metrics, function(x) x$Poly_Type, character(1))

  p_mafs  <- vapply(p_metrics, function(x) x$MAF, numeric(1))
  r_mafs  <- vapply(r_metrics, function(x) x$MAF, numeric(1))

  p_entropies <- vapply(p_metrics, function(x) x$Entropy, numeric(1))
  r_entropies <- vapply(r_metrics, function(x) x$Entropy, numeric(1))

  len_p <- ncol(mat_p)
  len_r <- ncol(mat_r)

  table_11 <- data.frame(
    Genomic_Compartment         = c("Plastome", "Nuclear_45S_rDNA", "Combined_Total"),
    Aligned_Length_bp           = c(len_p, len_r, len_p + len_r),
    Candidate_Polymorphic_SNPs  = c(length(p_cand), length(r_cand), length(p_cand) + length(r_cand)),
    Polymorphism_Density_per_kb = c(
      round(length(p_cand) / (len_p / 1000), 3),
      round(length(r_cand) / (len_r / 1000), 3),
      round((length(p_cand) + length(r_cand)) / ((len_p + len_r) / 1000), 3)
    ),
    Biallelic_SNPs_Count        = c(sum(p_types == "Biallelic"), sum(r_types == "Biallelic"), sum(c(p_types, r_types) == "Biallelic")),
    Biallelic_Ratio_Pct         = c(
      round(mean(p_types == "Biallelic") * 100, 2),
      round(mean(r_types == "Biallelic") * 100, 2),
      round(mean(c(p_types, r_types) == "Biallelic") * 100, 2)
    ),
    Multiallelic_SNPs_Count     = c(sum(p_types == "Multiallelic"), sum(r_types == "Multiallelic"), sum(c(p_types, r_types) == "Multiallelic")),
    Multiallelic_Ratio_Pct      = c(
      round(mean(p_types == "Multiallelic") * 100, 2),
      round(mean(r_types == "Multiallelic") * 100, 2),
      round(mean(c(p_types, r_types) == "Multiallelic") * 100, 2)
    ),
    Mean_Candidate_MAF          = c(round(mean(p_mafs), 4), round(mean(r_mafs), 4), round(mean(c(p_mafs, r_mafs)), 4)),
    Mean_Candidate_Entropy      = c(round(mean(p_entropies), 4), round(mean(r_entropies), 4), round(mean(c(p_entropies, r_entropies)), 4)),
    Core_Markers_Selected       = c(
      if (exists("sel_p")) length(sel_p) else NA_integer_,
      if (exists("sel_r")) length(sel_r) else NA_integer_,
      if (exists("sel_p") && exists("sel_r")) length(sel_p) + length(sel_r) else NA_integer_
    ),
    stringsAsFactors            = FALSE
  )

  write.csv(table_11, file.path(outdir, "11_Genome_Polymorphism_Summary.csv"), row.names = FALSE)
  cat(" -> [Table 11 Exported] 11_Genome_Polymorphism_Summary.csv\n")
}

# ------------------------------------------------------------------------------
# 2. Figure 1: Two-Tier Convergence Trajectories (Vector PDF)
# ------------------------------------------------------------------------------
cat(" -> Generating Figure 1: Two-Tier Convergence Trajectories...\n")

t1_file <- file.path(outdir, "06_Tier1_Greedy_Trajectory.csv")
t2_file <- file.path(outdir, "07_Tier2_Residual_Refinement_Trajectory.csv")

if (file.exists(t1_file) && file.exists(t2_file)) {
  h1 <- read.csv(t1_file, stringsAsFactors = FALSE)
  h2 <- read.csv(t2_file, stringsAsFactors = FALSE)

  pdf_traj_path <- file.path(outdir, "Figure_1_Two_Tier_Convergence_Trajectory.pdf")
  pdf(file = pdf_traj_path, width = 11, height = 5, bg = "white")

  par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.5, 1.5))

  # Panel A: Tier 1 Haplotype Saturation
  y1_col <- if ("Unique_Profiles" %in% names(h1)) "Unique_Profiles" else names(h1)[grep("unique|haplo", names(h1), ignore.case = TRUE)[1]]
  plot(
    h1$Step, h1[[y1_col]],
    type = "b", pch = 19, col = "#1F77B4", lwd = 2.2, cex = 1.1,
    xlab = "Greedy Selection Step (Tier 1 Plastome)",
    ylab = "Cumulative Unique Haplotypes",
    main = "A: Tier 1 Maternal Backbone Saturation",
    ylim = c(0, max(h1[[y1_col]]) * 1.15),
    xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = seq(1, nrow(h1), by = max(1, round(nrow(h1) / 10))))
  axis(2, las = 1)
  grid(nx = NA, ny = NULL, col = "#EBEBEB", lty = "solid")
  points(nrow(h1), tail(h1[[y1_col]], 1), pch = 17, col = "#D62728", cex = 1.8)
  text(nrow(h1) * 0.75, tail(h1[[y1_col]], 1) * 0.92, 
       labels = paste0("Saturated (Haps = ", tail(h1[[y1_col]], 1), ")"), 
       col = "#D62728", font = 2, cex = 0.85)

  # Panel B: Tier 2 Pairwise Collision Minimization
  y2_col <- if ("Collision_Pairs" %in% names(h2)) "Collision_Pairs" else names(h2)[grep("pair|collision", names(h2), ignore.case = TRUE)[1]]
  pair_vals <- h2[[y2_col]]

  plot(
    h2$Step, pair_vals,
    type = "b", pch = 19, col = "#FF7F0E", lwd = 2.2, cex = 1.1,
    xlab = "Greedy Selection Step (Tier 2 Nuclear rDNA)",
    ylab = "Global Pairwise Collisions C(S)",
    main = "B: Tier 2 Residual Collision Minimization",
    ylim = c(0, max(pair_vals) * 1.2),
    xaxt = "n", yaxt = "n", bty = "l"
  )
  axis(1, at = seq(1, nrow(h2), by = max(1, round(nrow(h2) / 8))))
  axis(2, las = 1)
  grid(nx = NA, ny = NULL, col = "#EBEBEB", lty = "solid")
  points(nrow(h2), tail(pair_vals, 1), pch = 17, col = "#D62728", cex = 1.8)
  text(nrow(h2) * 0.8, tail(pair_vals, 1) + (max(pair_vals) * 0.1), 
       labels = paste0("Locked (Pairs = ", tail(pair_vals, 1), ")"), 
       col = "#D62728", font = 2, cex = 0.85)

  invisible(dev.off())
  cat(" -> [Figure 1 Exported] Figure_1_Two_Tier_Convergence_Trajectory.pdf\n")
}

# ------------------------------------------------------------------------------
# 3. Figure 2 & Table 08: Pairwise Hamming Distance Matrix & Heatmap
# ------------------------------------------------------------------------------
cat(" -> Deriving Pairwise Hamming Distance Matrix & Heatmap (Figure 2)...\n")

db_file <- file.path(outdir, "01_Core_Molecular_ID_Database.csv")
if (file.exists(db_file)) {
  db_df <- read.csv(db_file, stringsAsFactors = FALSE)
  sample_names <- db_df$Sample_ID
  profiles     <- db_df$Combined_Barcode

  n_samples <- length(profiles)
  dist_matrix <- matrix(0L, nrow = n_samples, ncol = n_samples, dimnames = list(sample_names, sample_names))

  for (i in seq_len(n_samples)) {
    for (j in i:n_samples) {
      d <- calculate_hamming_dist(profiles[i], profiles[j])
      dist_matrix[i, j] <- d
      dist_matrix[j, i] <- d
    }
  }

  dist_out_df <- data.frame(Sample_ID = sample_names, dist_matrix, check.names = FALSE)
  write.csv(dist_out_df, file.path(outdir, "08_Pairwise_Hamming_Distance_Matrix.csv"), row.names = FALSE)
  cat(" -> [Table 08 Exported] 08_Pairwise_Hamming_Distance_Matrix.csv\n")

  pdf_heatmap_path <- file.path(outdir, "Figure_2_Hamming_Distance_Heatmap.pdf")
  pdf(file = pdf_heatmap_path, width = 10, height = 9.5, bg = "white")

  # High-contrast publication palette (Light blue to navy)
  heatmap_palette <- colorRampPalette(c("#FFFFFF", "#C6DBEF", "#4292C6", "#084594"))(100)

  pheatmap::pheatmap(
    dist_matrix,
    color            = heatmap_palette,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    treeheight_row   = 35,
    treeheight_col   = 35,
    fontsize_row     = if (n_samples > 80) 5 else 6.5,
    fontsize_col     = if (n_samples > 80) 5 else 6.5,
    border_color     = NA,
    main             = "Pairwise Hamming Distance Matrix across Core SNP Space"
  )
  invisible(dev.off())
  cat(" -> [Figure 2 Exported] Figure_2_Hamming_Distance_Heatmap.pdf\n")
}

cat("\n============================================================\n")
cat(" [MODULE 5 COMPLETE] Visualizations and Summary Tables Ready.\n")
cat("============================================================\n")