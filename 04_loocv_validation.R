# ==============================================================================
# Pipeline: CeraSNP - Plastome & rDNA Core SNP Discovery for Cerasus Cultivars
# Module:   04_loocv_validation.R (Strict Zero-Leakage LOOCV & Marker Stability)
# Version:  1.0.0 (GitHub Release Edition)
# License:  MIT License
# ==============================================================================

cat("\n============================================================\n")
cat(" [MODULE 4] Executing Strict Zero-Leakage LOOCV & Stability Analysis...\n")
cat("============================================================\n")
flush.console()

# ------------------------------------------------------------------------------
# 0. Defensive Checks and Environment Assertions
# ------------------------------------------------------------------------------
if (!exists("mat_p") || !exists("mat_r") || !exists("common")) {
  stop("Fatal Error: Aligned sequence matrices (mat_p, mat_r) not found. Run Module 3 first.")
}

n <- length(common)
loocv        <- vector("list", n)
lo_p         <- vector("list", n)
lo_r         <- vector("list", n)
jaccard_all  <- rep(NA_real_, n)

# Baseline discovery set identifiers with cross-genome namespace tags
disc_tagged <- tag_sites(sel_p, sel_r)

cat(sprintf(" -> Launching %d independent LOOCV discovery folds (Strict Isolation)...\n", n))
flush.console()

t_start <- Sys.time()

# ------------------------------------------------------------------------------
# 1. Zero-Leakage Cross-Validation Loop
# ------------------------------------------------------------------------------
for (i in seq_len(n)) {
  tr <- setdiff(seq_len(n), i)
  te <- i

  tp  <- mat_p[tr, , drop = FALSE]
  trd <- mat_r[tr, , drop = FALSE]
  ep  <- mat_p[te, , drop = FALSE]
  er  <- mat_r[te, , drop = FALSE]

  # Independent polymorphism screening strictly within current training fold
  cp <- find_polymorphic_sites(tp, plastome_complete_rate)
  cr <- find_polymorphic_sites(trd, rdna_complete_rate)

  # Tier 1 selection in training fold
  ps     <- integer(0)
  tpp    <- rep("", length(tr))
  p_stop <- "No_Plastome_Candidates"
  if (length(cp) > 0) {
    f1     <- adaptive_greedy_tier1(tp, cp, max_p_trials, verbose = FALSE)
    ps     <- f1$sites
    tpp    <- make_profile(tp, ps)
    p_stop <- f1$stop_reason
  }

  # Tier 2 selection conditioned on training plastome background
  rs     <- integer(0)
  trp    <- rep("", length(tr))
  r_stop <- "No_rDNA_Candidates"
  if (length(cr) > 0) {
    f2     <- adaptive_greedy_tier2(trd, cr, tpp, max_r_trials, verbose = FALSE)
    rs     <- f2$sites
    trp    <- make_profile(trd, rs)
    r_stop <- f2$stop_reason
  }

  lo_p[[i]] <- ps
  lo_r[[i]] <- rs

  # Evaluate locus set stability against baseline discovery set
  fold_tagged <- tag_sites(ps, rs)
  jaccard_all[i] <- calc_jaccard(disc_tagged, fold_tagged)

  if (length(ps) == 0 && length(rs) == 0) {
    loocv[[i]] <- data.frame(
      Held_Out_Sample          = common[i],
      Plastome_Markers         = 0L,
      rDNA_Markers             = 0L,
      Total_Markers            = 0L,
      Training_Unique_Barcodes = 0L,
      Barcode_Collision        = NA,
      Nearest_SNP_Hamming      = NA_real_,
      Unseen_Alleles_Count     = NA_integer_,
      Jaccard_Similarity_Pct   = NA_real_,
      Tier1_Stop_Reason        = p_stop,
      Tier2_Stop_Reason        = r_stop,
      stringsAsFactors         = FALSE
    )
    next
  }

  train_bar <- paste(tpp, trp, sep = "_rDNA_")

  # Project held-out blind sample using loci selected purely from training fold
  epp <- make_profile(ep, ps)
  erp <- make_profile(er, rs)
  test_bar <- paste(epp, erp, sep = "_rDNA_")

  # Evaluate collision burden
  collision <- test_bar %in% train_bar

  # Calculate nearest Hamming distance to training barcode space
  test_str   <- paste0(epp, erp)
  train_strs <- paste0(tpp, trp)
  d <- vapply(train_strs, function(ts) calculate_hamming_dist(test_str, ts), integer(1))
  nearest_d <- if (all(is.na(d))) NA_real_ else min(d, na.rm = TRUE)

  # Monitor novel unseen alleles in blind sample
  unseen_cnt <- 0L
  all_sites <- c(ps, rs)
  all_train_mat <- cbind(tp[, ps, drop = FALSE], trd[, rs, drop = FALSE])
  all_test_vec  <- c(ep[, ps, drop = TRUE], er[, rs, drop = TRUE])

  if (length(all_sites) > 0) {
    for (s_idx in seq_along(all_sites)) {
      train_alleles <- unique(all_train_mat[, s_idx])
      if (!all_test_vec[s_idx] %in% train_alleles) {
        unseen_cnt <- unseen_cnt + 1L
      }
    }
  }

  loocv[[i]] <- data.frame(
    Held_Out_Sample          = common[i],
    Plastome_Markers         = length(ps),
    rDNA_Markers             = length(rs),
    Total_Markers            = length(ps) + length(rs),
    Training_Unique_Barcodes = length(unique(train_bar)),
    Barcode_Collision        = collision,
    Nearest_SNP_Hamming      = nearest_d,
    Unseen_Alleles_Count     = unseen_cnt,
    Jaccard_Similarity_Pct   = round(jaccard_all[i] * 100, 2),
    Tier1_Stop_Reason        = p_stop,
    Tier2_Stop_Reason        = r_stop,
    stringsAsFactors         = FALSE
  )

  # Console progress monitor
  pct <- round(i / n * 100, 1)
  current_collisions <- sum(vapply(loocv[1:i], function(x) x$Barcode_Collision, logical(1)), na.rm = TRUE)
  cat(sprintf("\r -> [LOOCV Progress] Fold: %2d/%2d (%5.1f%%) | Sample: %-15s | Collisions: %d", 
              i, n, pct, substr(common[i], 1, 15), current_collisions))
  flush.console()
}

cat("\n -> LOOCV execution completed in ", round(difftime(Sys.time(), t_start, units = "secs"), 1), " seconds.\n\n")
flush.console()

loocv_df <- do.call(rbind, loocv)
write.csv(loocv_df, file.path(outdir, "04_Marker_Selection_LOOCV.csv"), row.names = FALSE)
cat(" -> [Table 04 Exported] 04_Marker_Selection_LOOCV.csv\n")

# ------------------------------------------------------------------------------
# 2. Marker Selection Retention Frequency & Jaccard Stability
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(" [STAGE 6/7] Evaluating Marker Retention & Set Stability Indices...\n")
cat("============================================================\n")
flush.console()

calculate_retention <- function(discovery, folds, prefix, type, nfolds) {
  if (!length(discovery)) return(data.frame())
  f <- vapply(discovery, function(s) {
    sum(vapply(folds, function(x) s %in% x, logical(1)))
  }, integer(1))

  data.frame(
    Marker_ID                = sprintf("%s_%02d", prefix, seq_along(discovery)),
    Genome_Type              = type,
    Align_Pos                = discovery,
    LOOCV_Retention_Freq_Pct = round(f / nfolds * 100, 2),
    Retention_Classification = ifelse(f / nfolds >= 0.90, "High_Retention (>=90%)", "Lower_Retention (<90%)"),
    stringsAsFactors         = FALSE
  )
}

stab_p <- calculate_retention(sel_p, lo_p, "Plastome_SNP", "Plastome", n)
stab_r <- calculate_retention(sel_r, lo_r, "rDNA_SNP", "45S_rDNA", n)
stab   <- rbind(stab_p, stab_r)

mean_retention   <- mean(stab$LOOCV_Retention_Freq_Pct)
mean_jaccard_all <- mean(jaccard_all, na.rm = TRUE) * 100
high_ret_count   <- sum(stab$LOOCV_Retention_Freq_Pct >= 90)

cat(sprintf(" -> Mean marker retention across LOOCV folds: %.2f%%\n", mean_retention))
cat(sprintf(" -> Mean Jaccard locus set stability: %.2f%%\n", mean_jaccard_all))
cat(sprintf(" -> High-retention core markers (>= 90%%): %d / %d\n\n", high_ret_count, nrow(stab)))
flush.console()

write.csv(stab, file.path(outdir, "05_Marker_Selection_Stability.csv"), row.names = FALSE)
cat(" -> [Table 05 Exported] 05_Marker_Selection_Stability.csv\n")

# ------------------------------------------------------------------------------
# 3. Independent External Validation (If Specified)
# ------------------------------------------------------------------------------
ext_df <- data.frame()
if (!is.null(external_plastome_file) && !is.null(external_rdna_file)) {
  cat("============================================================\n")
  cat(" [STAGE 7/7] Evaluating Independent Held-Out Validation Panel...\n")
  cat("============================================================\n")
  flush.console()

  xp <- read_clean_fasta(external_plastome_file)
  xr <- read_clean_fasta(external_rdna_file)
  ec <- intersect(names(xp), names(xr))
  if (!length(ec)) stop("Fatal Error: No common taxa in external validation panel.")
  
  leak_samples <- intersect(ec, common)
  if (length(leak_samples) > 0) {
    stop("Fatal Leakage Error: External validation set contains discovery taxa: ", paste(leak_samples, collapse = ", "))
  }

  xp <- xp[match(ec, names(xp))]
  xr <- xr[match(ec, names(xr))]
  check_fixed_alignment(xp, "External plastome FASTA")
  check_fixed_alignment(xr, "External 45S rDNA FASTA")

  xpm <- seq_to_matrix(xp)
  xrm <- seq_to_matrix(xr)

  if (ncol(xpm) != ncol(mat_p)) stop("Coordinate Error: External plastome alignment width mismatch.")
  if (ncol(xrm) != ncol(mat_r)) stop("Coordinate Error: External rDNA alignment width mismatch.")

  x_p_profile <- make_profile(xpm, sel_p)
  x_r_profile <- make_profile(xrm, sel_r)
  xb <- paste(x_p_profile, x_r_profile, sep = "_rDNA_")

  ext_df <- do.call(rbind, lapply(seq_along(ec), function(j) {
    match_count <- sum(combined == xb[j])
    concordance_status <- if (match_count == 1) "Unique_Barcode_Concordance"
                          else if (match_count > 1) "Multiple_Barcode_Concordance"
                          else "Novel_Unseen_Barcode"
    data.frame(
      External_Sample    = ec[j],
      Assigned_Barcode   = xb[j],
      Discovery_Matches  = match_count,
      Concordance_Status = concordance_status,
      stringsAsFactors   = FALSE
    )
  }))
  write.csv(ext_df, file.path(outdir, "08_Independent_Heldout_Validation.csv"), row.names = FALSE)
  cat(sprintf(" -> Evaluated %d external blind test samples.\n\n", nrow(ext_df)))
}

# ------------------------------------------------------------------------------
# 4. Executive Summary Table (Table 00)
# ------------------------------------------------------------------------------
resolution <- resolved / n * 100
valid_loocv <- !is.na(loocv_df$Barcode_Collision)
loocv_non_collision <- if (any(valid_loocv)) mean(!loocv_df$Barcode_Collision[valid_loocv]) * 100 else NA_real_
median_dist <- if (any(!is.na(loocv_df$Nearest_SNP_Hamming))) median(loocv_df$Nearest_SNP_Hamming, na.rm = TRUE) else NA_real_
total_unseen_alleles <- sum(loocv_df$Unseen_Alleles_Count, na.rm = TRUE)

summary_exec <- data.frame(
  Discovery_Samples                   = n,
  Plastome_Candidate_SNPs             = length(poly_p),
  rDNA_Candidate_SNPs                 = length(poly_r),
  Plastome_Core_SNPs                  = length(sel_p),
  rDNA_Core_SNPs                      = length(sel_r),
  Total_Core_Markers                  = length(sel_p) + length(sel_r),
  Unique_Resolved_Samples             = resolved,
  Discovery_Resolution_Rate_Pct       = round(resolution, 2),
  Tier1_Discovery_Stop_Reason         = t1_res$stop_reason,
  Tier2_Discovery_Stop_Reason         = t2_res$stop_reason,
  LOOCV_Valid_Folds                   = sum(valid_loocv),
  LOOCV_NonCollision_Rate_Pct         = round(loocv_non_collision, 2),
  Mean_Marker_Retention_Frequency_Pct = round(mean_retention, 2),
  Mean_LOOCV_Jaccard_Stability_Pct    = round(mean_jaccard_all, 2),
  Median_Nearest_SNP_Hamming          = round(median_dist, 2),
  Total_Unseen_Alleles_Encountered    = total_unseen_alleles,
  Discovery_Collision_Pairs           = collision_pairs(combined),
  Discovery_Collision_Samples         = collision_samples(combined),
  stringsAsFactors                    = FALSE
)

write.csv(summary_exec, file.path(outdir, "00_Executive_Summary.csv"), row.names = FALSE)
cat(" -> [Table 00 Exported] 00_Executive_Summary.csv\n")
cat(" -> [MODULE 4 COMPLETE] Cross-validation and stability assessment finalized.\n\n")