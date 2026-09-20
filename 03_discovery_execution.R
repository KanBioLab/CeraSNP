# ==============================================================================
# Pipeline: CeraSNP - Plastome & rDNA Core SNP Discovery for Cerasus Cultivars
# Module:   03_discovery_execution.R (Discovery Execution & Unified Database Export)
# Version:  1.0.0 (Unified Production Edition for GitHub)
# License:  MIT License
# ==============================================================================

cat("============================================================\n")
cat(" [STAGE 1/7] Loading FASTA Alignments & Harmonizing Sample Cohorts...\n")
cat("============================================================\n")
flush.console()

# ------------------------------------------------------------------------------
# 1. Load Sequences and Align Matrices
# ------------------------------------------------------------------------------
seqs_p <- read_clean_fasta(plastome_file)
seqs_r <- read_clean_fasta(rdna_file)

raw_p_names <- names(seqs_p)
raw_r_names <- names(seqs_r)
common <- intersect(raw_p_names, raw_r_names)

if (length(common) < 2) {
  stop("Fatal Error: Shared taxa between plastome and 45S rDNA alignments are fewer than 2.")
}

p_only <- setdiff(raw_p_names, common)
r_only <- setdiff(raw_r_names, common)

cat(sprintf(" -> Synchronized %d shared cultivars across both compartments.\n", length(common)))
if (length(p_only) > 0) cat(sprintf(" -> Audit Note: %d taxa present only in plastome alignment omitted.\n", length(p_only)))
if (length(r_only) > 0) cat(sprintf(" -> Audit Note: %d taxa present only in 45S rDNA alignment omitted.\n", length(r_only)))

seqs_p <- seqs_p[match(common, names(seqs_p))]
seqs_r <- seqs_r[match(common, names(seqs_r))]

check_fixed_alignment(seqs_p, "Plastome Alignment")
check_fixed_alignment(seqs_r, "45S rDNA Alignment")

mat_p <- seq_to_matrix(seqs_p)
mat_r <- seq_to_matrix(seqs_r)

sample_names <- common
n_samples    <- length(common)
seq_len_p    <- ncol(mat_p)
seq_len_r    <- ncol(mat_r)

cat(sprintf(" -> Aligned Plastome Width: %d bp | 45S rDNA Width: %d bp\n\n", seq_len_p, seq_len_r))
flush.console()

# ------------------------------------------------------------------------------
# 2. Polymorphic Candidate Pool Screening (100% Complete Calls)
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(" [STAGE 2/7] Screening Zero-Missing Candidate SNP Pools...\n")
cat("============================================================\n")
flush.console()

poly_p <- find_polymorphic_sites(mat_p, plastome_complete_rate)
poly_r <- find_polymorphic_sites(mat_r, rdna_complete_rate)

if (!length(poly_p)) stop("Fatal Error: No qualified polymorphic sites found in plastome alignment.")
if (!length(poly_r)) stop("Fatal Error: No qualified polymorphic sites found in 45S rDNA alignment.")

cat(sprintf(" -> Plastome qualified candidate SNPs (100%% complete): %d loci\n", length(poly_p)))
cat(sprintf(" -> 45S rDNA qualified candidate SNPs (100%% complete): %d loci\n\n", length(poly_r)))
flush.console()

# ------------------------------------------------------------------------------
# 3. Two-Tier Cascaded Adaptive Greedy Discovery
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(" [STAGE 3/7] Executing Two-Tier Cascaded Adaptive Greedy Search...\n")
cat("============================================================\n")

# --- Tier 1: Plastome Maternal Backbone Optimization ---
cat(">>> Launching Tier 1: Maternal Plastome Backbone Forward Selection...\n")
flush.console()

t1_start <- Sys.time()
t1_res   <- adaptive_greedy_tier1(mat_p, poly_p, max_p_trials, verbose = TRUE)
t1_end   <- Sys.time()

sel_p       <- t1_res$sites
p_prof      <- make_profile(mat_p, sel_p)
n_core_p    <- length(sel_p)
unique_haps <- length(unique(p_prof))

# 构建按频率排序的标准母系单倍型编号 (Hap_01 ~ Hap_xx)
p_counts     <- table(p_prof)
haplo_levels <- names(sort(p_counts, decreasing = TRUE))
hap_map      <- setNames(sprintf("Hap_%02d", seq_along(haplo_levels)), haplo_levels)
hap_codes    <- unname(hap_map[p_prof])

cat(sprintf("\n -> Tier 1 Complete in %.2f s: %d core plastome SNPs converged into %d maternal haplotypes.\n\n", 
            as.numeric(difftime(t1_end, t1_start, units = "secs")), n_core_p, unique_haps))

# --- Tier 2: Nuclear rDNA Residual Collision Decoupling ---
cat(">>> Launching Tier 2: Nuclear 45S rDNA Residual Collision Minimization...\n")
flush.console()

t2_start <- Sys.time()
t2_res   <- adaptive_greedy_tier2(mat_r, poly_r, p_prof, max_r_trials, verbose = TRUE)
t2_end   <- Sys.time()

sel_r    <- t2_res$sites
n_core_r <- length(sel_r)
r_prof   <- make_profile(mat_r, sel_r)
combined <- paste(p_prof, r_prof, sep = "_rDNA_")

counts   <- table(combined)
is_uniquely_resolved <- as.logical(counts[combined] == 1)
resolved <- sum(is_uniquely_resolved)
res_rate <- round((resolved / n_samples) * 100, 2)

cat(sprintf("\n -> Tier 2 Complete in %.2f s: %d core rDNA SNPs converged.\n", 
            as.numeric(difftime(t2_end, t2_start, units = "secs")), n_core_r))
cat(sprintf(" -> Total Core Markers: %d (Plastome %d + rDNA %d)\n", n_core_p + n_core_r, n_core_p, n_core_r))
cat(sprintf(" -> Final Identification: %d / %d cultivars uniquely discriminated (%.2f%%).\n\n", 
            resolved, n_samples, res_rate))
flush.console()

# ------------------------------------------------------------------------------
# 4. In-Silico KASP Flanking QC & Compatibility
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(" [STAGE 4/7] In-Silico KASP Flanking QC & Assay Design Screening...\n")
cat("============================================================\n")
flush.console()

kp <- generate_kasp_qc(mat_p, sel_p, poly_p, "Plastome_SNP", "Plastome", flank_len)
kr <- generate_kasp_qc(mat_r, sel_r, poly_r, "rDNA_SNP", "45S_rDNA", flank_len)
kasp_qc_df <- rbind(kp, kr)

cat(sprintf(" -> Evaluated %d Core Markers for KASP assay compatibility.\n", nrow(kasp_qc_df)))
cat(sprintf(" -> Mean Flanking GC Content: %.2f%%\n", mean(kasp_qc_df$Flank_GC_Pct, na.rm = TRUE)))
cat(sprintf(" -> Markers with Secondary Flanking Variations: %d\n", sum(kasp_qc_df$Flank_Secondary_SNP_Warn == "WARN_SECONDARY_POLY")))
cat(sprintf(" -> Markers with Homopolymer Runs (>=5 bp): %d\n\n", sum(kasp_qc_df$Homopolymer_Warn == "WARN_POLY_RUN")))
flush.console()

# ------------------------------------------------------------------------------
# 5. Database Synthesis & Standardized File Export
# ------------------------------------------------------------------------------
cat(" -> Exporting Core Discovery Databases to: ", outdir, "\n")

# 5.1 Table 01: Core Molecular ID Database (展开具体 Pos 列)
barcode_p_mat <- mat_p[, sel_p, drop = FALSE]
colnames(barcode_p_mat) <- sprintf("P_Pos_%d", sel_p)

barcode_r_mat <- mat_r[, sel_r, drop = FALSE]
colnames(barcode_r_mat) <- sprintf("R_Pos_%d", sel_r)

db_df <- data.frame(
  Sample_ID             = sample_names,
  Maternal_Haplotype    = hap_codes,
  Is_Unique_Resolved    = ifelse(is_uniquely_resolved, "YES", "NO"),
  Plastome_Barcode      = p_prof,
  rDNA_Barcode          = r_prof,
  Combined_Barcode      = combined,
  as.data.frame(barcode_p_mat),
  as.data.frame(barcode_r_mat),
  stringsAsFactors      = FALSE,
  check.names           = FALSE
)
write.csv(db_df, file.path(outdir, "01_Core_Molecular_ID_Database.csv"), row.names = FALSE)
cat(" -> [Table 01 Exported] 01_Core_Molecular_ID_Database.csv\n")

# 5.2 Table 07: Maternal Haplotype Diversity Summary (Haplotype-level)
haplo_stat_list <- lapply(haplo_levels, function(h_prof) {
  sub_df <- db_df[db_df$Plastome_Barcode == h_prof, ]
  h_name <- sub_df$Maternal_Haplotype[1]
  cnt    <- nrow(sub_df)
  data.frame(
    Haplotype_ID      = h_name,
    Barcode_Length    = n_core_p,
    Plastome_Barcode  = h_prof,
    Sample_Count      = cnt,
    Frequency_Pct     = round((cnt / n_samples) * 100, 2),
    Is_Unique         = ifelse(cnt == 1, "YES", "NO"),
    Sample_Members    = paste(sub_df$Sample_ID, collapse = ", "),
    stringsAsFactors  = FALSE
  )
})
haplo_summary_df <- do.call(rbind, haplo_stat_list)
write.csv(haplo_summary_df, file.path(outdir, "07_Maternal_Haplotype_Diversity_Summary.csv"), row.names = FALSE)
cat(sprintf(" -> [Table 07 Exported] 07_Maternal_Haplotype_Diversity_Summary.csv (Total: %d Haplotypes: %d Shared + %d Unique)\n",
            nrow(haplo_summary_df), sum(haplo_summary_df$Is_Unique == "NO"), sum(haplo_summary_df$Is_Unique == "YES")))

# 5.3 Table 02: Residual Unresolved Clusters (Shared Barcodes after Tier 2)
unres <- db_df[db_df$Is_Unique_Resolved == "NO", , drop = FALSE]
if (nrow(unres)) {
  sp <- split(unres$Sample_ID, unres$Combined_Barcode)
  unresolved_agg <- data.frame(
    Combined_Barcode = names(sp),
    Cluster_Size     = vapply(sp, length, integer(1)),
    Cluster_Members  = vapply(sp, paste, collapse = ", ", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  )
  unresolved_agg <- unresolved_agg[order(unresolved_agg$Cluster_Size, decreasing = TRUE), ]
} else {
  unresolved_agg <- data.frame(Combined_Barcode = character(), Cluster_Size = integer(), Cluster_Members = character(), stringsAsFactors = FALSE)
}
write.csv(unresolved_agg, file.path(outdir, "02_Unresolved_Clusters.csv"), row.names = FALSE)
cat(" -> [Table 02 Exported] 02_Unresolved_Clusters.csv\n")

# 5.4 Table 03: KASP QC
write.csv(kasp_qc_df, file.path(outdir, "03_KASP_Marker_QC.csv"), row.names = FALSE)
cat(" -> [Table 03 Exported] 03_KASP_Marker_QC.csv\n")

# 5.5 Trajectory Tables
write.csv(t1_res$history, file.path(outdir, "06_Tier1_Greedy_Trajectory.csv"), row.names = FALSE)
write.csv(t2_res$history, file.path(outdir, "07_Tier2_Residual_Refinement_Trajectory.csv"), row.names = FALSE)
cat(" -> [Trajectories Exported] CSV 06 & 07 Trajectory Tables.\n")

# 5.6 Tables 09 & 10: Full Candidate Pools with Detailed Genetic Metrics
p_metrics_list <- lapply(poly_p, function(p) get_site_genetic_metrics(mat_p[, p]))
cand_p_df <- data.frame(
  Align_Pos    = poly_p,
  Poly_Type    = vapply(p_metrics_list, function(x) x$Poly_Type, character(1)),
  Alleles      = vapply(p_metrics_list, function(x) x$Alleles, character(1)),
  Major_Allele = vapply(p_metrics_list, function(x) x$Major_Allele, character(1)),
  Major_Freq   = vapply(p_metrics_list, function(x) x$Major_Freq, numeric(1)),
  MAF          = vapply(p_metrics_list, function(x) x$MAF, numeric(1)),
  Entropy      = vapply(p_metrics_list, function(x) x$Entropy, numeric(1)),
  Is_Core_SNP  = ifelse(poly_p %in% sel_p, "YES", "NO"),
  stringsAsFactors = FALSE
)
write.csv(cand_p_df, file.path(outdir, "09_Discovery_Plastome_Candidate_SNPs.csv"), row.names = FALSE)

r_metrics_list <- lapply(poly_r, function(p) get_site_genetic_metrics(mat_r[, p]))
cand_r_df <- data.frame(
  Align_Pos    = poly_r,
  Poly_Type    = vapply(r_metrics_list, function(x) x$Poly_Type, character(1)),
  Alleles      = vapply(r_metrics_list, function(x) x$Alleles, character(1)),
  Major_Allele = vapply(r_metrics_list, function(x) x$Major_Allele, character(1)),
  Major_Freq   = vapply(r_metrics_list, function(x) x$Major_Freq, numeric(1)),
  MAF          = vapply(r_metrics_list, function(x) x$MAF, numeric(1)),
  Entropy      = vapply(r_metrics_list, function(x) x$Entropy, numeric(1)),
  Is_Core_SNP  = ifelse(poly_r %in% sel_r, "YES", "NO"),
  stringsAsFactors = FALSE
)
write.csv(cand_r_df, file.path(outdir, "10_Discovery_rDNA_Candidate_SNPs.csv"), row.names = FALSE)
cat(" -> [Tables 09 & 10 Exported] Detailed Candidate Pool Matrices.\n")

# 5.7 Table 00: Unified Executive Summary (Two-column Key-Value format)
exec_summary <- data.frame(
  Metric = c(
    "Discovery_Cohort_Samples",
    "Plastome_Aligned_Length_bp",
    "rDNA_Aligned_Length_bp",
    "Plastome_Candidate_ZeroMissing_SNPs",
    "rDNA_Candidate_ZeroMissing_SNPs",
    "Plastome_Core_SNPs_Tier1",
    "rDNA_Core_SNPs_Tier2",
    "Total_Core_Markers",
    "Maternal_Haplotypes_Tier1",
    "Unique_Plastome_Haplotypes",
    "Shared_Plastome_Haplotypes",
    "Uniquely_Resolved_Cultivars_Tier2",
    "Final_Identification_Rate_Pct",
    "Tier1_Greedy_Stop_Reason",
    "Tier2_Greedy_Stop_Reason",
    "Tier1_Execution_Time_Secs",
    "Tier2_Execution_Time_Secs"
  ),
  Value = c(
    as.character(n_samples),
    as.character(seq_len_p),
    as.character(seq_len_r),
    as.character(length(poly_p)),
    as.character(length(poly_r)),
    as.character(n_core_p),
    as.character(n_core_r),
    as.character(n_core_p + n_core_r),
    as.character(unique_haps),
    as.character(sum(haplo_summary_df$Is_Unique == "YES")),
    as.character(sum(haplo_summary_df$Is_Unique == "NO")),
    as.character(resolved),
    as.character(res_rate),
    t1_res$stop_reason,
    t2_res$stop_reason,
    as.character(round(as.numeric(difftime(t1_end, t1_start, units = "secs")), 2)),
    as.character(round(as.numeric(difftime(t2_end, t2_start, units = "secs")), 2))
  ),
  stringsAsFactors = FALSE
)
write.csv(exec_summary, file.path(outdir, "00_Executive_Summary.csv"), row.names = FALSE)
cat(" -> [Table 00 Exported] 00_Executive_Summary.csv\n")
cat(" -> [MODULE 3 COMPLETE] Two-tier discovery and comprehensive databases exported.\n\n")