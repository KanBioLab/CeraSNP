# ==============================================================================
# Pipeline: CeraSNP - Plastome & rDNA Core SNP Discovery for Cerasus Cultivars
# Module:   01_config_and_utils.R (Environment Setup & Core Utilities)
# Version:  1.0.0 (GitHub Release Edition)
# License:  MIT License
# ==============================================================================

cat("\n============================================================\n")
cat(" [STAGE 0/5] Initializing CeraSNP Dependencies & Utilities...\n")
cat("============================================================\n")
flush.console()

# ------------------------------------------------------------------------------
# 0. Global Configuration and Parameters
# ------------------------------------------------------------------------------
work_dir <- getwd()
outdir   <- file.path(work_dir, "results")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Standard input FASTA alignment files (Place these in the working directory)
plastome_file <- file.path(work_dir, "plastome_aligned.fasta")
rdna_file     <- file.path(work_dir, "rdna_aligned.fasta")

# Independent validation panel (Leave as NULL if not applicable)
external_plastome_file <- NULL
external_rdna_file     <- NULL

# Computational iteration caps for greedy selection
max_p_trials <- 100
max_r_trials <- 50

# Flanking sequence physical half-window (+/- bp) for KASP probe QC
flank_len <- 50

# Strict quality control thresholds (1.00 = 100% complete canonical calls)
plastome_complete_rate <- 1.00
rdna_complete_rate     <- 1.00

# ------------------------------------------------------------------------------
# 1. Dependency Management
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

required_cran <- c("pheatmap")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

required_bioc <- c("Biostrings")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat(" -> Core dependencies (Biostrings, pheatmap) loaded successfully.\n")
cat(sprintf(" -> Designated output path: %s\n\n", normalizePath(outdir)))
flush.console()

# ------------------------------------------------------------------------------
# 2. Sequence Sanitization & Alignment Validation Functions
# ------------------------------------------------------------------------------

#' Read and Sanitize Multi-FASTA File
#' @param file Character path to input FASTA.
#' @return Sanitized DNAStringSet object.
#' @export
read_clean_fasta <- function(file) {
  if (!file.exists(file)) {
    stop("Fatal Error: FASTA file not found: ", file)
  }
  x <- Biostrings::readDNAStringSet(file)
  if (length(x) == 0) {
    stop("Fatal Error: FASTA file contains no sequences: ", file)
  }
  # Keep taxon ID before first whitespace
  nm <- vapply(strsplit(names(x), "\\s+"), function(z) z[1], character(1))
  names(x) <- nm
  
  if (any(!nzchar(names(x)))) stop("Fatal Error: Empty identifier detected in: ", file)
  if (anyDuplicated(names(x))) stop("Fatal Error: Duplicate identifier detected in: ", file)
  x
}

#' Assert Fixed Alignment Width
#' @param x DNAStringSet object.
#' @param label Character descriptive label.
#' @export
check_fixed_alignment <- function(x, label) {
  w <- Biostrings::width(x)
  if (length(unique(w)) != 1) {
    stop(sprintf("Fatal Error: %s sequences have variable lengths. Full alignment required.", label))
  }
  invisible(TRUE)
}

#' Convert DNAStringSet to 2D Character Matrix
#' @param x DNAStringSet object.
#' @return 2D character matrix (taxa x sites).
#' @export
seq_to_matrix <- function(x) {
  if (length(unique(Biostrings::width(x))) != 1) {
    stop("seq_to_matrix: Aligned sequences must have equal widths.")
  }
  z <- lapply(as.character(x), function(s) {
    strsplit(toupper(s), "", fixed = TRUE)[[1]]
  })
  m <- do.call(rbind, z)
  rownames(m) <- names(x)
  if (is.vector(m)) m <- matrix(m, nrow = 1, dimnames = list(names(x), NULL))
  m
}

#' Identify Zero-Missing Polymorphic Candidate Loci
#' @param mat 2D character matrix.
#' @param complete_rate Completeness threshold (default: 1.00).
#' @return Integer vector of valid polymorphic column indices.
#' @export
find_polymorphic_sites <- function(mat, complete_rate = 1.00) {
  n <- nrow(mat)
  min_valid <- ceiling(complete_rate * n)
  good <- apply(mat, 2, function(col) {
    v <- col[toupper(col) %in% c("A", "C", "G", "T")]
    (length(v) >= min_valid) && (length(unique(v)) >= 2)
  })
  which(good)
}

# ------------------------------------------------------------------------------
# 3. Barcode Assembly and Objective Function Metrics
# ------------------------------------------------------------------------------

#' Assemble Barcode Profile
#' @param mat 2D character matrix.
#' @param sites Integer vector of selected column indices.
#' @return Character vector of concatenated profiles.
#' @export
make_profile <- function(mat, sites) {
  if (!length(sites)) return(rep("", nrow(mat)))
  apply(mat[, sites, drop = FALSE], 1, paste0, collapse = "")
}

#' Calculate Scalar Hamming Distance
#' @param a Character barcode string.
#' @param b Character barcode string.
#' @return Integer distance or NA.
#' @export
calculate_hamming_dist <- function(a, b) {
  if (length(a) != 1 || length(b) != 1 || is.na(a) || is.na(b)) return(NA_integer_)
  if (nchar(a) != nchar(b) || nchar(a) == 0) return(NA_integer_)
  sum(utf8ToInt(a) != utf8ToInt(b))
}

#' Compute Global Pairwise Collision Burden
#' Objective function for Tier 2: C(S) = sum [n_k * (n_k - 1) / 2]
#' @param profile Character vector of barcodes.
#' @return Numeric count of pairwise collisions.
#' @export
collision_pairs <- function(profile) {
  tb <- table(profile)
  sum(tb * (tb - 1) / 2)
}

#' Count Unresolved Samples in Colliding Clusters
#' @param profile Character vector of barcodes.
#' @return Integer count of unresolved samples.
#' @export
collision_samples <- function(profile) {
  tb <- table(profile)
  if (!length(tb)) return(0L)
  sum(tb[tb > 1])
}

#' Cross-Genome Namespaced Tagging
#' @param p_sites Integer vector of plastome coordinates.
#' @param r_sites Integer vector of rDNA coordinates.
#' @return Character vector of tagged loci.
#' @export
tag_sites <- function(p_sites, r_sites) {
  c(
    if (length(p_sites)) paste0("P:", p_sites) else character(0),
    if (length(r_sites)) paste0("R:", r_sites) else character(0)
  )
}

#' Compute Strict Jaccard Set Similarity Index
#' @param set_a Vector.
#' @param set_b Vector.
#' @return Numeric Jaccard index in [0, 1].
#' @export
calc_jaccard <- function(set_a, set_b) {
  u <- length(union(set_a, set_b))
  if (u == 0) return(NA_real_)
  length(intersect(set_a, set_b)) / u
}

cat(" -> [MODULE 1 LOADED] Environment, parameters, and utilities ready.\n\n")