# ==============================================================================
# Pipeline: CeraSNP - Plastome & rDNA Core SNP Discovery for Cerasus Cultivars
# Module:   02_greedy_engine.R (Two-Tier Optimization Engines & KASP Assay QC)
# Version:  1.0.0 (GitHub Release Edition)
# License:  MIT License
# ==============================================================================

cat("\n============================================================\n")
cat(" [MODULE 2] Loading Two-Tier Adaptive Engines & In-Silico QC...\n")
cat("============================================================\n")
flush.console()

# ------------------------------------------------------------------------------
# 1. Helper: Site-Level Population Genetics Metrics
# ------------------------------------------------------------------------------
get_site_genetic_metrics <- function(col_vec) {
  valid_bases <- c("A", "C", "G", "T")
  v <- toupper(col_vec)
  v <- v[v %in% valid_bases]
  n_valid <- length(v)
  
  if (n_valid == 0) {
    return(list(
      Alleles = "", Allele_Count = 0, Poly_Type = "Empty",
      Major_Allele = "", Major_Freq = 0, MAF = 0, Entropy = 0
    ))
  }
  
  tbl <- sort(table(v), decreasing = TRUE)
  alleles <- names(tbl)
  allele_str <- paste(alleles, collapse = "/")
  n_alleles <- length(alleles)
  
  poly_type <- if (n_alleles == 2) "Biallelic" else if (n_alleles > 2) "Multiallelic" else "Monomorphic"
  freqs <- as.numeric(tbl) / n_valid
  
  major_allele <- alleles[1]
  major_freq   <- round(freqs[1], 4)
  maf          <- if (n_alleles >= 2) round(min(freqs), 4) else 0.0000
  
  raw_entropy <- -sum(freqs * log2(freqs + 1e-9))
  entropy     <- round(max(0, raw_entropy), 4)
  
  list(
    Alleles      = allele_str,
    Allele_Count = n_alleles,
    Poly_Type    = poly_type,
    Major_Allele = major_allele,
    Major_Freq   = major_freq,
    MAF          = maf,
    Entropy      = entropy
  )
}

# ------------------------------------------------------------------------------
# 2. Tier 1: Maternal Plastome Backbone Forward Selection Engine
# ------------------------------------------------------------------------------
adaptive_greedy_tier1 <- function(mat, candidate_sites, max_trials = 40, verbose = FALSE) {
  n <- nrow(mat)
  selected <- integer(0)
  remaining <- candidate_sites
  current <- rep("", n)
  current_unique <- length(unique(current))
  hist <- list()
  stop_reason <- NULL

  if (!length(remaining)) {
    return(list(sites = selected, profile = current, history = data.frame(), stop_reason = "No_Candidate_Sites"))
  }
  if (current_unique == n) {
    return(list(sites = selected, profile = current, history = data.frame(), stop_reason = "All_Profiles_Unique_Initial"))
  }

  n_steps <- min(max_trials, length(remaining))

  for (step in seq_len(n_steps)) {
    score_vec <- vapply(remaining, function(pos) {
      trial_profile <- paste0(current, "_", mat[, pos])
      length(unique(trial_profile))
    }, integer(1))

    max_val <- max(score_vec)
    gain <- max_val - current_unique

    # Adaptive termination: saturation reached
    if (gain <= 0) {
      stop_reason <- "Marginal_Gain_Nonpositive (Adaptive_Convergence)"
      if (verbose) cat(sprintf("    [Tier 1 Convergence] Step %d: Plastome gain saturated (Delta Gain <= 0).\n", step))
      break
    }

    # Shannon entropy tie-breaking
    best_candidates <- remaining[which(score_vec == max_val)]
    if (length(best_candidates) > 1) {
      entropies <- vapply(best_candidates, function(pos) {
        p <- as.numeric(table(mat[, pos])) / n
        -sum(p * log2(p + 1e-9))
      }, numeric(1))
      site <- best_candidates[which.max(entropies)]
    } else {
      site <- best_candidates[1]
    }

    metrics <- get_site_genetic_metrics(mat[, site])

    selected <- c(selected, site)
    current <- paste0(current, "_", mat[, site])
    current_unique <- max_val
    remaining <- setdiff(remaining, site)

    hist[[step]] <- data.frame(
      Step            = step,
      Align_Pos       = site,
      Alleles         = metrics$Alleles,
      Poly_Type       = metrics$Poly_Type,
      Major_Allele    = metrics$Major_Allele,
      Major_Freq      = metrics$Major_Freq,
      MAF             = metrics$MAF,
      Shannon_Entropy = metrics$Entropy,
      Unique_Profiles = current_unique,
      Marginal_Gain   = gain,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat(sprintf("    [Step %02d] Selected Plastome SNP: Pos_%-6d | Type: %-12s | Alleles: %-5s | MAF: %5.3f | Haplotypes: %2d (+%d)\n", 
                  step, site, metrics$Poly_Type, metrics$Alleles, metrics$MAF, current_unique, gain))
      flush.console()
    }

    if (current_unique == n) {
      stop_reason <- "All_Profiles_Unique"
      break
    }
  }

  if (is.null(stop_reason)) {
    stop_reason <- if (!length(remaining)) "Candidate_Pool_Exhausted" else "Max_Trials_Reached"
  }

  list(
    sites       = selected,
    profile     = current,
    history     = if (length(hist)) do.call(rbind, hist) else data.frame(),
    stop_reason = stop_reason
  )
}

# ------------------------------------------------------------------------------
# 3. Tier 2: Nuclear rDNA Residual Collision Minimization Engine
# ------------------------------------------------------------------------------
adaptive_greedy_tier2 <- function(mat_r, candidate_sites, background_profile, max_trials = 20, verbose = FALSE) {
  n <- nrow(mat_r)
  selected <- integer(0)
  remaining <- candidate_sites
  current_r <- rep("", n)
  current_combined <- paste0(background_profile, "_rDNA_", current_r)

  calc_pairs <- function(prof) {
    tb <- table(prof)
    sum(tb * (tb - 1) / 2)
  }
  calc_samples <- function(prof) {
    tb <- table(prof)
    if (!length(tb)) return(0L)
    sum(tb[tb > 1])
  }

  current_pairs <- calc_pairs(current_combined)
  current_collisions <- calc_samples(current_combined)
  hist <- list()
  stop_reason <- NULL

  if (!length(remaining)) {
    return(list(sites = selected, rDNA_profile = current_r, combined_profile = current_combined, 
                history = data.frame(), stop_reason = "No_Candidate_Sites"))
  }
  if (current_pairs == 0) {
    return(list(sites = selected, rDNA_profile = current_r, combined_profile = current_combined, 
                history = data.frame(), stop_reason = "Zero_Initial_Collisions"))
  }

  n_steps <- min(max_trials, length(remaining))

  for (step in seq_len(n_steps)) {
    score_pairs <- vapply(remaining, function(pos) {
      trial_r <- paste0(current_r, "_", mat_r[, pos])
      trial_combined <- paste0(background_profile, "_rDNA_", trial_r)
      calc_pairs(trial_combined)
    }, numeric(1))

    min_val <- min(score_pairs)
    pair_reduction <- current_pairs - min_val

    # Adaptive termination: collisions irreducible
    if (pair_reduction <= 0) {
      stop_reason <- "Marginal_Collision_Reduction_Nonpositive (Adaptive_Convergence)"
      if (verbose) cat(sprintf("    [Tier 2 Convergence] Step %d: Residual collisions irreducible (Delta C <= 0).\n", step))
      break
    }

    best_candidates <- remaining[which(score_pairs == min_val)]
    if (length(best_candidates) > 1) {
      entropies <- vapply(best_candidates, function(pos) {
        p <- as.numeric(table(mat_r[, pos])) / n
        -sum(p * log2(p + 1e-9))
      }, numeric(1))
      site <- best_candidates[which.max(entropies)]
    } else {
      site <- best_candidates[1]
    }

    metrics <- get_site_genetic_metrics(mat_r[, site])

    selected <- c(selected, site)
    current_r <- paste0(current_r, "_", mat_r[, site])
    current_combined <- paste0(background_profile, "_rDNA_", current_r)
    current_pairs <- min_val
    current_collisions <- calc_samples(current_combined)
    remaining <- setdiff(remaining, site)

    hist[[step]] <- data.frame(
      Step                              = step,
      Align_Pos                         = site,
      Alleles                           = metrics$Alleles,
      Poly_Type                         = metrics$Poly_Type,
      Major_Allele                      = metrics$Major_Allele,
      Major_Freq                        = metrics$Major_Freq,
      MAF                               = metrics$MAF,
      Shannon_Entropy                   = metrics$Entropy,
      Collision_Pairs                   = current_pairs,
      Collision_Samples                 = current_collisions,
      Marginal_Collision_Pair_Reduction = pair_reduction,
      stringsAsFactors                  = FALSE
    )

    if (verbose) {
      cat(sprintf("    [Step %02d] Selected rDNA SNP: Pos_%-6d | Type: %-12s | Alleles: %-5s | MAF: %5.3f | Residual Pairs: %2d (-%d)\n", 
                  step, site, metrics$Poly_Type, metrics$Alleles, metrics$MAF, current_pairs, pair_reduction))
      flush.console()
    }

    if (current_pairs == 0) {
      stop_reason <- "All_Collisions_Resolved"
      break
    }
  }

  if (is.null(stop_reason)) {
    stop_reason <- if (!length(remaining)) "Candidate_Pool_Exhausted" else "Max_Trials_Reached"
  }

  list(
    sites            = selected,
    rDNA_profile     = current_r,
    combined_profile = current_combined,
    history          = if (length(hist)) do.call(rbind, hist) else data.frame(),
    stop_reason      = stop_reason
  )
}

# ------------------------------------------------------------------------------
# 4. In-Silico KASP Quality Control Assay
# ------------------------------------------------------------------------------
generate_kasp_qc <- function(mat, sites, candidate_all, prefix, genome_type, flank = 50) {
  if (!length(sites)) return(data.frame())
  L <- ncol(mat)
  valid_bases <- c("A", "C", "G", "T")
  rows <- vector("list", length(sites))

  for (i in seq_along(sites)) {
    pos <- sites[i]
    st <- max(1, pos - flank)
    ed <- min(L, pos + flank)

    col_bases <- toupper(mat[, pos])
    col_bases <- col_bases[col_bases %in% valid_bases]
    alleles <- sort(unique(col_bases))
    allele_str <- paste(alleles, collapse = "/")
    is_biallelic <- length(alleles) == 2

    flank_indices <- setdiff(st:ed, pos)
    sec_snps <- sum(flank_indices %in% candidate_all)

    get_consensus <- function(mm) {
      if (!ncol(mm)) return("")
      paste0(apply(mm, 2, function(col) {
        v <- toupper(col)
        v <- v[v %in% valid_bases]
        if (!length(v)) return("N")
        tbl <- sort(table(v), decreasing = TRUE)
        if (length(tbl) > 1 && tbl[1] == tbl[2]) "N" else names(tbl)[1]
      }), collapse = "")
    }

    left  <- if (pos > 1) get_consensus(mat[, st:(pos - 1), drop = FALSE]) else ""
    right <- if (pos < L) get_consensus(mat[, (pos + 1):ed, drop = FALSE]) else ""
    left  <- gsub("-", "", left, fixed = TRUE)
    right <- gsub("-", "", right, fixed = TRUE)

    chars_l <- if (nchar(left)) strsplit(left, "", fixed = TRUE)[[1]] else character(0)
    chars_r <- if (nchar(right)) strsplit(right, "", fixed = TRUE)[[1]] else character(0)
    all_chars <- c(chars_l, chars_r)

    valid_gc <- all_chars[all_chars %in% valid_bases]
    gc <- if (length(valid_gc)) round(mean(valid_gc %in% c("G", "C")) * 100, 1) else NA_real_
    
    poly_left  <- if (length(chars_l)) any(rle(chars_l)$lengths >= 5) else FALSE
    poly_right <- if (length(chars_r)) any(rle(chars_r)$lengths >= 5) else FALSE
    has_poly   <- poly_left || poly_right

    rows[[i]] <- data.frame(
      Marker_ID                   = sprintf("%s_%02d", prefix, i),
      Genome_Type                 = genome_type,
      Align_Pos                   = pos,
      Alleles                     = allele_str,
      Allele_Count                = length(alleles),
      Platform_Compatibility      = ifelse(is_biallelic, "Biallelic_Standard_KASP", "Multiallelic_Custom_Assay"),
      Flank_Window_bp             = flank,
      Consensus_Left_Flank        = left,
      Consensus_Right_Flank       = right,
      Flank_GC_Pct                = gc,
      Flank_Secondary_SNPs_Count  = sec_snps,
      Flank_Secondary_SNP_Warn    = ifelse(sec_snps > 0, "WARN_SECONDARY_POLY", "PASS"),
      Homopolymer_Warn            = ifelse(has_poly, "WARN_POLY_RUN", "PASS"),
      stringsAsFactors            = FALSE
    )
  }
  do.call(rbind, rows)
}

cat(" -> [MODULE 2 COMPILED] Optimization engines and KASP QC functions ready.\n\n")