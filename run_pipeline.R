# ==============================================================================
# Pipeline: CeraSNP - Master Pipeline Controller
# Description: Automated execution of discovery, cross-validation, and plotting
# ==============================================================================

rm(list = ls())
cat("\n==================================================================\n")
cat("            CeraSNP Pipeline: Cultivar Identification             \n")
cat("==================================================================\n")

base_dir <- getwd()

# 1. Load Configurations & Utilities
cat("\n>>> [MODULE 1/5] Loading Configurations & Utility Functions...\n")
source(file.path(base_dir, "01_config_and_utils.R"), local = FALSE)

# 2. Compile Optimization Engines
cat("\n>>> [MODULE 2/5] Compiling Two-Tier Greedy Engines...\n")
source(file.path(base_dir, "02_greedy_engine.R"), local = FALSE)

# 3. Discovery Execution & Core Marker Database Construction
cat("\n>>> [MODULE 3/5] Executing Discovery Pipeline & Database Export...\n")
source(file.path(base_dir, "03_discovery_execution.R"), local = FALSE)

# 4. Strict Zero-Leakage LOOCV Validation
cat("\n>>> [MODULE 4/5] Executing Leakage-Free LOOCV & Marker Stability...\n")
source(file.path(base_dir, "04_loocv_validation.R"), local = FALSE)

# 5. Publication Visualizations & Genome Summary
cat("\n>>> [MODULE 5/5] Generating Summary Tables & Figures...\n")
source(file.path(base_dir, "05_visualization.R"), local = FALSE)

cat("\n==================================================================\n")
cat(" [SUCCESS] CeraSNP pipeline finished. All tables & plots in /results\n")
cat("==================================================================\n")