# =============================================================================
# install_symsim_patched.R
#
# Clones SymSim from YosefLab/SymSim, applies the BiocGenerics rank() patch,
# and installs from patched source.
#
# WHY THIS PATCH IS NEEDED:
# BiocGenerics (loaded via scran/scater) registers an S4 generic for rank().
# SymSim's Get_params() uses alply() + do.call(c, temp) which produces an S4
# list-like object. When rank() is called on this object, S4 method dispatch
# validates the argument before calling the method and fails with:
#   "x must be an array of at least two dimensions"
# The fix replaces alply/do.call(c,...) with as.numeric(t(as.matrix(X))) —
# mathematically identical but producing a plain numeric vector that base::rank
# handles correctly. Both rank() calls are qualified as base::rank() to bypass
# S4 dispatch entirely.
#
# PATCH LOCATIONS in simulation_functions.R:
#   Line ~96:  Get_params()         — alply/do.call(c)/rank block
#   Line ~666: SimulateTrueCounts() — rank(rowSums(...)) for HGE genes
#
# USAGE: Rscript R/01_simulation/install_symsim_patched.R
# =============================================================================

dest <- "vendor/SymSim_patched"

if (!dir.exists(dest)) {
  cat("Cloning SymSim source...\n")
  system2("git", c("clone", "https://github.com/YosefLab/SymSim.git", dest))
} else {
  cat("vendor/SymSim_patched already exists — skipping clone.\n")
}

target <- file.path(dest, "R", "simulation_functions.R")
content <- readLines(target, warn=FALSE)
content_orig <- content

# Fix 1: replace alply/do.call(c,temp)/rank block in Get_params
alply_line <- grep("temp <- alply(X, 1, function(Y){Y})", content, fixed=TRUE)
dc_line    <- grep("values <- do.call(c,temp)",           content, fixed=TRUE)
rank_line1 <- grep("^    ranks <- rank(values)",           content, fixed=TRUE)

if (length(alply_line) == 1 && length(dc_line) == 1 && length(rank_line1) == 1) {
  content[alply_line] <- "    values <- as.numeric(t(as.matrix(X)))"
  content[dc_line]    <- ""
  content[rank_line1] <- "    ranks <- base::rank(values)"
  cat("Fix 1 applied: Get_params alply/do.call(c)/rank block\n")
} else if (!any(grepl("base::rank(values)", content, fixed=TRUE))) {
  stop("Fix 1: target lines not found — check SymSim version")
} else {
  cat("Fix 1: already applied\n")
}

# Fix 2: rank_sum in SimulateTrueCounts
rank_sum_line <- grep("rank_sum <- rank(rowSums(", content, fixed=TRUE)

if (length(rank_sum_line) == 1) {
  content[rank_sum_line] <- paste0(
    "    rank_sum <- base::rank(as.numeric(rowSums(",
    "as.matrix(params[[3]][chosen_hge, , drop=FALSE]))))")
  cat("Fix 2 applied: SimulateTrueCounts rank_sum line\n")
} else if (!any(grepl("base::rank(as.numeric", content, fixed=TRUE))) {
  stop("Fix 2: target line not found — check SymSim version")
} else {
  cat("Fix 2: already applied\n")
}

if (!identical(content, content_orig)) {
  writeLines(content, target)
  cat("Patched file written.\n")
}

cat("\nInstalling patched SymSim...\n")
install.packages(dest, repos=NULL, type="source")
cat("\nDone. Verify with: library(SymSim)\n")
