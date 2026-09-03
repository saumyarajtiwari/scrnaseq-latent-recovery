# Step 1.6 Phase A — Null-control replication manifest
# Builds data/simulated/null_control_grid.csv: 45 rows tracking 3 replicates
# (1 existing + 2 new) for each of 5 sparsity levels x 3 simulators.
# Replicate seeds offset by +100000 (rep2) / +200000 (rep3), consistent with
# the SymSim batch-effect seed-offset convention (run_id + 10000/batch) from Step 1.5.

param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors = FALSE)
null_rows  <- param_grid[param_grid$is_null_control == TRUE, ]

simulators <- c("splatter", "scdesign3", "symsim")

manifest <- do.call(rbind, lapply(simulators, function(sim) {
  do.call(rbind, lapply(seq_len(nrow(null_rows)), function(i) {
    row <- null_rows[i, ]
    do.call(rbind, lapply(1:3, function(rep) {
      seed <- if (rep == 1) row$run_id else row$run_id + (rep - 1) * 100000
      is_existing <- (rep == 1)
      file_path <- if (is_existing) {
        sprintf("data/simulated/%s/%s_run_%d.rds", sim, sim, row$run_id)
      } else {
        sprintf("data/simulated/null_control/%s/%s_null_sparsity%s_rep%d.rds",
                sim, sim, sprintf("%.2f", row$sparsity), rep)
      }
      data.frame(
        simulator   = sim,
        base_run_id = row$run_id,
        sparsity    = row$sparsity,
        replicate   = rep,
        seed        = seed,
        file_path   = file_path,
        is_new      = !is_existing,
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

stopifnot(nrow(manifest) == 45)
stopifnot(sum(manifest$is_new) == 30)
stopifnot(sum(!manifest$is_new) == 15)

dir.create("data/simulated/null_control", showWarnings = FALSE)
for (sim in simulators) {
  dir.create(file.path("data/simulated/null_control", sim), showWarnings = FALSE, recursive = TRUE)
}

out_path <- "data/simulated/null_control_grid.csv"
write.csv(manifest, out_path, row.names = FALSE)

cat("=== Null-control manifest summary ===\n")
cat("Total rows       :", nrow(manifest), "\n")
cat("Existing (rep 1) :", sum(!manifest$is_new), "\n")
cat("New (rep 2/3)    :", sum(manifest$is_new), "\n")
cat("Manifest written to:", out_path, "\n")
