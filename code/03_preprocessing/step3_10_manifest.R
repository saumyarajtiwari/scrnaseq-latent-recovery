suppressPackageStartupMessages({ library(parallel) })

METHOD_MAP <- list(
  pca_raw            = "rawpca",
  pca_libnorm        = "libnormpca",
  pca_log            = "logpca",
  pca_shiftedlog     = "shiftedlogpca",
  pca_sctransform_v2 = "sctpca",
  pca_glmpca         = "glmpca"
)

N_WORKERS   <- 4
BATCH_SIZE  <- 2000
OUT_DIR     <- "data/processed"
MANIFEST_CSV   <- file.path(OUT_DIR, "step3_10_embedding_manifest.csv")
CHECKPOINT_RDS <- file.path(OUT_DIR, "step3_10_checkpoint.rds")
PROGRESS_LOG   <- file.path(OUT_DIR, "step3_10_progress.log")

param_grid <- read.csv("data/simulated/param_grid.csv", stringsAsFactors = FALSE)
real_inv   <- read.csv("data/real_data_inventory.csv", stringsAsFactors = FALSE)

log_msg <- function(...) {
  msg <- paste(format(Sys.time()), "-", ...)
  cat(msg, "\n")
  write(msg, PROGRESS_LOG, append = file.exists(PROGRESS_LOG))
}

build_file_list <- function() {
  rows <- list()
  for (method_dir in names(METHOD_MAP)) {
    methodtag <- METHOD_MAP[[method_dir]]

    real_files <- list.files(file.path("data/processed", method_dir, "real"),
                              pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
    for (f in real_files) {
      bn <- basename(f)
      pat <- paste0("^(.+)_harmonized_", methodtag, "\\.rds$")
      m <- regmatches(bn, regexec(pat, bn))[[1]]
      src <- if (length(m) == 2) m[2] else NA_character_
      rows[[length(rows) + 1]] <- data.frame(
        file_path = f, method = method_dir, data_type = "real",
        source = src, run_id = NA_integer_, stringsAsFactors = FALSE
      )
    }

    for (sim in c("scdesign3", "splatter", "symsim")) {
      sim_files <- list.files(file.path("data/processed", method_dir, "simulated", sim),
                               pattern = "\\.rds$", full.names = TRUE)
      for (f in sim_files) {
        bn <- basename(f)
        pat <- paste0("^", sim, "_sce_run_([0-9]{5})_", methodtag, "\\.rds$")
        m <- regmatches(bn, regexec(pat, bn))[[1]]
        rid <- if (length(m) == 2) as.integer(m[2]) else NA_integer_
        rows[[length(rows) + 1]] <- data.frame(
          file_path = f, method = method_dir, data_type = "simulated",
          source = sim, run_id = rid, stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

log_msg("Building raw file list...")
file_list <- build_file_list()
log_msg("Total files found:", nrow(file_list))

n_unparsed <- sum(is.na(file_list$source) | (file_list$data_type == "simulated" & is.na(file_list$run_id)))
if (n_unparsed > 0) {
  log_msg("WARNING:", n_unparsed, "files failed filename parsing - inspect before proceeding")
  write.csv(file_list[is.na(file_list$source) | (file_list$data_type == "simulated" & is.na(file_list$run_id)), ],
            file.path(OUT_DIR, "step3_10_UNPARSED_files.csv"), row.names = FALSE)
}

check_one <- function(fp) {
  tryCatch({
    d <- readRDS(fp)
    list(
      status = "ok",
      has_embedding = !is.null(d$embedding),
      has_loadings = !is.null(d$loadings),
      embedding_dim = if (!is.null(d$embedding)) paste(dim(d$embedding), collapse = "x") else NA,
      loadings_dim = if (!is.null(d$loadings)) paste(dim(d$loadings), collapse = "x") else NA,
      nv_used = if (!is.null(d$nv_used)) d$nv_used else NA,
      any_na_embedding = if (!is.null(d$embedding)) anyNA(d$embedding) else NA,
      # Item 5 fix: distinguish GLM-PCA's chunked fixed-loadings projection
      # (pbmc68k, tabula_sapiens_lung -- validated only on Baron, ARI
      # 0.652->0.421 degradation) from every other file's native fit.
      # internal_method carries the .rds file's own `method` field verbatim
      # ("glmpca_chunked_projection" vs "glmpca_nb" vs other methods' own
      # tags); NA-safe for methods without this field.
      internal_method = if (!is.null(d$method)) d$method else NA_character_,
      has_validation_note = !is.null(d$validation_note),
      error_msg = NA_character_
    )
  }, error = function(e) {
    list(status = "error", has_embedding = NA, has_loadings = NA,
         embedding_dim = NA, loadings_dim = NA, nv_used = NA,
         any_na_embedding = NA, internal_method = NA_character_,
         has_validation_note = NA, error_msg = conditionMessage(e))
  })
}

n_total <- nrow(file_list)
n_batches <- ceiling(n_total / BATCH_SIZE)

completed_batches <- if (file.exists(CHECKPOINT_RDS)) readRDS(CHECKPOINT_RDS) else integer(0)
log_msg("Resuming:", length(completed_batches), "of", n_batches, "batches already done")

for (b in seq_len(n_batches)) {
  if (b %in% completed_batches) next
  idx_start <- (b - 1) * BATCH_SIZE + 1
  idx_end <- min(b * BATCH_SIZE, n_total)
  batch_paths <- file_list$file_path[idx_start:idx_end]

  results <- mclapply(batch_paths, check_one, mc.cores = N_WORKERS, mc.preschedule = FALSE)

  batch_df <- file_list[idx_start:idx_end, ]
  batch_df$status <- vapply(results, function(r) r$status, character(1))
  batch_df$has_embedding <- vapply(results, function(r) as.character(r$has_embedding), character(1))
  batch_df$has_loadings <- vapply(results, function(r) as.character(r$has_loadings), character(1))
  batch_df$embedding_dim <- vapply(results, function(r) ifelse(is.na(r$embedding_dim), "", r$embedding_dim), character(1))
  batch_df$loadings_dim <- vapply(results, function(r) ifelse(is.na(r$loadings_dim), "", r$loadings_dim), character(1))
  batch_df$nv_used <- vapply(results, function(r) ifelse(is.na(r$nv_used), NA_character_, as.character(r$nv_used)), character(1))
  batch_df$any_na_embedding <- vapply(results, function(r) as.character(r$any_na_embedding), character(1))
  batch_df$internal_method <- vapply(results, function(r) ifelse(is.na(r$internal_method), "", r$internal_method), character(1))
  batch_df$has_validation_note <- vapply(results, function(r) as.character(r$has_validation_note), character(1))
  batch_df$error_msg <- vapply(results, function(r) ifelse(is.na(r$error_msg), "", r$error_msg), character(1))

  write.table(batch_df, MANIFEST_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(MANIFEST_CSV), append = file.exists(MANIFEST_CSV))

  completed_batches <- c(completed_batches, b)
  saveRDS(completed_batches, CHECKPOINT_RDS)

  n_errors_batch <- sum(batch_df$status == "error")
  log_msg(sprintf("Batch %d/%d done (%d files). Errors this batch: %d. Total processed: %d/%d",
                   b, n_batches, length(batch_paths), n_errors_batch, idx_end, n_total))
}

log_msg("=== INTEGRITY SCAN COMPLETE ===")

manifest <- read.csv(MANIFEST_CSV, stringsAsFactors = FALSE)

real_inv_renamed <- real_inv
names(real_inv_renamed)[names(real_inv_renamed) == "file_path"] <- "source_sce_file_path"
real_inv_renamed$dataset_name[real_inv_renamed$dataset_name == "ts_lung"] <- "tabula_sapiens_lung"

sim_rows <- manifest[manifest$data_type == "simulated", ]
sim_joined <- merge(sim_rows, param_grid, by = "run_id", all.x = TRUE)

real_rows <- manifest[manifest$data_type == "real", ]
real_joined <- merge(real_rows, real_inv_renamed, by.x = "source", by.y = "dataset_name", all.x = TRUE)

all_cols <- union(names(sim_joined), names(real_joined))
for (col in setdiff(all_cols, names(sim_joined))) sim_joined[[col]] <- NA
for (col in setdiff(all_cols, names(real_joined))) real_joined[[col]] <- NA
final_manifest <- rbind(sim_joined[, all_cols], real_joined[, all_cols])

FINAL_OUT <- file.path(OUT_DIR, "embedding_manifest.csv")
write.csv(final_manifest, FINAL_OUT, row.names = FALSE)

n_errors_total <- sum(manifest$status == "error")
log_msg("Final manifest written to", FINAL_OUT, "-", nrow(final_manifest), "rows,", n_errors_total, "total errors/corrupted files")
if (n_errors_total > 0) {
  write.csv(manifest[manifest$status == "error", ], file.path(OUT_DIR, "step3_10_CORRUPTED_files.csv"), row.names = FALSE)
  log_msg("Corrupted/unreadable file list written to step3_10_CORRUPTED_files.csv - REVIEW BEFORE PROCEEDING TO STEP 4")
}
