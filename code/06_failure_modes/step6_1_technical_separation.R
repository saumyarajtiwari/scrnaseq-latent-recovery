# Step 6.1 — Detect Technical Separation
# AMI(cluster_labels, batch_id) and AMI(cluster_labels, umi_quartile)
# Flag: AMI > 0.5 against either covariate = Technical Separation confirmed.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(aricode)
  library(parallel)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args

set.seed(42)

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)

## FIX (found during Step 6.1 dry-run diagnostics): real-data rows have
## is_null_control = NA, not FALSE. `x %in% FALSE` silently drops NA rows,
## which was excluding ALL real-data rows from every prior test run.
## Use !(x %in% TRUE) so NA (not-applicable, i.e. real data) correctly passes.
work <- manifest[!(manifest$is_null_control %in% TRUE) & manifest$status == "ok", ]
cat(sprintf("Total candidate rows (non-null-control, status ok): %d\n", nrow(work)))
cat("Breakdown by data_type:\n")
print(table(work$data_type))

if (dry_run) {
  set.seed(1)
  n_sample <- 500
  large_real <- work[work$source %in% c("pbmc68k", "tabula_sapiens_lung"), ]
  n_large <- min(6, nrow(large_real))
  large_real_sample <- large_real[sample(nrow(large_real), n_large), ]
  remaining <- work[!(rownames(work) %in% rownames(large_real_sample)), ]
  rest_sample <- remaining[sample(nrow(remaining), n_sample - n_large), ]
  work <- rbind(large_real_sample, rest_sample)
  cat(sprintf("DRY RUN: subset to %d rows (incl. %d large-real rows)\n",
              nrow(work), n_large))
}

## ---- Stage 1: per-unique-source_file UMI quartile cache ----
umi_cache_dir <- "data/processed/step6_umi_cache"
dir.create(umi_cache_dir, showWarnings = FALSE, recursive = TRUE)

get_umi_quartile <- function(source_file) {
  cache_key <- gsub("[/\\\\]", "_", source_file)
  cache_path <- file.path(umi_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  src <- readRDS(source_file)
  if (!is(src, "SingleCellExperiment")) stop("source_file is not an SCE: ", source_file)
  counts <- assay(src, "counts")
  umi <- Matrix::colSums(counts)
  q <- quantile(umi, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  quart <- cut(umi, breaks = unique(q), include.lowest = TRUE, labels = FALSE)
  result <- list(umi = umi, quartile = quart, cell_id = colData(src)$cell_id)
  saveRDS(result, cache_path)
  rm(src, counts); gc()
  result
}

## ---- Stage 2: main per-row computation ----
process_row <- function(i) {
  r <- work[i, ]
  out <- list(
    data_type = r$data_type, source = r$source, run_id = r$run_id,
    method = r$method, file_path = r$file_path,
    n_cells = NA, n_labeled = NA, n_unlabeled_excluded = NA, k_used = NA,
    n_batches = r$n_batches, batch_ami = NA, batch_ami_note = NA,
    umi_quartile_ami = NA, technical_separation_flag = NA,
    status = "ok", error_msg = NA
  )

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    tg  <- obj$true_group
    bid <- obj$batch_id
    src_file <- obj$source_file

    n_total <- length(tg)
    labeled <- !is.na(tg) & tg != ""
    n_labeled <- sum(labeled)
    n_unlabeled <- n_total - n_labeled

    emb_l <- emb[labeled, , drop = FALSE]
    tg_l  <- tg[labeled]
    bid_l <- bid[labeled]

    k <- length(unique(tg_l))
    if (k < 2) stop("fewer than 2 true_group levels after filtering; cannot cluster")

    km <- kmeans(emb_l, centers = k, nstart = 25)$cluster

    umi_info <- get_umi_quartile(src_file)
    if (length(umi_info$umi) != n_total) {
      stop(sprintf("alignment failure: source_file has %d cells, embedding has %d",
                    length(umi_info$umi), n_total))
    }
    umi_quart_l <- umi_info$quartile[labeled]

    batch_ami <- NA_real_
    batch_note <- NA_character_
    n_batch_levels <- length(unique(bid_l))
    if (n_batch_levels <= 1) {
      batch_note <- "excluded_no_batch_variation"
    } else {
      bid_int <- as.integer(factor(bid_l))
      batch_ami <- aricode::AMI(km, bid_int)
      if (n_batch_levels > 50) batch_note <- "high_cardinality_batch_covariate"
    }

    umi_int <- as.integer(factor(umi_quart_l))
    umi_ami <- aricode::AMI(km, umi_int)

    tech_sep <- isTRUE(batch_ami > 0.5) || isTRUE(umi_ami > 0.5)

    list(n_cells = n_total, n_labeled = n_labeled, n_unlabeled_excluded = n_unlabeled,
         k_used = k, batch_ami = batch_ami, batch_ami_note = batch_note,
         umi_quartile_ami = umi_ami, technical_separation_flag = tech_sep,
         status = "ok", error_msg = NA_character_)
  }, error = function(e) {
    list(n_cells = NA, n_labeled = NA, n_unlabeled_excluded = NA, k_used = NA,
         batch_ami = NA, batch_ami_note = NA, umi_quartile_ami = NA,
         technical_separation_flag = NA, status = "error",
         error_msg = conditionMessage(e))
  })

  out[names(res)] <- res
  out
}

cat("Starting mclapply (4 workers, mc.preschedule=FALSE)...\n")
t0 <- Sys.time()
results_list <- mclapply(seq_len(nrow(work)), process_row,
                          mc.cores = 4, mc.preschedule = FALSE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

results_df <- do.call(rbind.data.frame, results_list)

n_errors <- sum(results_df$status == "error")
cat(sprintf("Rows processed: %d | errors: %d\n", nrow(results_df), n_errors))
if (n_errors > 0) {
  cat("Error summary:\n")
  print(table(results_df$error_msg))
}

out_path <- if (dry_run) "data/processed/step6_1_technical_separation_DRYRUN.csv" else
                          "data/processed/step6_1_technical_separation.csv"
write.csv(results_df, out_path, row.names = FALSE)
cat("Written to:", out_path, "\n")
