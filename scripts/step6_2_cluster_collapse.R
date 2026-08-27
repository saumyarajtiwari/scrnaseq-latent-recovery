# Step 6.2 — Detect Cluster Collapse
# Pairwise centroid distances between true biological populations, compared
# against mean intra-cluster distance. silhouette < 0.2 (Step 4.6) is the
# operational trigger; centroid-pair detail identifies WHICH populations collapsed.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% args

set.seed(42)

manifest <- read.csv("data/processed/embedding_manifest.csv", stringsAsFactors = FALSE)
work <- manifest[!(manifest$is_null_control %in% TRUE) & manifest$status == "ok", ]
cat(sprintf("Total candidate rows: %d\n", nrow(work)))

master <- read.csv("data/processed/step4_master_results_table.csv", stringsAsFactors = FALSE)
master_key <- paste(master$data_type, master$source, master$run_id, master$method)
sil_lookup <- setNames(master$silhouette, master_key)

if (dry_run) {
  set.seed(1)
  n_sample <- 500
  large_real <- work[work$source %in% c("pbmc68k", "tabula_sapiens_lung", "tasic2018"), ]
  n_large <- min(9, nrow(large_real))
  large_real_sample <- large_real[sample(nrow(large_real), n_large), ]
  remaining <- work[!(rownames(work) %in% rownames(large_real_sample)), ]
  rest_sample <- remaining[sample(nrow(remaining), n_sample - n_large), ]
  work <- rbind(large_real_sample, rest_sample)
  cat(sprintf("DRY RUN: subset to %d rows (incl. %d large-real rows)\n", nrow(work), n_large))
}

process_row <- function(i) {
  r <- work[i, ]
  key <- paste(r$data_type, r$source, r$run_id, r$method)
  sil <- sil_lookup[[key]]

  summary_out <- list(
    data_type = r$data_type, source = r$source, run_id = r$run_id, method = r$method,
    file_path = r$file_path, n_cells = NA, n_labeled = NA, n_groups_used = NA,
    silhouette = ifelse(is.null(sil), NA, sil),
    silhouette_trigger = NA, n_pairs_total = NA, n_pairs_collapsed = NA,
    min_margin = NA, worst_pair = NA_character_,
    status = "ok", error_msg = NA_character_
  )
  detail_rows <- list()

  res <- tryCatch({
    obj <- readRDS(r$file_path)
    emb <- obj$embedding
    tg  <- obj$true_group

    n_total <- length(tg)
    labeled <- !is.na(tg) & tg != ""
    emb_l <- emb[labeled, , drop = FALSE]
    tg_l  <- as.character(tg[labeled])
    groups <- unique(tg_l)
    k <- length(groups)
    if (k < 2) stop("fewer than 2 true_group levels after filtering; cannot compare pairs")

    centroids <- t(sapply(groups, function(g) colMeans(emb_l[tg_l == g, , drop = FALSE])))
    rownames(centroids) <- groups

    intra_dist <- sapply(groups, function(g) {
      pts <- emb_l[tg_l == g, , drop = FALSE]
      n_g <- nrow(pts)
      if (n_g <= 1) return(0)
      cen <- centroids[g, ]
      diffs <- sweep(pts, 2, cen, "-")
      mean(sqrt(rowSums(diffs^2)))
    })
    names(intra_dist) <- groups
    group_sizes <- table(tg_l)

    pairs <- combn(groups, 2, simplify = FALSE)
    n_pairs <- length(pairs)
    collapsed_list <- list()
    margins <- numeric(n_pairs)

    for (p_idx in seq_along(pairs)) {
      g1 <- pairs[[p_idx]][1]; g2 <- pairs[[p_idx]][2]
      cd <- sqrt(sum((centroids[g1, ] - centroids[g2, ])^2))
      thresh <- mean(c(intra_dist[[g1]], intra_dist[[g2]]))
      margin <- cd - thresh
      margins[p_idx] <- margin
      if (cd < thresh) {
        collapsed_list[[length(collapsed_list) + 1]] <- data.frame(
          data_type = r$data_type, source = r$source, run_id = r$run_id, method = r$method,
          group1 = g1, group2 = g2, n_group1 = as.integer(group_sizes[g1]),
          n_group2 = as.integer(group_sizes[g2]), centroid_dist = cd,
          intra_dist_1 = intra_dist[[g1]], intra_dist_2 = intra_dist[[g2]],
          threshold = thresh, margin = margin
        )
      }
    }

    n_collapsed <- length(collapsed_list)
    worst_idx <- which.min(margins)
    worst_pair_str <- paste(pairs[[worst_idx]][1], pairs[[worst_idx]][2], sep = "|")

    list(n_cells = n_total, n_labeled = sum(labeled), n_groups_used = k,
         silhouette_trigger = isTRUE(sil < 0.2),
         n_pairs_total = n_pairs, n_pairs_collapsed = n_collapsed,
         min_margin = margins[worst_idx], worst_pair = worst_pair_str,
         status = "ok", error_msg = NA_character_,
         detail = if (n_collapsed > 0) do.call(rbind, collapsed_list) else NULL)
  }, error = function(e) {
    list(n_cells = NA, n_labeled = NA, n_groups_used = NA, silhouette_trigger = NA,
         n_pairs_total = NA, n_pairs_collapsed = NA, min_margin = NA, worst_pair = NA_character_,
         status = "error", error_msg = conditionMessage(e), detail = NULL)
  })

  detail_df <- res$detail
  res$detail <- NULL
  summary_out[names(res)] <- res
  list(summary = summary_out, detail = detail_df)
}

cat("Starting mclapply (4 workers, mc.preschedule=FALSE)...\n")
t0 <- Sys.time()
results_list <- mclapply(seq_len(nrow(work)), process_row, mc.cores = 4, mc.preschedule = FALSE)
cat("Elapsed:", format(Sys.time() - t0), "\n")

summary_df <- do.call(rbind.data.frame, lapply(results_list, `[[`, "summary"))
detail_list <- lapply(results_list, `[[`, "detail")
detail_list <- detail_list[!sapply(detail_list, is.null)]
detail_df <- if (length(detail_list) > 0) do.call(rbind, detail_list) else NULL

n_errors <- sum(summary_df$status == "error")
cat(sprintf("Rows processed: %d | errors: %d\n", nrow(summary_df), n_errors))
if (n_errors > 0) print(table(summary_df$error_msg))
cat(sprintf("Rows with silhouette_trigger TRUE: %d\n", sum(summary_df$silhouette_trigger, na.rm=TRUE)))
cat(sprintf("Total collapsed pairs found: %d\n", if (is.null(detail_df)) 0 else nrow(detail_df)))

suffix <- if (dry_run) "_DRYRUN" else ""
write.csv(summary_df, paste0("data/processed/step6_2_cluster_collapse", suffix, ".csv"), row.names = FALSE)
if (!is.null(detail_df)) {
  write.csv(detail_df, paste0("data/processed/step6_2_collapsed_pairs_detail", suffix, ".csv"), row.names = FALSE)
}
cat("Written.\n")
