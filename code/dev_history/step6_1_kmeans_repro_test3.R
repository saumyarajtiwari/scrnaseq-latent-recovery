suppressPackageStartupMessages(library(aricode))

manifest <- read.csv("data/processed/embedding_manifest.csv")
master   <- read.csv("data/processed/step4_master_results_table.csv")

test_master <- master[!is.na(master$ari), ]
set.seed(1)
sim_sample  <- test_master[test_master$data_type == "simulated", ][sample(1:sum(test_master$data_type=="simulated"), 5), ]
real_sample <- test_master[test_master$data_type == "real", ][sample(1:sum(test_master$data_type=="real"), 5), ]
test_rows <- rbind(sim_sample, real_sample)

results <- data.frame()

for (i in seq_len(nrow(test_rows))) {
  r <- test_rows[i, ]
  res <- tryCatch({
    if (r$data_type == "simulated") {
      mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                        manifest$run_id == r$run_id & manifest$is_null_control == FALSE, ]
    } else {
      mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                        manifest$data_type == "real", ]
    }
    if (nrow(mrow) == 0) stop("NO MATCH in manifest")

    obj <- readRDS(mrow$file_path[1])
    emb <- obj$embedding
    tg  <- obj$true_group

    # Handle BOTH missing-label conventions found in the data:
    # real NA (muraro, segerstolpe) and empty string "" (tasic2018)
    labeled <- !is.na(tg) & tg != ""
    n_unlabeled <- sum(!labeled)
    emb_l <- emb[labeled, , drop = FALSE]
    tg_l  <- tg[labeled]

    cat(sprintf("row %d [%s/%s]: n_total=%d n_unlabeled=%d n_after_filter=%d n_unique_groups=%d any_na_tg_l=%s any_na_emb_l=%s\n",
                i, r$source, r$method, length(tg), n_unlabeled, length(tg_l),
                length(unique(tg_l)), any(is.na(tg_l)), any(is.na(emb_l))))

    k <- length(unique(tg_l))
    set.seed(42)
    km1 <- kmeans(emb_l, centers = k)$cluster
    cat("  km1 cluster table: ", paste(table(km1), collapse=","),
        " | any NA in km1:", any(is.na(km1)), "\n")

    ari_1 <- aricode::ARI(km1, tg_l)

    set.seed(42)
    km25 <- kmeans(emb_l, centers = k, nstart = 25)$cluster
    ari_25 <- aricode::ARI(km25, tg_l)

    data.frame(source = r$source, method = r$method, run_id = r$run_id,
               n_unlabeled_skipped = n_unlabeled, stored_ari = r$ari,
               recomputed_ari_nstart1 = ari_1, recomputed_ari_nstart25 = ari_25)
  }, error = function(e) {
    cat(sprintf("*** row %d [%s/%s] FAILED: %s\n", i, r$source, r$method, conditionMessage(e)))
    NULL
  })
  if (!is.null(res)) results <- rbind(results, res)
}

cat("\n=== Results (successful rows only) ===\n")
print(results, digits = 4)
