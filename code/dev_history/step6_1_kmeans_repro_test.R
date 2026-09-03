suppressPackageStartupMessages(library(aricode))

manifest <- read.csv("data/processed/embedding_manifest.csv")
master   <- read.csv("data/processed/step4_master_results_table.csv")

# Pick 5 simulated + 5 real rows with a stored ARI to test against
test_master <- master[!is.na(master$ari), ]
set.seed(1)
sim_sample  <- test_master[test_master$data_type == "simulated", ][sample(1:sum(test_master$data_type=="simulated"), 5), ]
real_sample <- test_master[test_master$data_type == "real", ][sample(1:sum(test_master$data_type=="real"), 5), ]
test_rows <- rbind(sim_sample, real_sample)

results <- data.frame()

for (i in seq_len(nrow(test_rows))) {
  r <- test_rows[i, ]
  if (r$data_type == "simulated") {
    mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                      manifest$run_id == r$run_id & manifest$is_null_control == FALSE, ]
  } else {
    mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                      manifest$data_type == "real", ]
  }
  if (nrow(mrow) == 0) { cat("NO MATCH for row", i, "\n"); next }

  obj <- readRDS(mrow$file_path[1])
  emb <- obj$embedding
  tg  <- obj$true_group
  k   <- length(unique(tg))

  set.seed(42)
  km_nstart1  <- kmeans(emb, centers = k)$cluster
  set.seed(42)
  km_nstart25 <- kmeans(emb, centers = k, nstart = 25)$cluster

  ari_1  <- aricode::ARI(km_nstart1, tg)
  ari_25 <- aricode::ARI(km_nstart25, tg)

  results <- rbind(results, data.frame(
    source = r$source, method = r$method, run_id = r$run_id,
    stored_ari = r$ari, recomputed_ari_nstart1 = ari_1, recomputed_ari_nstart25 = ari_25
  ))
}

print(results, digits = 4)
