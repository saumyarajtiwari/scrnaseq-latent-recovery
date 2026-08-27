manifest <- read.csv("data/processed/embedding_manifest.csv")
master   <- read.csv("data/processed/step4_master_results_table.csv")

test_master <- master[!is.na(master$ari), ]
set.seed(1)
sim_sample  <- test_master[test_master$data_type == "simulated", ][sample(1:sum(test_master$data_type=="simulated"), 5), ]
real_sample <- test_master[test_master$data_type == "real", ][sample(1:sum(test_master$data_type=="real"), 5), ]
test_rows <- rbind(sim_sample, real_sample)

for (i in seq_len(nrow(test_rows))) {
  r <- test_rows[i, ]
  if (r$data_type == "simulated") {
    mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                      manifest$run_id == r$run_id & manifest$is_null_control == FALSE, ]
  } else {
    mrow <- manifest[manifest$source == r$source & manifest$method == r$method &
                      manifest$data_type == "real", ]
  }
  if (nrow(mrow) == 0) { cat("row", i, ": NO MATCH\n"); next }

  obj <- readRDS(mrow$file_path[1])
  tg  <- obj$true_group
  emb <- obj$embedding

  cat(sprintf("\n--- row %d: source=%s method=%s run_id=%s ---\n",
              i, r$source, r$method, as.character(r$run_id)))
  cat("class(true_group):", class(tg), "\n")
  cat("length:", length(tg), " | n real NA:", sum(is.na(tg)),
      " | n literal 'NA' string:", sum(tg == "NA", na.rm = TRUE), "\n")
  cat("unique values (up to 20): ", paste(head(unique(tg), 20), collapse = " | "), "\n")
  cat("any NA in embedding matrix:", any(is.na(emb)), "\n")
  cat("dim(embedding):", paste(dim(emb), collapse=" x "), "\n")
}
