suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(glmpca); library(cluster)
})

sce <- readRDS("data/real/baron/baron_harmonized.rds")
counts_mat <- counts(sce)
gene_totals <- Matrix::rowSums(counts_mat)
counts_mat <- counts_mat[gene_totals > 0, ]

MAX_GENES <- 3000
if (nrow(counts_mat) > MAX_GENES) {
  gene_means <- Matrix::rowMeans(counts_mat)
  top_genes <- order(gene_means, decreasing=TRUE)[1:MAX_GENES]
  counts_mat <- counts_mat[top_genes, ]
}

set.seed(42)
n_cells <- ncol(counts_mat)
train_idx <- sample(seq_len(n_cells), min(2000, floor(n_cells * 0.7)))
test_idx <- setdiff(seq_len(n_cells), train_idx)

train_mat <- as.matrix(counts_mat[, train_idx])
set.seed(42)
fit <- glmpca(train_mat, L=30, fam="nb", ctl=list(verbose=FALSE, maxIter=25, minIter=2, lr=0.0008))
V <- as.matrix(fit$loadings)
coefX_vec <- as.numeric(as.matrix(fit$coefX)[,1])
train_offset_mean <- mean(log(Matrix::colSums(train_mat)))
gf <- fit$glmpca_family

project_one_cell <- function(y_vec, lib_size, V, coefX_vec, offset_mean, gf, L) {
  offset_i <- log(lib_size) - offset_mean
  deviance_fn <- function(u) gf$dev_func(y_vec, offset_i + coefX_vec + as.numeric(V %*% u))
  grad_fn <- function(u) {
    eta <- offset_i + coefX_vec + as.numeric(V %*% u)
    ig <- gf$infograd(y_vec, eta)
    as.numeric(-2 * t(V) %*% ig$grad)
  }
  optim(par = rep(0, L), fn = deviance_fn, gr = grad_fn, method = "BFGS", control = list(maxit = 100))$par
}

test_mat <- as.matrix(counts_mat[, test_idx])
test_lib_sizes <- Matrix::colSums(test_mat)
cat("Projecting", ncol(test_mat), "held-out test cells...\n")
projected_U <- t(sapply(seq_len(ncol(test_mat)), function(i) {
  project_one_cell(test_mat[, i], test_lib_sizes[i], V, coefX_vec, train_offset_mean, gf, 30)
}))
rownames(projected_U) <- colnames(test_mat)

orig <- readRDS("data/processed/pca_glmpca/real/baron/baron_harmonized_glmpca.rds")
test_cell_names <- colnames(test_mat)
orig_sub <- orig$embedding[test_cell_names, ]
true_labels <- sce$true_group[test_idx]

evaluate_embedding <- function(emb, labels, label_name) {
  k <- length(unique(labels))
  set.seed(42)
  emb <- scale(emb)
  km <- kmeans(emb, centers = k, nstart = 25)
  ari <- mclust::adjustedRandIndex(km$cluster, labels)
  sil <- mean(silhouette(km$cluster, dist(emb))[, 3])
  cat(label_name, "- ARI:", round(ari, 3), " Silhouette:", round(sil, 3), "\n")
}

if (!requireNamespace("mclust", quietly=TRUE)) install.packages("mclust", repos="https://cloud.r-project.org")

cat("\n=== Downstream metric comparison (same held-out cells) ===\n")
evaluate_embedding(projected_U, true_labels, "Projected (chunked)")
evaluate_embedding(orig_sub, true_labels, "Original (full-scale fit)")

saveRDS(list(projected_U=projected_U, orig_sub=orig_sub, true_labels=true_labels, test_idx=test_idx),
        "results/glmpca_projection_downstream_result.rds")
