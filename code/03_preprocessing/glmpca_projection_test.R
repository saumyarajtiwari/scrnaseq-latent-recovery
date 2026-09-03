suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix); library(glmpca)
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
test_idx <- setdiff(seq_len(n_cells), train_idx)[1:20]

train_mat <- as.matrix(counts_mat[, train_idx])
set.seed(42)
fit <- glmpca(train_mat, L=30, fam="nb", ctl=list(verbose=FALSE, maxIter=25, minIter=2, lr=0.0008))

V <- as.matrix(fit$loadings)
coefX_vec <- as.numeric(as.matrix(fit$coefX)[,1])
train_offset_mean <- mean(log(Matrix::colSums(train_mat)))
gf <- fit$glmpca_family

project_one_cell <- function(y_vec, lib_size, V, coefX_vec, offset_mean, gf, L) {
  offset_i <- log(lib_size) - offset_mean
  deviance_fn <- function(u) {
    eta <- offset_i + coefX_vec + as.numeric(V %*% u)
    gf$dev_func(y_vec, eta)
  }
  grad_fn <- function(u) {
    eta <- offset_i + coefX_vec + as.numeric(V %*% u)
    ig <- gf$infograd(y_vec, eta)
    as.numeric(-2 * t(V) %*% ig$grad)
  }
  opt <- optim(par = rep(0, L), fn = deviance_fn, gr = grad_fn, method = "BFGS", control = list(maxit = 100))
  opt$par
}

test_mat <- as.matrix(counts_mat[, test_idx])
test_lib_sizes <- Matrix::colSums(test_mat)
projected_U <- t(sapply(seq_len(ncol(test_mat)), function(i) {
  project_one_cell(test_mat[, i], test_lib_sizes[i], V, coefX_vec, train_offset_mean, gf, 30)
}))
rownames(projected_U) <- colnames(test_mat)

# --- Compare directly against the trusted original full-scale fit ---
orig <- readRDS("data/processed/pca_glmpca/real/baron/baron_harmonized_glmpca.rds")
test_cell_names <- colnames(test_mat)
orig_sub <- orig$embedding[test_cell_names, ]

cat("Original factors range:", range(orig_sub), "\n")
cat("Projected U range:", range(projected_U), "\n")
cat("\nOriginal factors, cell 1, first 5 components:\n")
print(round(orig_sub[1, 1:5], 4))
cat("\nProjected U, cell 1, first 5 components:\n")
print(round(projected_U[1, 1:5], 4))
cat("\nRatio (projected/original), cell 1, first 5:\n")
print(round(projected_U[1, 1:5] / orig_sub[1, 1:5], 2))

saveRDS(list(V=V, coefX=coefX_vec, offset_mean=train_offset_mean, projected_U=projected_U,
             orig_sub=orig_sub, test_idx=test_idx, train_idx=train_idx),
        "results/glmpca_projection_test_result.rds")
