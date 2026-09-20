scale01 <- function(x) {
  if (diff(range(x)) == 0) rep(0, length(x)) else (x - min(x)) / diff(range(x))
}

trait_importance <- function(x, y, n_boot = 200L, trees = 1000L, seed = 1324L) {
  keep <- is.finite(y) & apply(x, 1, function(z) all(is.finite(z)))
  x <- x[keep, , drop = FALSE]; y <- y[keep]
  if (nrow(x) < 10 || sd(y) == 0) stop("Trait requires at least 10 complete, varying observations")
  if (n_boot < 2) stop("Use at least two bootstrap replicates")
  original_names <- colnames(x)
  # Standardized coefficients are comparable across nodes; constant nodes remain zero.
  varying <- apply(x, 2, sd) > 1e-12
  en <- rf <- matrix(0, n_boot, ncol(x), dimnames = list(NULL, original_names))
  if (!any(varying)) stop("All node predictors are constant")
  z <- scale(x[, varying, drop = FALSE])
  set.seed(seed)
  for (b in seq_len(n_boot)) {
    for (attempt in seq_len(100)) {
      idx <- sample.int(nrow(z), replace = TRUE)
      if (length(unique(idx)) >= 5 && sd(y[idx]) > 0) break
    }
    if (length(unique(idx)) < 5 || sd(y[idx]) == 0) stop("Degenerate bootstrap draw")
    ids <- unique(idx)
    # Copies of one observation must not occur in both training and validation.
    fold_map <- sample(rep(seq_len(5), length.out = length(ids)))
    folds <- fold_map[match(idx, ids)]
    design <- z[idx, , drop = FALSE]
    if (ncol(design) == 1) design <- cbind(design, padding = 0)
    fit <- glmnet::cv.glmnet(design, y[idx], alpha = 0.5, foldid = folds)
    en[b, varying] <- as.matrix(coef(fit, s = "lambda.min"))[-1, 1][seq_len(sum(varying))]
    # Sample original rows once per tree inside ranger. This avoids duplicate
    # original samples leaking into nominal OOB permutation evaluations.
    forest <- ranger::ranger(x = as.data.frame(z), y = y, num.trees = trees,
                             importance = "permutation", min.node.size = 3,
                             mtry = max(1, floor(sqrt(ncol(z)))), replace = FALSE,
                             sample.fraction = 0.632, seed = seed + b, num.threads = 1)
    rf[b, varying] <- forest$variable.importance[colnames(z)]
  }
  frequency <- colMeans(en != 0); en_score <- frequency * colMeans(abs(en))
  rf_mean <- colMeans(rf)
  result <- data.frame(node = original_names, selection_frequency = frequency,
    mean_abs_coefficient = colMeans(abs(en)), mean_coefficient = colMeans(en),
    elastic_net_score = en_score, rf_mean = rf_mean, rf_sd = apply(rf, 2, sd),
    consensus = (scale01(en_score) + scale01(rf_mean)) / 2)
  result$rank <- rank(-result$consensus, ties.method = "min")
  result$spearman_rho <- vapply(seq_len(ncol(x)), function(j)
    if (sd(x[, j]) == 0) NA_real_ else cor(x[, j], y, method = "spearman"), numeric(1))
  result$p_value <- vapply(seq_len(ncol(x)), function(j)
    if (sd(x[, j]) == 0) NA_real_ else cor.test(x[, j], y, method = "spearman", exact = FALSE)$p.value, numeric(1))
  result$q_value <- p.adjust(result$p_value, method = "BH")
  list(summary = result[order(result$rank), ], elastic_net_replicates = en,
       forest_replicates = rf, n_complete = length(y))
}

save_importance <- function(x, soil, config, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (trait in config$traits) {
    if (!trait %in% names(soil)) stop("Missing soil trait: ", trait)
    result <- trait_importance(x, soil[[trait]], config$bootstrap_replicates,
                               config$forest_trees, config$seed)
    write.csv(result$summary, file.path(directory, paste0(trait, ".csv")), row.names = FALSE)
    saveRDS(result, file.path(directory, paste0(trait, ".rds")))
    pdf(file.path(directory, paste0(trait, ".pdf")), width = 10, height = 5)
    par(mar = c(8, 4, 3, 1))
    barplot(result$summary$consensus, names.arg = result$summary$node,
            las = 2, col = "steelblue", main = trait, ylab = "Consensus importance")
    dev.off()
  }
}
