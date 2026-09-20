aggregate_modules <- function(x, modules, statistic = c("sum", "mean")) {
  statistic <- match.arg(statistic)
  if (!setequal(names(modules), colnames(x))) stop("OTU IDs do not match module assignments")
  modules <- modules[colnames(x)]
  ids <- sort(unique(modules))
  answer <- vapply(ids, function(m) {
    z <- x[, modules == m, drop = FALSE]
    if (statistic == "sum") rowSums(z) else rowMeans(z)
  }, numeric(nrow(x)))
  colnames(answer) <- paste0("M", ids); answer
}

select_predictors <- function(x, target, penalty = "legacy", seed = 1324L) {
  candidates <- setdiff(seq_len(ncol(x)), target)
  candidates <- candidates[apply(x[, candidates, drop = FALSE], 2, sd) > 1e-12]
  if (!length(candidates) || sd(x[, target]) < 1e-12) return(integer())
  # glmnet requires at least two predictors. Pad a one-predictor problem with a
  # zero column, which cannot be selected, to preserve the same selection rule.
  predictors <- x[, candidates, drop = FALSE]
  if (ncol(predictors) == 1L) predictors <- cbind(predictors, padding = 0)
  set.seed(seed); folds <- sample(rep(seq_len(5), length.out = nrow(x)))
  ridge <- glmnet::cv.glmnet(predictors, x[, target], alpha = 0, foldid = folds)
  beta <- abs(as.matrix(coef(ridge, s = "lambda.1se"))[-1, 1])
  weights <- switch(penalty, legacy = sqrt(pmax(beta, 1e-8)),
                    adaptive = 1 / sqrt(pmax(beta, 1e-8)), stop("Unknown penalty mode"))
  lasso <- glmnet::cv.glmnet(predictors, x[, target], alpha = 1,
                            penalty.factor = weights, foldid = folds)
  beta <- as.matrix(coef(lasso, s = "lambda.1se"))[-1, 1]
  candidates[which(abs(beta[seq_along(candidates)]) > 0)]
}

fit_network <- function(x, habitat, points = 100L, penalty = "legacy", seed = 1324L) {
  if (ncol(x) < 1L || length(habitat) != nrow(x)) stop("Invalid network input")
  grid <- seq(min(habitat), max(habitat), length.out = points)
  smooth <- vapply(seq_len(ncol(x)), function(j) fit_power(habitat, x[, j], grid), numeric(points))
  colnames(smooth) <- colnames(x)
  scaled_grid <- 2 * (grid - min(grid)) / diff(range(grid)) - 1
  # For order-one Legendre derivatives, the original RK expression reduces
  # exactly to a cumulative left-interval integral of each driving curve.
  integral <- apply(smooth, 2, function(z) c(0, cumsum(diff(scaled_grid) * head(z, -1))))
  integral <- matrix(integral, nrow = points, dimnames = dimnames(smooth))
  edges <- data.frame(source = character(), target = character(), weight = numeric(),
                      signed_weight = numeric(), edge_type = integer())
  nodes <- data.frame(node = colnames(x), self_effect = NA_real_, incoming = 0L,
                      fit_rmse = NA_real_, convergence = 0L)
  curves <- vector("list", ncol(x))
  for (j in seq_len(ncol(x))) {
    dep <- select_predictors(smooth, j, penalty, seed + j)
    design <- integral[, c(j, dep), drop = FALSE]; observed <- smooth[, j]
    objective <- function(p) {
      self <- observed[1] + design[, 1] * p[1]
      fitted <- observed[1] + as.vector(design %*% p)
      sum((observed - fitted)^2) + sum(self < 0) * sum(pmax(0, -self))
    }
    opt <- if (ncol(design) == 1L) {
      # Exact least-squares coefficient is admissible for nonnegative power curves.
      p <- sum(design[, 1] * (observed - observed[1])) / max(sum(design[, 1]^2), 1e-20)
      list(par = p, convergence = 0L)
    } else optim(rep(0.001, ncol(design)), objective, method = "Nelder-Mead",
                 control = list(maxit = 50000))
    effects <- sweep(design, 2, opt$par, "*")
    effects[, 1] <- effects[, 1] + observed[1]
    nodes$self_effect[j] <- mean(effects[, 1]); nodes$incoming[j] <- length(dep)
    nodes$fit_rmse[j] <- sqrt(mean((observed - rowSums(effects))^2))
    nodes$convergence[j] <- opt$convergence
    if (opt$convergence != 0) warning("Network optimization did not converge for ", colnames(x)[j])
    if (length(dep)) {
      weight <- colMeans(effects[, -1, drop = FALSE])
      e <- data.frame(source = colnames(x)[dep], target = colnames(x)[j],
                      weight = abs(weight), signed_weight = weight,
                      edge_type = ifelse(weight > 0, 1L, 2L))
      edges <- rbind(edges, e[e$weight > 0, , drop = FALSE])
    }
    curves[[j]] <- list(habitat = grid, observed = observed, effects = effects,
                        fitted = rowSums(effects), coefficients = opt$par)
  }
  names(curves) <- colnames(x)
  list(edges = edges, nodes = nodes, curves = curves, penalty = penalty)
}

save_network <- function(net, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  write.csv(net$edges, file.path(directory, "edges.csv"), row.names = FALSE)
  write.csv(net$nodes, file.path(directory, "nodes.csv"), row.names = FALSE)
  saveRDS(net, file.path(directory, "network.rds"))
  pdf(file.path(directory, "decomposition.pdf"), width = 8, height = 6)
  on.exit(dev.off())
  for (id in names(net$curves)) {
    c <- net$curves[[id]]
    matplot(c$habitat, cbind(c$observed, c$fitted, c$effects), type = "l",
            lty = c(1, 2, rep(3, ncol(c$effects))), xlab = "Habitat index",
            ylab = "Effect", main = id)
    legend("topleft", c("Observed smooth curve", "Total fitted", "Self / incoming effects"),
           lty = 1:3, bty = "n")
  }
}
