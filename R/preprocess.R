read_inputs <- function(config) {
  read_one <- function(name) {
    d <- openxlsx::read.xlsx(file.path(config$data_dir, name), check.names = FALSE)
    required <- c("sample", "Treatment", "Month")
    if (!all(required %in% names(d))) stop(name, ": missing metadata columns")
    d$sample <- as.character(d$sample)
    if (anyNA(d$sample) || anyDuplicated(d$sample)) stop(name, ": invalid sample IDs")
    d
  }
  f <- read_one("fungi_otu.xlsx")
  b <- read_one("bacteria_out.xlsx")
  s <- read_one("soil_properties.xlsx")
  align <- function(d) {
    if (!setequal(f$sample, d$sample)) stop("Sample IDs differ across input tables")
    d <- d[match(f$sample, d$sample), , drop = FALSE]
    if (!identical(as.character(f$Treatment), as.character(d$Treatment)) ||
        !identical(as.character(f$Month), as.character(d$Month))) stop("Metadata mismatch")
    d
  }
  b <- align(b); s <- align(s)
  keep <- f$Month %in% config$months & f$Treatment %in% config$treatments
  if (!any(keep)) stop("No samples selected")
  list(fungi = f[keep, ], bacteria = b[keep, ], soil = s[keep, ])
}

prepare_data <- function(config) {
  tables <- read_inputs(config)
  transform_one <- function(d, kingdom) {
    cols <- grep("^OTU", names(d), value = TRUE)
    if (!length(cols)) stop("No OTU columns found")
    x <- as.matrix(d[, cols, drop = FALSE]); storage.mode(x) <- "double"
    if (any(!is.finite(x)) || any(x < 0)) stop("Counts must be finite and nonnegative")
    x <- x[, colSums(x) >= config$abundance_thresholds[[kingdom]], drop = FALSE]
    if (ncol(x) < 2L || any(rowSums(x) == 0)) stop("Insufficient counts after filtering")
    rownames(x) <- d$sample
    # Prevent silent deletion of sparse samples or features.
    if (any(x == 0)) x <- zCompositions::cmultRepl(x, method = "GBM", output = "p-counts", z.delete = FALSE)
    x <- as.matrix(x)
    if (any(!is.finite(x)) || any(x <= 0)) stop("Zero replacement failed")
    x <- sweep(log(x), 1L, rowMeans(log(x)), "-")
    colnames(x) <- paste0(kingdom, "_", colnames(x))
    x
  }
  clr <- cbind(transform_one(tables$fungi, "fungi"), transform_one(tables$bacteria, "bacteria"))
  varying <- apply(clr, 2, sd) > 1e-12
  dropped <- colnames(clr)[!varying]
  x <- scale(clr[, varying, drop = FALSE])
  x <- x - min(x)
  groups <- setNames(lapply(config$treatments, function(tr) {
    idx <- which(tables$soil$Treatment == tr)
    if (length(idx) < 5) stop("Each treatment needs at least five samples")
    idx <- idx[order(rowSums(x[idx, , drop = FALSE]), tables$soil$sample[idx])]
    h <- rowSums(x[idx, , drop = FALSE])
    if (min(h) <= 0 || diff(range(h)) <= 0) stop("Invalid habitat index")
    list(x = h, abundance = x[idx, , drop = FALSE], soil = tables$soil[idx, , drop = FALSE])
  }), config$treatments)
  list(groups = groups, feature_ids = colnames(x), dropped_constant_features = dropped,
       samples = tables$soil, config = config)
}

# Fit y = a * (x / geometric_mean(x))^b without random perturbations.
fit_power <- function(x, y, new_x = x) {
  if (any(!is.finite(c(x, y, new_x))) || any(x <= 0) || any(new_x <= 0) || any(y < 0))
    stop("Power fitting requires positive x and nonnegative y")
  if (max(y) == 0) return(rep(0, length(new_x)))
  center <- exp(mean(log(x))); lx <- log(x / center)
  objective <- function(b) {
    z <- exp(b * lx); a <- sum(y * z) / sum(z * z)
    sum((y - a * z)^2)
  }
  opt <- optimize(objective, c(-20, 20))
  if (abs(opt$minimum) > 19.99) warning("Power exponent reached its bound")
  z <- exp(opt$minimum * lx); a <- sum(y * z) / sum(z * z)
  a * exp(opt$minimum * log(new_x / center))
}
