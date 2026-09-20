# Required clustering stage: load the study-selected 40-module result.
# All downstream analyses must use these validated memberships.
# Clustering is not refitted; original module IDs are preserved.
load_modules <- function(filename, feature_ids, expected_k = 40L) {
  environment <- new.env(parent = emptyenv())
  load(filename, envir = environment)
  if (!exists("return_object", envir = environment, inherits = FALSE))
    stop("Expected return_object in the supplied RData file")
  result <- environment$return_object
  d <- result$clustered_data
  if (is.null(dim(d)) || is.null(rownames(d)) || anyDuplicated(rownames(d)))
    stop("clustered_data must have unique OTU row names")
  if (!"cluster" %in% colnames(d)) stop("Missing named cluster column")
  module <- d[, "cluster"]
  if (any(!is.finite(module)) || any(module != as.integer(module))) stop("Invalid module IDs")
  module <- setNames(as.integer(module), rownames(d))
  if (!setequal(unique(module), seq_len(expected_k))) stop("Expected modules 1 through ", expected_k)
  if (!setequal(feature_ids, names(module)))
    stop("OTU IDs differ between preprocessing and k40.RData; do not align by position")
  list(module = module[feature_ids], k = expected_k,
       source = normalizePath(filename, winslash = "/"),
       source_md5 = unname(tools::md5sum(filename)), legacy_BIC = result$BIC,
       legacy_LL = result$LL)
}
