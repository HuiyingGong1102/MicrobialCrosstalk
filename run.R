args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "help"
if (stage == "help") {
  cat("Usage: Rscript run.R <prepare|clustering|coarse|fine|importance|fine-importance|all> [config.R]\n")
  quit(status = 0)
}
if (stage == "modules") stage <- "clustering" # Compatibility with the previous entry point.
if (!stage %in% c("prepare", "clustering", "coarse", "fine", "importance", "fine-importance", "all")) stop("Unknown stage")
source(if (length(args) >= 2) args[2] else "config.R")
for (f in c("preprocess", "clustering", "network", "importance")) source(file.path("R", paste0(f, ".R")))
packages <- c("openxlsx", "zCompositions", "glmnet", "ranger")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
path <- function(name) file.path(config$output_dir, name)
set.seed(config$seed)
saveRDS(config, path(paste0("config_", stage, ".rds")))
writeLines(capture.output(sessionInfo()), path(paste0("session_", stage, ".txt")))
if (stage %in% c("prepare", "all")) {
  prepared <- prepare_data(config)
  saveRDS(prepared, path("prepared.rds"))
  write.csv(prepared$samples, path("selected_samples.csv"), row.names = FALSE)
  message("Prepared ", nrow(prepared$samples), " samples and ", length(prepared$feature_ids), " OTUs")
}
if (stage != "prepare") {
  prepared <- readRDS(path("prepared.rds"))
  if (!identical(prepared$config$months, config$months) ||
      !identical(prepared$config$treatments, config$treatments) ||
      !identical(prepared$config$abundance_thresholds, config$abundance_thresholds))
    stop("Preprocessing config changed; rerun prepare")
}
if (stage != "prepare") {
  # Required for every analysis stage, including all: never bypass clustering membership.
  fit <- load_modules(config$cluster_file, prepared$feature_ids, config$k)
  saveRDS(fit, path("modules.rds"))
  write.csv(data.frame(otu = names(fit$module), module = fit$module), path("modules.csv"), row.names = FALSE)
  modules <- fit$module
  fine_ids <- if (identical(config$fine_modules, "all")) sort(unique(modules)) else config$fine_modules
  if (!all(fine_ids %in% modules)) stop("Unknown fine-scale module IDs")
  for (tr in names(prepared$groups)) {
    g <- prepared$groups[[tr]]
    if (stage %in% c("coarse", "all")) {
      x <- aggregate_modules(g$abundance, modules, "sum")
      save_network(fit_network(x, g$x, config$network_points, config$network_penalty, config$seed),
                   path(file.path("networks", "coarse", tr)))
    }
    if (stage %in% c("importance", "all"))
      save_importance(aggregate_modules(g$abundance, modules, "mean"), g$soil, config,
                       path(file.path("importance", "coarse", tr)))
    if (stage %in% c("fine", "fine-importance", "all")) for (m in fine_ids) {
      ids <- names(modules)[modules == m]
      x <- g$abundance[, ids, drop = FALSE]
      if (stage %in% c("fine", "all"))
        save_network(fit_network(x, g$x, config$network_points, config$network_penalty, config$seed),
                     path(file.path("networks", "fine", tr, paste0("M", m))))
      if (stage == "fine-importance")
        save_importance(x, g$soil, config, path(file.path("importance", "fine", tr, paste0("M", m))))
    }
  }
}
message("Completed stage: ", stage)
