# All paths are relative to the repository root unless absolute paths are supplied.
config <- list(
  data_dir = "data", output_dir = "results", seed = 1324L,
  months = c("April", "May", "June", "July", "August", "September"),
  treatments = c("CKN", "CKR", "T30N", "T30R", "T45N", "T45R"),
  abundance_thresholds = c(fungi = 3055, bacteria = 1848),
  # Copy of the supplied best result; an absolute path can also be used.
  cluster_file = "data/k40.RData", k = 40L, network_points = 100L,
  network_penalty = "legacy", # "legacy" or "adaptive"; see docs/METHODS.md
  fine_modules = "all", # "all" or an integer vector; no module is privileged
  traits = c("CMF", "NMF", "PMF"),
  bootstrap_replicates = 200L, forest_trees = 1000L
)
