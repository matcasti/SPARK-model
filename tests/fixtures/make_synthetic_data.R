# Writes synthetic CSV recordings (tools/prep-data.R layout) that mimic the three datasets.
source("config.R"); suppressPackageStartupMessages(library(data.table))
source("R/load.R"); load_project_code(c("core", "model"))
root <- Sys.getenv("SPARK_TEST_ROOT", unset = "tests/_run")
unlink(file.path(root, "data"), recursive = TRUE)   # always start from fresh files
tr <- write_synthetic_data(file.path(root, "data"), n_prcp = 4, n_ex = 6, n_fant = 5, seed = 7,
                          prcp_protocols = c("step", "repeated", "repeated"))
saveRDS(tr$truth, file.path(root, "synthetic_truth.rds")); print(tr$truth[, .N, by = Dataset])
