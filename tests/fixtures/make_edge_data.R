# Edge-case datasets (EDGE_CASE): nobp_step (no blood pressure anywhere, single-step tilt protocol),
# noangle (no tilt angle: events from heart rate), noevent (stationary PRCP recordings, nothing to detect).
# Fantasia never has respiration or blood pressure in these cases.
suppressPackageStartupMessages(library(data.table))
source("config.R"); source("R/load.R"); load_project_code(c("core", "model"))
root <- Sys.getenv("SPARK_TEST_ROOT"); stopifnot(nzchar(root)); unlink(file.path(root, "data"), recursive = TRUE)
case <- Sys.getenv("EDGE_CASE", "nobp_step")            # nobp_step | noangle | nofantasia_cov
for (d in c("prcp", "fantasia", "tmst")) dir.create(file.path(root, "data", d), recursive = TRUE)
for (i in 1:3) {
  x <- if (case == "noevent") sim_rest(40 + i, "young", 420) else sim_tilt_angle(40 + i, 420, "step")   # noevent: stationary recordings, no angle
  write_prcp_csv(x, sprintf("prcp%02d", i), file.path(root, "data", "prcp"), has_bp = FALSE, has_angle = !case %in% c("noangle", "noevent"), seed = i)
}
for (g in c("young", "older")) for (i in 1:3) {
  x <- sim_rest(60 + i + (g == "older") * 10, g, 360, m_p = rsa_drive(0.25, 0.04))
  write_fantasia_csv(x, sprintf("f1%s%02d", substr(g, 1, 1), i), file.path(root, "data", "fantasia"), has_bp = FALSE, has_resp = FALSE, seed = i)
}
for (i in 1:3) write_tmst_csv(sim_exercise(70 + i, 700, h = 0.02), sprintf("ex%02d", i), file.path(root, "data", "tmst"))
