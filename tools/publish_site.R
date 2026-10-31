# ==============================================================================
# publish_site.R -- put the rendered manuscript on GitHub Pages (make publish).
#
#   Rscript tools/publish_site.R             assemble the site and force-push it to the `site` branch of `origin`
#   Rscript tools/publish_site.R --dry-run   assemble it in a temporary folder and list it; nothing is pushed
#
# Reads   manuscript/_site/   (make render)   and   .github/workflows/pages.yml
# Writes  branch `site` of origin: public/ (the pages) + the workflow that deploys it. One commit, replaced on every publish,
#         so the repository does not grow with every version of the PDFs and HTML. The workflow publishes it
#         (https://<owner>.github.io/<repo>/).
# SITE_REMOTE=<url or path> pushes to another repository instead of origin (to try it out).
# Refuses to publish a site older than the manuscript sources or the results it was rendered from.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
unknown <- setdiff(args, "--dry-run")
if (length(unknown)) stop("Unknown argument: ", paste(unknown, collapse = ", "), call. = FALSE)

site_dir <- "manuscript/_site"
pages <- c("index.html", "manuscript.html", "supplementary.html", "manuscript.pdf", "supplementary.pdf", "model_code.R", "search.json", "site_libs")
workflow <- ".github/workflows/pages.yml"

missing <- c(file.path(site_dir, pages), workflow)[!file.exists(c(file.path(site_dir, pages), workflow))]
if (length(missing)) stop("Missing: ", paste(missing, collapse = ", "), "\nRun `make render` first.", call. = FALSE)

# The site must be at least as new as everything it is built from
built <- file.mtime(file.path(site_dir, "manuscript.html"))
sources <- c(list.files("manuscript", pattern = "\\.(qmd|yml)$", full.names = TRUE), list.files("manuscript/support", full.names = TRUE),
             list.files("manuscript/templates", full.names = TRUE), list.files("results", pattern = "\\.rds$", full.names = TRUE),
             list.files("R/report", full.names = TRUE))
newer <- sources[file.mtime(sources) > built]
if (length(newer)) stop("The rendered site is older than: ", paste(head(newer, 5), collapse = ", "), if (length(newer) > 5) ", ...", "\nRun `make render` first.", call. = FALSE)

run_git <- function(dir, ...) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(dir), shQuote(c(...))), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"), call. = FALSE)
  out
}

# ---- Assemble -----------------------------------------------------------------------------------------------------------
tmp <- tempfile("spark-site-")
public <- file.path(tmp, "public")
dir.create(public, recursive = TRUE)
dir.create(file.path(tmp, dirname(workflow)), recursive = TRUE)
invisible(file.copy(file.path(site_dir, pages), public, recursive = TRUE))
invisible(file.copy(workflow, file.path(tmp, workflow)))

files <- list.files(public, recursive = TRUE, all.files = TRUE)
size_mb <- sum(file.info(file.path(public, files))$size) / 1e6
cat(sprintf("Site: %d files, %.1f MB (%s)\n", length(files), size_mb, public))
page_mb <- function(page) {
  path <- file.path(public, page)
  all_files <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE) else path
  sum(file.info(all_files)$size) / 1e6
}
cat(sprintf("  %-20s %6.1f MB\n", pages, vapply(pages, page_mb, numeric(1))), sep = "")

if (dry_run) {
  cat("Dry run: nothing pushed. Inspect the folder above (open public/index.html).\n")
  quit(save = "no", status = 0)
}

# ---- Commit and push ------------------------------------------------------------------------------------------------------
remote <- Sys.getenv("SITE_REMOTE", run_git(".", "remote", "get-url", "origin")[1])      # SITE_REMOTE: another repository (tests)
sha <- run_git(".", "rev-parse", "--short", "HEAD")[1]
user <- tryCatch(run_git(".", "config", "user.name")[1], error = function(e) "SPARK site")
email <- tryCatch(run_git(".", "config", "user.email")[1], error = function(e) "noreply@users.noreply.github.com")
dirty <- length(run_git(".", "status", "--porcelain", "--", "manuscript", "R", "analysis", "engine", "config.R", "pipeline.R")) > 0
if (dirty) cat("Note: the working tree has uncommitted changes; the site shows them, the commit message names the last commit only.\n")

invisible(run_git(tmp, "init", "-q", "-b", "site"))
invisible(run_git(tmp, "add", "-A"))
invisible(run_git(tmp, "-c", paste0("user.name=", user), "-c", paste0("user.email=", email), "commit", "-q", "-m", sprintf("Site rendered from %s (%s)", sha, format(Sys.time(), "%Y-%m-%d %H:%M"))))
cat(sprintf("Pushing to %s, branch `site` (replaced, not appended)\n", remote))
invisible(run_git(tmp, "push", "--force", remote, "site:site"))

if (grepl("github\\.com", remote)) {
  slug <- sub("\\.git$", "", sub("^.*github\\.com[:/]", "", remote))
  owner_repo <- strsplit(slug, "/", fixed = TRUE)[[1]]
  cat(sprintf("Pushed. The workflow deploys it: https://github.com/%s/actions\nSite: https://%s.github.io/%s/\n", slug, owner_repo[1], owner_repo[2]))
} else cat("Pushed.\n")
