# ==============================================================================
# R/core/utils.R -- small helpers used everywhere.
#   logging           log_msg()
#   results store     write_result(), read_result(), ensure_dir()
#   reproducibility   id_seed(), with_seed()
#   parallel map      map_recordings(), is_err()
#   summaries         median_range(), format_median_range()
#   subsetting        pick_subset()
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })


# ---- Logging ----------------------------------------------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}


# ---- Results store (results/<name>.rds) -----------------------------------------
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

write_result <- function(object, name, cfg) {
  file <- file.path(ensure_dir(cfg$paths$results), paste0(name, ".rds"))
  saveRDS(object, file)
  invisible(object)
}

#' Read results/<name>.rds; stops with a hint if missing unless `required = FALSE`.
read_result <- function(name, cfg, required = TRUE) {
  file <- file.path(cfg$paths$results, paste0(name, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  if (required) {
    stop(sprintf("Missing results/%s.rds: run the script that creates it first.", name), call. = FALSE)
  }
  NULL
}


# ---- Reproducibility ------------------------------------------------------------
#' Stable integer seed for each element of `id`, combined with a base seed. A polynomial hash of the characters modulo the prime
#' 2^31 - 1 (every intermediate value is exact in double precision); the earlier position-weighted sum of character codes gave only
#' about 37 distinct seeds for 100 consecutive indices, so simulated "independent" recordings were repeated.
id_seed <- function(id, base = 1L) {
  modulus <- 2147483647
  seed_of <- function(x) {
    hash <- (as.numeric(base) * 7919) %% modulus
    for (code in utf8ToInt(x)) hash <- (hash * 131 + code) %% modulus
    as.integer(hash)
  }
  vapply(as.character(id), seed_of, integer(1), USE.NAMES = FALSE)
}

#' Evaluate `expr` under a temporary seed and restore the RNG state afterwards.
with_seed <- function(seed, expr) {
  previous <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
  on.exit(if (!is.null(previous)) assign(".Random.seed", previous, globalenv()), add = TRUE)
  set.seed(seed)
  force(expr)
}


# ---- Parallel map over recordings ---------------------------------------------------
#' lapply over `X`, in parallel when `cores > 1` (not on Windows). An error in one
#' item does not stop the others: that item becomes a "spark_error" string (see is_err()).
map_recordings <- function(X, FUN, cores = 1L, ...) {
  safe_FUN <- function(x) {
    tryCatch(FUN(x, ...), error = function(e) structure(conditionMessage(e), class = "spark_error"))
  }
  if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(X, safe_FUN, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(X, safe_FUN)
  }
}

is_err <- function(x) inherits(x, "spark_error") || inherits(x, "try-error")


# ---- Summaries across recordings --------------------------------------------------------
#' Median and 2.5-97.5 percentile range across recordings (between-recording dispersion).
median_range <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
  q <- quantile(x, c(probs[1], 0.5, probs[2]), names = FALSE)
  c(med = q[2], lo = q[1], hi = q[3])
}

#' "median [lo, hi]" as text, e.g. "0.42 ms [0.30, 0.65]".
format_median_range <- function(x, digits = 3, unit = "") {
  r <- median_range(x)
  if (anyNA(r)) return("NA")
  template <- paste0("%.", digits, "f%s [%.", digits, "f, %.", digits, "f]")
  sprintf(template, r["med"], unit, r["lo"], r["hi"])
}


# ---- Subsetting ------------------------------------------------------------------------------
#' Deterministic subset of at most `n` recordings for compute-heavy analyses,
#' balanced across groups when `index` has a Group column.
pick_subset <- function(index, n, seed = 1L) {
  ids <- index$SubjectID
  if (length(ids) <= n) return(ids)
  with_seed(seed, {
    if ("Group" %in% names(index) && !anyNA(index$Group)) {
      per_group <- ceiling(n / uniqueN(index$Group))
      chosen <- unlist(lapply(split(index$SubjectID, index$Group),
                              function(x) sample(x, min(length(x), per_group))))
      sort(head(chosen, n))
    } else {
      sort(sample(ids, n))
    }
  })
}


# ---- Combining per-recording results ---------------------------------------------------------------
#' Stack the `field` table of every per-recording result; NULL (never an empty
#' zero-column table) when no recording contributes it.
bind_results <- function(results, field) {
  stacked <- rbindlist(lapply(results, `[[`, field), fill = TRUE)
  if (nrow(stacked)) stacked else NULL
}

#' Drop failed items (see map_recordings()) and report them.
drop_failed <- function(results, label = "recording") {
  failed <- vapply(results, is_err, logical(1))
  if (any(failed)) {
    log_msg("failed %s(s): %s", label, paste(as.character(results[failed]), collapse = " | "))
  }
  results[!failed]
}

#' Does a states table have column `name` with enough finite values?
column_available <- function(states, name, min_coverage) {
  name %in% names(states) && mean(is.finite(states[[name]])) >= min_coverage
}

#' x[[field]] or NA when x is NULL (used for step-response fits that may fail).
field_or_na <- function(x, field) if (is.null(x)) NA_real_ else x[[field]]
