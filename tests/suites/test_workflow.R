# Checks of the project's structure and bookkeeping, not of any number: that what pipeline.R declares is what the scripts do,
# that the code is laid out as R/load.R says, and that the documentation and deposited code have not drifted.
#   Rscript tests/suites/test_workflow.R      (seconds; no data needed)
source("config.R"); source("pipeline.R"); source("R/load.R")
n_pass <- 0L
check <- function(name, cond, detail = NULL) {
  if (!isTRUE(cond)) stop("FAILED: ", name, if (length(detail)) paste0("\n    ", paste(detail, collapse = "\n    ")), call. = FALSE)
  n_pass <<- n_pass + 1L; cat(sprintf("ok   %s\n", name))
}
read_text <- function(file) paste(readLines(file, warn = FALSE), collapse = "\n")

# ---- 1. Layout of R/ -------------------------------------------------------------------------------------------------------
r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE)
check("R/ has exactly the layers of R/load.R (+ load.R itself)", setequal(unique(dirname(r_files)), c(PROJECT_LAYERS, ".")) && identical(r_files[dirname(r_files) == "."], "load.R"),
      setdiff(unique(dirname(r_files)), c(PROJECT_LAYERS, ".")))
check("every R file starts with a header naming itself", all(vapply(setdiff(r_files, "load.R"), function(f) grepl(paste0("^# R/", f, " --"), readLines(file.path("R", f), n = 2)[2]), NA)),
      setdiff(r_files, "load.R")[!vapply(setdiff(r_files, "load.R"), function(f) grepl(paste0("^# R/", f, " --"), readLines(file.path("R", f), n = 2)[2]), NA)])

defined <- do.call(rbind, lapply(r_files, function(f) {
  expressions <- parse(file.path("R", f), keep.source = FALSE)
  names <- vapply(expressions, function(e) if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") && is.name(e[[2]])) as.character(e[[2]]) else NA_character_, "")
  data.frame(file = f, name = names[!is.na(names)])
}))
duplicated_names <- unique(defined$name[duplicated(defined$name)])
check("no function or constant is defined twice in R/", !length(duplicated_names),
      vapply(duplicated_names, function(n) sprintf("%s: %s", n, paste(defined$file[defined$name == n], collapse = ", ")), ""))

# ---- 2. pipeline.R against the scripts -------------------------------------------------------------------------------------
by_number <- function(x) x[order(basename(x))]
scripts <- by_number(setdiff(list.files("analysis", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE, recursive = TRUE), "analysis/00_setup.R"))
registered <- vapply(PIPELINE, `[[`, "", "script")
check("pipeline.R lists every numbered script of analysis/ (except 00_setup.R) once, in order", identical(by_number(registered), scripts) && identical(registered, by_number(registered)) && !anyDuplicated(pipeline_ids()),
      c(setdiff(scripts, registered), setdiff(registered, scripts)))
check("step ids are the numbers of the script names", all(vapply(PIPELINE, function(s) startsWith(basename(s$script), paste0(s$id, "_")), NA)))

# The functions in R/steps/ that a script calls are named in their headers ("Called by analysis/<group>/NN_*.R")
step_files <- list.files("R/steps", full.names = TRUE)
called_by <- lapply(step_files, function(f) regmatches(paste(readLines(f, n = 12), collapse = " "), gregexpr("analysis/[a-z]+/[0-9]{2}_[A-Za-z_]+\\.R", paste(readLines(f, n = 12), collapse = " ")))[[1]])
names(called_by) <- step_files
check("every file of R/steps/ names the script(s) that call it", all(lengths(called_by) > 0), names(called_by)[lengths(called_by) == 0])

datasets <- c("prcp", "fantasia", "exercise")
# Names given to read_result() / write_result() in a piece of code, found in the parsed code (not by pattern). A name built as
# paste0("fit_", x) stands for the three fit_<dataset> results; a name held in a variable cannot be known and is ignored.
calls_of <- function(expr, fn) {
  found <- list()
  walk <- function(e) if (is.call(e)) {
    if (identical(e[[1]], as.name(fn))) found[[length(found) + 1L]] <<- e
    for (i in seq_along(e)) if (!identical(e[[i]], quote(expr = ))) walk(e[[i]])
  }
  for (e in expr) walk(e)
  found
}
result_names <- function(files, fn, argument) {
  names <- unlist(lapply(files, function(f) lapply(calls_of(parse(f, keep.source = FALSE), fn), function(call) {
    value <- call[[argument + 1L]]
    if (is.character(value)) value else if (is.call(value) && identical(value[[1]], as.name("paste0")) && identical(value[[2]], "fit_")) paste0("fit_", datasets) else character()
  })))
  unique(names)
}
for (step in PIPELINE) {
  if (identical(step$reads, "*")) next                                  # the figures read whichever results exist
  own <- c(step$script, step_files[vapply(called_by, function(x) identical(x, step$script), NA)])          # files only this script uses
  shared <- step_files[vapply(called_by, function(x) step$script %in% x && length(x) > 1L, NA)]           # R/steps files called by several scripts
  files <- c(own, shared)
  writes <- result_names(files, "write_result", 2); reads <- setdiff(result_names(own, "read_result", 1), writes)   # a step may read back what it has just written
  declared <- pipeline_reads(step)
  check(sprintf("step %s: writes in the code = writes in pipeline.R", step$id), setequal(writes, step$writes), c(sprintf("code: %s", paste(writes, collapse = ", ")), sprintf("pipeline.R: %s", paste(step$writes, collapse = ", "))))
  check(sprintf("step %s: reads in the code = reads in pipeline.R", step$id), all(reads %in% declared) && all(declared %in% c(reads, setdiff(result_names(shared, "read_result", 1), writes))),
        c(sprintf("code: %s", paste(reads, collapse = ", ")), sprintf("pipeline.R: %s", paste(declared, collapse = ", "))))
  check(sprintf("step %s: its optional reads are read with required = FALSE", step$id),
        all(vapply(pipeline_reads(step, "optional"), function(r) !grepl(sprintf('read_result\\("%s", cfg\\)', r), paste(vapply(files, read_text, ""), collapse = "\n")), NA)))
}
check("a step reads only results that an earlier step writes", all(vapply(seq_along(PIPELINE), function(i)
  all(pipeline_reads(PIPELINE[[i]]) %in% unlist(lapply(PIPELINE[seq_len(i - 1L)], `[[`, "writes"))), NA)))

# ---- 3. Documentation and deposited code -----------------------------------------------------------------------------------
check("docs/DATAFLOW.md is what tools/make_docs.R generates", system2("Rscript", c("tools/make_docs.R", "--check"), stdout = FALSE, stderr = FALSE) == 0,
      "run: Rscript tools/make_docs.R")
referenced <- unique(unlist(lapply(c("README.md", "docs/DATAFLOW.md"), function(f) regmatches(read_text(f), gregexpr("`((R|analysis|tools|tests|engine|manuscript)/[A-Za-z0-9_./-]+\\.(R|qmd|cpp|md))`", read_text(f)))[[1]])))
referenced <- gsub("`", "", referenced)
check("every source file named in README.md and docs/DATAFLOW.md exists", all(file.exists(referenced)), referenced[!file.exists(referenced)])
check("the engine is byte-identical to the recorded one", unname(tools::md5sum(cfg$paths$engine)) == trimws(readLines("tests/fixtures/engine.md5", n = 1)))
check("manuscript/model_code.R (deposited code) is the engine", identical(unname(tools::md5sum("manuscript/model_code.R")), unname(tools::md5sum(cfg$paths$engine))))
supplement_scripts <- regmatches(read_text("manuscript/supplementary.qmd"), gregexpr("(?<=\\| `)[0-9]{2}_[A-Za-z_]+\\.R(?=` \\|)", read_text("manuscript/supplementary.qmd"), perl = TRUE))[[1]]
check("the table of scripts in the supplement lists the scripts of pipeline.R", identical(supplement_scripts, basename(registered)), setdiff(basename(registered), supplement_scripts))
check("the manuscript and the supplement source files that exist", all(file.exists(unique(unlist(lapply(c("manuscript/manuscript.qmd", "manuscript/supplementary.qmd"), function(f)
  file.path("manuscript", regmatches(read_text(f), gregexpr('(?<=source\\(")[^"]+', read_text(f), perl = TRUE))[[1]])))))))

# ---- 4. Journal limits that do not depend on the results (Chaos, Solitons & Fractals guide for authors) --------------------
highlights <- sub("^- ", "", grep("^- ", readLines("manuscript/support/highlights.md", warn = FALSE), value = TRUE))
check("highlights: 3 to 5 bullets of at most 85 characters", length(highlights) %in% 3:5 && all(nchar(highlights) <= 85), highlights[nchar(highlights) > 85])
keywords <- strsplit(sub("^\\*\\*Keywords\\*\\*: ", "", grep("^\\*\\*Keywords\\*\\*:", readLines("manuscript/manuscript.qmd", warn = FALSE), value = TRUE)), ", ")[[1]]
check("keywords: 1 to 7", length(keywords) %in% 1:7, keywords)

cat(sprintf("\n%d workflow checks passed.\n", n_pass))
