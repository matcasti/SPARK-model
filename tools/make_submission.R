# ==============================================================================
# tools/make_submission.R -- gathers the files the journal asks for (the PDF is also what bioRxiv takes) in submission/ (not under version control)
#   make submission          renders the LaTeX route (keeping the .tex; the PDF itself comes from `make render`, i.e. Typst) and then runs this script
# Writes: manuscript.tex and supplementary.tex (single column, compile with pdfLaTeX), Figure_N.pdf / FigureS_N.pdf (vector, in order of appearance),
#         highlights.docx (separate editable file, as the journal asks), manuscript.pdf and supplementary.pdf (Typst); compiles the .tex files as a check.
# ==============================================================================
out <- "submission"; src <- "manuscript"
if (!file.exists(file.path(src, "manuscript.tex")) || !file.exists(file.path(src, "supplementary.tex"))) stop("render the LaTeX route first (`make submission`): manuscript.tex / supplementary.tex not found")
unlink(out, recursive = TRUE); dir.create(file.path(out, "figures"), recursive = TRUE)

# copies <doc>.tex into out/ with its PNG figures swapped for numbered vector PDFs (order of appearance); a figure already
# numbered for another document keeps that file, so nothing is copied twice
used <- character()   # named: figure PNG -> file in submission/
localize <- function(doc, label) {
  lines <- readLines(file.path(src, paste0(doc, ".tex")), encoding = "UTF-8")
  png <- unique(unlist(regmatches(lines, gregexpr("figures/[A-Za-z0-9_]+\\.png", lines))))
  for (p in png) {
    if (is.na(used[p])) {
      pdf <- file.path(src, sub("\\.png$", ".pdf", p))
      if (!file.exists(pdf)) stop("vector version missing: ", pdf)
      used[p] <<- sprintf("figures/%s_%d.pdf", label, sum(startsWith(used, paste0("figures/", label, "_"))) + 1)
      file.copy(pdf, file.path(out, used[p]), overwrite = TRUE)
    }
    lines <- gsub(p, used[p], lines, fixed = TRUE)
  }
  writeLines(lines, file.path(out, paste0(doc, ".tex")), useBytes = TRUE)
}
localize("manuscript", "Figure"); localize("supplementary", "FigureS")

status <- system2("quarto", c("pandoc", file.path(src, "support", "highlights.md"), "-o", file.path(out, "highlights.docx")))
if (status != 0) warning("highlights.docx not written")

# check: each copy compiles on its own (twice, for the cross-references)
for (doc in c("manuscript", "supplementary")) {
  for (k in 1:2) status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "-output-directory", out, file.path(out, paste0(doc, ".tex"))), stdout = FALSE)
  if (status != 0) stop(doc, ".tex does not compile with pdflatex; see ", file.path(out, paste0(doc, ".log")))
  unlink(file.path(out, paste0(doc, c(".aux", ".log", ".out", ".toc"))))
}
# the PDF to upload (bioRxiv, journal submission system) is the Typst one, not the LaTeX check build
for (doc in c("manuscript", "supplementary")) file.copy(file.path(src, "_site", paste0(doc, ".pdf")), file.path(out, paste0(doc, ".pdf")), overwrite = TRUE)
unlink(file.path(src, c("manuscript-latex.pdf", "supplementary-latex.pdf")))
message("submission/ ready: ", paste(list.files(out, recursive = TRUE), collapse = ", "))
