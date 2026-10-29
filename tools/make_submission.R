# ==============================================================================
# tools/make_submission.R -- gathers the files the journal asks for (the PDF is also what bioRxiv takes) in submission/ (not under version control)
#   make submission          renders the LaTeX route (keeping the .tex; the PDF itself comes from `make render`, i.e. Typst) and then runs this script
# Writes: manuscript.tex (single column, compiles with pdfLaTeX), Figure_N.pdf (vector, in order of appearance),
#         highlights.docx (separate editable file, as the journal asks), manuscript.pdf and supplementary.pdf (Typst); compiles the .tex as a check.
# ==============================================================================
out <- "submission"; src <- "manuscript"
tex <- file.path(src, "manuscript.tex")
if (!file.exists(tex)) stop("render the PDF first (`make submission`): ", tex, " not found")
unlink(out, recursive = TRUE); dir.create(file.path(out, "figures"), recursive = TRUE)

lines <- readLines(tex, encoding = "UTF-8")
png <- unique(unlist(regmatches(lines, gregexpr("figures/[A-Za-z0-9_]+\\.png", lines))))   # order of appearance
for (i in seq_along(png)) {
  pdf <- file.path(src, sub("\\.png$", ".pdf", png[i]))
  if (!file.exists(pdf)) stop("vector version missing: ", pdf)
  new <- sprintf("figures/Figure_%d.pdf", i)
  file.copy(pdf, file.path(out, new), overwrite = TRUE)
  lines <- gsub(png[i], new, lines, fixed = TRUE)
}
writeLines(lines, file.path(out, "manuscript.tex"), useBytes = TRUE)

status <- system2("quarto", c("pandoc", file.path(src, "highlights.md"), "-o", file.path(out, "highlights.docx")))
if (status != 0) warning("highlights.docx not written")
sup <- file.path(src, "_site", "supplementary.pdf")
if (file.exists(sup)) file.copy(sup, file.path(out, "supplementary.pdf"), overwrite = TRUE)

# check: the copy compiles on its own (twice, for the cross-references)
for (k in 1:2) status <- system2("pdflatex", c("-interaction=nonstopmode", "-halt-on-error", "-output-directory", out, file.path(out, "manuscript.tex")), stdout = FALSE)
if (status != 0) stop("manuscript.tex does not compile with pdflatex; see ", file.path(out, "manuscript.log"))
unlink(file.path(out, c("manuscript.aux", "manuscript.log", "manuscript.out", "manuscript.toc")))
# the PDF to upload (bioRxiv, journal submission system) is the Typst one, not the LaTeX check build
file.copy(file.path(src, "_site", "manuscript.pdf"), file.path(out, "manuscript.pdf"), overwrite = TRUE)
unlink(file.path(src, "_site", "manuscript-latex.pdf"))
message("submission/ ready: ", paste(list.files(out, recursive = TRUE), collapse = ", "))
