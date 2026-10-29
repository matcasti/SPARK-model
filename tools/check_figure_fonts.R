# ==============================================================================
# check_figure_fonts.R -- which font the figures use and whether it has every symbol of the labels.
#   Rscript tools/check_figure_fonts.R
# Prints, for each candidate family, the symbols (Greek letters, minus, operators) it lacks, the family chosen by
# figure_font() and the packages that draw the text. If symbols are missing or empty boxes appear in a figure, install a
# font with full Greek and mathematical coverage (DejaVu Sans, Arial, Noto Sans) or set cfg$figures$font in config.R.
# ==============================================================================
source("config.R"); source("R/report/symbols.R")
cat("Packages: ")
for (p in c("systemfonts", "ragg", "ggtext", "gridtext", "svglite")) cat(sprintf("%s %s  ", p, if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "NOT INSTALLED"))
cat("\n\nSymbols checked (", length(figure_glyphs()), "): ", paste(figure_glyphs(), collapse = " "), "\n\n", sep = "")
if (!requireNamespace("systemfonts", quietly = TRUE)) stop("install.packages('systemfonts')")
available <- unique(systemfonts::system_fonts()$family)
for (family in unique(c(cfg$figures$font, "Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "DejaVu Sans", "FreeSans", "Noto Sans", "Arial Unicode MS"))) {
  if (!family %in% available) { cat(sprintf("%-18s not installed\n", family)); next }
  m <- font_missing_glyphs(family)
  cat(sprintf("%-18s %s\n", family, if (length(m) == 1L && is.na(m)) "could not be inspected" else if (!length(m)) "complete" else paste("lacks:", paste(m, collapse = " "))))
}
cat("\nFamily used by the figures: '", figure_font(cfg), "'\n", sep = "")
