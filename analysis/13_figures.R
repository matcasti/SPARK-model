# ==============================================================================
# 13_figures.R -- regenerate all manuscript figures from the saved results.
#
# Reads   results/*.rds (whichever exist; a figure whose inputs are missing is skipped)
# Writes  manuscript/figures/<name>.<format>   (names in FIGURE_REGISTRY, R/report/figures.R)
#         results/figure_info.rds   which recording was showcased in the example figures, and why
# Look    cfg$figures in config.R: format, dpi, fonts, palettes, landscape settings and the
#         showcased recordings (cfg$figures$example$prcp / $exercise: "largest", "median" or an id).
#         To restyle or change the example, edit config.R and re-run only this script.
# ==============================================================================
source("analysis/00_setup.R")

# Symbols in labels (R/report/symbols.R) are HTML rendered by ggtext; ragg (when installed) draws them with full Unicode support
if (!requireNamespace("ragg", quietly = TRUE)) log_msg("Package 'ragg' not installed: using the default PNG device (install.packages('ragg') for the best rendering of symbols).")

figure_dir <- cfg$paths$figures
log_msg("Figure font: %s (symbols checked with figure_glyphs(); run tools/checks/check_figure_fonts.R if a symbol is missing)", figure_font(cfg))
log_msg("If a figure shows markup such as <i>x</i><sub>p</sub> instead of symbols, run tools/checks/check_figure_rendering.R and send its output.")
results <- function(name) read_result(name, cfg, required = FALSE)
format <- cfg$figures$format

figure_info <- list()
for (name in names(FIGURE_REGISTRY)) {
  tryCatch({
    figure <- FIGURE_REGISTRY[[name]](results, cfg)
    ggsave(file.path(figure_dir, paste0(name, ".", format)), figure$plot, width = figure$width, height = figure$height,
           dpi = cfg$figures$dpi, bg = "white", device = figure_device(format))
    if (!is.null(cfg$figures$vector_format) && cfg$figures$vector_format == "pdf") {
      ggsave(file.path(figure_dir, paste0(name, ".pdf")), figure$plot, width = figure$width, height = figure$height, device = grDevices::cairo_pdf, bg = "white")
    }
    if (length(figure$info)) figure_info[[name]] <- figure$info
    shown_in <- min(cfg$figures$journal_width_in, cfg$figures$text_width_in)           # axis text is base_size - 1 (theme_spark)
    final_pt <- (cfg$figures$base_size - 1) * min(1, shown_in / figure$width)
    log_msg("  wrote %s.%s%s", name, format, if (final_pt < cfg$figures$min_text_pt) sprintf("   [width %.1f in: axis text %.1f pt at %.2f in, below %g pt]", figure$width, final_pt, shown_in, cfg$figures$min_text_pt) else "")
  }, error = function(e) log_msg("  SKIPPED %s: %s", name, conditionMessage(e)))
}
write_result(figure_info, "figure_info", cfg)
log_msg("Figures written to %s", figure_dir)
