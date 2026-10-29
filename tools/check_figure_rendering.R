# ==============================================================================
# check_figure_rendering.R -- do the symbols in figure labels render as symbols in THIS R installation?
#   Rscript tools/check_figure_rendering.R
# Builds small plots with markup (tex() labels) in each text element of theme_spark -- axis titles, tick labels (discrete and
# continuous axes), facet strips, legend, title -- also inside a patchwork, and reports for each whether the element was drawn as
# formatted text (gridtext) or as plain text, in which case the markup (<i>x</i><sub>p</sub>) is shown literally.
# Paste the whole output if something fails: it contains the versions of the packages involved.
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggtext); library(patchwork); library(data.table); library(grid) })
source("config.R"); for (f in c("R/report/plots.R", "R/report/symbols.R")) source(f)

cat("R", as.character(getRversion()), "| ", paste(sprintf("%s %s", c("ggplot2", "ggtext", "gridtext", "patchwork", "systemfonts", "ragg"),
    vapply(c("ggplot2", "ggtext", "gridtext", "patchwork", "systemfonts", "ragg"), function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "-", "")), collapse = ", "), "\n", sep = "")
cat("Figure font: '", figure_font(cfg), "'\n\n", sep = "")

# is there formatted text (a gridtext richtext grob) inside the grob(s) whose layout name matches `pattern`?
has_richtext <- function(g) {
  if (inherits(g, "richtext_grob")) return(TRUE)
  kids <- if (inherits(g, "gtable")) g$grobs else if (inherits(g, "gTree")) g$children else if (inherits(g, "gList")) unclass(g) else NULL
  any(vapply(kids, has_richtext, logical(1)))
}
rich_in <- function(g, pattern) {
  if (!inherits(g, "gtable")) return(NA)
  hit <- logical(0)
  for (i in seq_along(g$grobs)) {
    name <- g$layout$name[i]
    if (grepl(pattern, name)) hit <- c(hit, has_richtext(g$grobs[[i]]))
    if (inherits(g$grobs[[i]], "gtable")) hit <- c(hit, rich_in(g$grobs[[i]], pattern))
  }
  hit <- hit[!is.na(hit)]
  if (!length(hit)) NA else any(hit)
}
grob_of <- function(p) if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplotGrob(p)

# measure text on the same kind of device that draws the figures (ragg when installed), so that no warnings about missing glyphs metrics appear
device_opened <- requireNamespace("ragg", quietly = TRUE)
if (device_opened) ragg::agg_png(tempfile(fileext = ".png"))
lab <- tex_named(c(x_p = "$x_p$", x_s = "$x_s$"))
df <- data.table(g = factor(rep(c("x_p", "x_s"), each = 4), c("x_p", "x_s")), v = c(1:4, 2:5), t = rep(1:4, 2), k = rep(c("a", "b"), 4))
base <- function() ggplot(df, aes(g, v)) + geom_point()
cases <- list(
  list("axis title (x)", base() + labs(x = tex("time $t$"), y = NULL) + theme_spark(cfg), "xlab-b"),
  list("axis title (y)", base() + labs(y = tex("$\\tau = 1/\\kappa$ (s)"), x = NULL) + theme_spark(cfg), "ylab-l"),
  list("tick labels, discrete x axis", base() + scale_x_discrete(labels = lab) + labs(x = NULL, y = NULL) + theme_spark(cfg), "axis-b"),
  list("tick labels, discrete y axis", ggplot(df, aes(v, g)) + geom_point() + scale_y_discrete(labels = lab) + labs(x = NULL, y = NULL) + theme_spark(cfg), "axis-l"),
  list("tick labels, continuous axis", ggplot(df, aes(t, v)) + geom_point() + scale_x_continuous(breaks = 1:4, labels = tex(c("$1$", "$\\alpha$", "$x_3$", "$4$"))) + labs(x = NULL, y = NULL) + theme_spark(cfg), "axis-b"),
  list("facet strips", base() + facet_wrap(~g, labeller = as_labeller(lab)) + labs(x = NULL, y = NULL) + theme_spark(cfg), "strip"),
  list("legend labels", ggplot(df, aes(t, v, colour = g)) + geom_point() + scale_colour_discrete(labels = lab) + labs(x = NULL, y = NULL) + theme_spark(cfg), "guide-box"),
  list("subtitle", base() + labs(x = NULL, y = NULL, subtitle = tex("below $D$ critical")) + theme_spark(cfg), "subtitle"))
check <- function(plot, pattern) rich_in(grob_of(plot), pattern)
res <- data.table(element = vapply(cases, `[[`, "", 1), alone = vapply(cases, function(x) check(x[[2]], x[[3]]), NA))
in_pw_and <- function(x) { pw <- (x[[2]] | x[[2]]) + plot_layout(guides = "collect") & theme(legend.position = "bottom"); check(pw, x[[3]]) }
in_pw_annot <- function(x) { pw <- (x[[2]] | x[[2]]) + plot_layout(guides = "collect") + plot_annotation(theme = theme(legend.position = "bottom")); check(pw, x[[3]]) }
res[, `patchwork with &` := vapply(cases, in_pw_and, NA)]
res[, `patchwork with plot_annotation()` := vapply(cases, in_pw_annot, NA)]
fmt <- function(x) ifelse(is.na(x), "n/a", ifelse(x, "formatted", "LITERAL MARKUP"))
print(res[, lapply(.SD, function(c) if (is.logical(c)) fmt(c) else c)], row.names = FALSE)
bad <- sum(unlist(res[, -1]) %in% FALSE)
cat(if (bad == 0) "\nAll checked elements render their symbols.\n" else sprintf("\n%d check(s) show the markup literally: send this output.\n", bad))
if (bad > 0) {
  cat("\nWorkarounds tested for the discrete x tick labels (TRUE = formatted):\n")
  p0 <- base() + scale_x_discrete(labels = lab) + labs(x = NULL, y = NULL)
  v1 <- p0 + theme_spark(cfg) + theme(axis.text.x = element_markdown(size = 9))
  v2 <- p0 + theme_spark(cfg) + theme(axis.text.x.bottom = element_markdown(size = 9))
  cat("  theme(axis.text.x = element_markdown()) after theme_spark:        ", check(v1, "axis-b"), "\n")
  cat("  theme(axis.text.x.bottom = element_markdown()) after theme_spark: ", check(v2, "axis-b"), "\n")
}
if (device_opened) invisible(dev.off())
