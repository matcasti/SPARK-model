# EXPERIMENTAL. One simulated and one real exercise recording: the states under M3 (filtered), M3O+ (oscillator + RTS) and M3E (events estimated from the RR).
# Run from the project root:  Rscript experiments/m3plus/compare/compare_figure.R   -> experiments/m3plus/_out/compare_figure.png
src <- readLines("experiments/m3plus/compare/compare_sim.R"); source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
eval(parse(text = src[grep("^osc_draw", src):grep("^one <- function", src)[1] - 1]))
DUR <- 420; suppressMessages(library(ggplot2)); suppressMessages(library(patchwork))
models_for <- function(dy) { series <- prepare_series(dy); dy <- series$dy; ev <- events_from_rr(dy, "exercise", cfg$segmentation)
  f3 <- fit_ext(series); fo <- fit_ext(series, list(osc = TRUE)); fe <- if (is.null(ev)) f3 else fit_ext(series, list(inputs = list(events = ev)))
  list(t = cumsum(dy), ev = ev, M3 = state_estimates(f3, FALSE), `M3O+` = state_estimates(fo, TRUE), M3E = state_estimates(fe, FALSE), taus = c(M3 = f3$tau_p, `M3O+` = fo$tau_p, M3E = fe$tau_p)) }
long <- function(m, truth = NULL) { d <- rbindlist(lapply(c("M3", "M3O+", "M3E"), function(k) data.table(t = m$t, model = k, vagal = m[[k]]$xp, symp = m[[k]]$xs)))
  if (!is.null(truth)) d <- rbind(d, data.table(t = m$t, model = "truth", vagal = truth$xp, symp = truth$xs)); d }
cols <- c(truth = "black", M3 = "#b8860b", `M3O+` = "#2a9d8f", M3E = "#c1121f")
panel <- function(d, y, title, ylab) ggplot(d, aes(t, .data[[y]], colour = model, linewidth = model, alpha = model)) + geom_line() + scale_colour_manual(values = cols) +
  scale_linewidth_manual(values = c(truth = 0.9, M3 = 0.5, `M3O+` = 0.5, M3E = 0.6)) + scale_alpha_manual(values = c(truth = 1, M3 = 0.9, `M3O+` = 0.9, M3E = 0.9)) + labs(title = title, x = "Time (s)", y = ylab, colour = NULL) + guides(linewidth = "none", alpha = "none") + theme_minimal(base_size = 10) + theme(legend.position = "bottom")
# simulated exercise recording (truth known)
job <- list(scenario = "exercise", sim = 3, seed = id_seed("compare exercise 3", cfg$seed)); z <- build("exercise", job$seed); sm <- models_for(z$s$rr)
dsim <- long(sm, list(xp = z$s$x_p + z$s$z1, xs = z$s$x_s))
# real exercise recording closest to the median standardized response
d <- load_dataset("exercise", cfg); id <- "001"; rm <- models_for(d$rr[[id]]); dreal <- long(rm)
shade <- function(m) if (is.null(m$ev)) NULL else annotate("rect", xmin = m$ev$t0[1], xmax = m$ev$t1[1], ymin = -Inf, ymax = Inf, alpha = 0.07)
p <- (panel(dsim, "vagal", "Simulated exercise: vagal total (OU + oscillator + response)", "Hz") + shade(sm)) / (panel(dsim, "symp", "Simulated exercise: sympathetic total", "Hz") + shade(sm)) |
     (panel(dreal, "vagal", sprintf("Real exercise recording %s: vagal-like state", id), "Hz") + shade(rm)) / (panel(dreal, "symp", "Real: sympathetic-like state", "Hz") + shade(rm))
p <- p + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("experiments/m3plus/_out/compare_figure.png", p, width = 12, height = 7.5, dpi = 150); cat("tau_p simulated (true", round(z$tau_p, 1), "):", paste(names(sm$taus), round(sm$taus, 1), collapse = ", "), "| real:", paste(names(rm$taus), round(rm$taus, 1), collapse = ", "), "\n")
