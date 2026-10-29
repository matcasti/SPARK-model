# ==============================================================================
# R/report/symbols.R -- mathematical notation in figures and tables, from one source.
#
#   tex()          LaTeX-style text -> HTML for ggtext (figure titles, axes, strips, legends):
#                  "Sympathetic-like state $x_s$ and $\tau_p$ (s)" -> italic x with subscript s, Greek tau, ...
#   SYMBOL_TEX     name of a variable in the result tables (Tau_P, sd_x_p, ...) -> its LaTeX
#   sym()          the same as HTML for figures;   sym_tex()  as LaTeX for tables ($...$ in knitr::kable)
#   table_names()  column / row labels for knitr::kable written in LaTeX
#
# Conventions (as in the manuscript): variables and lower-case Greek letters in italics, subscripts and
# labels upright, upper-case Greek letters and functions (ln, SD) upright, a true minus sign (U+2212).
# Supported in math: \alpha ... \omega, \Gamma \Delta \Theta \Lambda \Sigma \Phi \Psi \Omega, _ and ^ with one
# character or {group}, \hat{}, \mathrm{}, \text{}, \ln \log \exp, \times \cdot \pm \approx \leq \geq \neq \to
# \rightarrow \infty \partial \, \; \%. Anything else is a bug in the label and is reported by tex().
# ==============================================================================


GREEK <- c(alpha = "\u03b1", beta = "\u03b2", gamma = "\u03b3", delta = "\u03b4", epsilon = "\u03b5", varepsilon = "\u03b5",
           zeta = "\u03b6", eta = "\u03b7", theta = "\u03b8", kappa = "\u03ba", lambda = "\u03bb", mu = "\u03bc", nu = "\u03bd",
           xi = "\u03be", pi = "\u03c0", rho = "\u03c1", sigma = "\u03c3", tau = "\u03c4", phi = "\u03c6", varphi = "\u03c6",
           chi = "\u03c7", psi = "\u03c8", omega = "\u03c9")
GREEK_UPPER <- c(Gamma = "\u0393", Delta = "\u0394", Theta = "\u0398", Lambda = "\u039b", Xi = "\u039e", Pi = "\u03a0",
                 Sigma = "\u03a3", Phi = "\u03a6", Psi = "\u03a8", Omega = "\u03a9")
MATH_SYMBOLS <- c(`in` = "\u2208", times = "\u00d7", cdot = "\u00b7", pm = "\u00b1", approx = "\u2248", leq = "\u2264", le = "\u2264", geq = "\u2265",
                  ge = "\u2265", neq = "\u2260", ne = "\u2260", to = "\u2192", rightarrow = "\u2192", leftarrow = "\u2190",
                  infty = "\u221e", partial = "\u2202", sim = "\u223c", propto = "\u221d", ldots = "\u2026", degree = "\u00b0")
RELATIONS <- c("in", "approx", "leq", "le", "geq", "ge", "neq", "ne", "to", "rightarrow", "leftarrow", "sim", "propto")
MATH_FUNCTIONS <- c("ln", "log", "exp", "max", "min", "SD", "Var", "Cov", "cor", "sin", "cos", "sup")

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("*", "&#42;", x, fixed = TRUE)          # not markdown emphasis
  gsub("_", "&#95;", x, fixed = TRUE)
}

#' Math (the inside of $...$) to HTML. Small recursive-descent parser over characters.
tex_math <- function(math) {
  chars <- strsplit(math, "")[[1]]
  n <- length(chars)
  position <- 1L
  peek <- function() if (position <= n) chars[position] else ""
  advance <- function() { value <- chars[position]; position <<- position + 1L; value }
  skip_spaces <- function() while (peek() == " ") position <<- position + 1L
  script_depth <- 0L                                   # > 0 inside a subscript or superscript: no operator spacing there
  previous_char <- function(from) {                     # last non-space source character before `from`
    i <- from - 1L
    while (i >= 1L && chars[i] == " ") i <- i - 1L
    if (i >= 1L) chars[i] else ""
  }

  read_command <- function() {                       # after the backslash
    if (grepl("[A-Za-z]", peek())) {
      name <- ""
      while (grepl("[A-Za-z]", peek())) name <- paste0(name, advance())
      name
    } else advance()
  }
  # A group {...} or one atom, as HTML; `upright` forces roman letters (subscripts, \mathrm)
  read_argument <- function(upright) {
    skip_spaces()
    if (peek() == "{") {
      advance()
      html <- read_until_close(upright)
      return(html)
    }
    read_atom(upright)
  }
  read_until_close <- function(upright) {
    html <- ""
    while (position <= n && peek() != "}") html <- paste0(html, read_item(upright))
    if (peek() == "}") advance()
    html
  }
  read_atom <- function(upright) {
    ch <- advance()
    if (ch == "\\") {
      command <- read_command()
      if (command %in% names(GREEK)) return(if (upright) GREEK[[command]] else paste0("<i>", GREEK[[command]], "</i>"))
      if (command %in% names(GREEK_UPPER)) return(GREEK_UPPER[[command]])
      if (command %in% names(MATH_SYMBOLS)) {
        symbol <- MATH_SYMBOLS[[command]]
        return(if (command %in% RELATIONS) paste0("\u2005", symbol, "\u2005") else symbol)       # relations get a space on both sides
      }
      if (command %in% MATH_FUNCTIONS) return(if (peek() %in% c("(", "_", "^", "")) command else paste0(command, "\u2009"))
      if (command == ",") return("\u2009")
      if (command %in% c(";", ":", " ")) return("\u2005")
      if (command == "quad") return("\u2003")
      if (command %in% c("%", "{", "}", "$", "&", "#")) return(html_escape(command))
      if (command == "sqrt") return(paste0("\u221a", read_argument(upright)))      # no overline in HTML: written as \u221a\u03bd, \u221a(N)
      if (command %in% c("hat", "tilde", "bar")) {      # combining accent goes inside the italic tag
        accent <- c(hat = "\u0302", tilde = "\u0303", bar = "\u0304")[[command]]
        base <- read_argument(upright)
        return(if (grepl("</i>$", base)) sub("</i>$", paste0(accent, "</i>"), base) else paste0(base, accent))
      }
      if (command %in% c("left", "right")) return("")
      if (command %in% c("mathrm", "text", "textrm", "operatorname")) return(read_argument(TRUE))
      if (command %in% c("mathit", "textit", "mathcal")) return(paste0("<i>", read_argument(TRUE), "</i>"))
      if (command %in% c("mathbf", "textbf")) return(paste0("<b>", read_argument(TRUE), "</b>"))
      stop(sprintf("tex(): unsupported command \\%s in '%s'", command, math), call. = FALSE)
    }
    if (ch == "{") return(read_until_close(upright))
    if (ch == "-") {                                                       # binary minus (after an operand) is spaced, unary is not
      binary <- script_depth == 0L && previous_char(position - 1L) %in% c(letters, LETTERS, 0:9, ")", "}", "]")
      return(if (binary) "\u2005\u2212\u2005" else "\u2212")
    }
    if (ch == ",") return(if (script_depth == 0L) ",\u2005" else ",")      # space after a comma inside math
    if (ch %in% c("=", "+", "<", ">")) {                                   # binary operators get a space on both sides
      return(if (script_depth == 0L) paste0("\u2005", html_escape(ch), "\u2005") else html_escape(ch))
    }
    if (grepl("[A-Za-z]", ch)) return(if (upright) ch else paste0("<i>", ch, "</i>"))
    html_escape(ch)
  }
  read_item <- function(upright) {
    skip_spaces()
    if (position > n) return("")
    atom <- if (peek() %in% c("_", "^")) "" else read_atom(upright)
    repeat {
      skip_spaces()
      if (peek() == "_") {
        advance()
        first <- if (peek() == "{" && position < n) chars[position + 1L] else peek()
        index <- first %in% c("i", "j", "k", "n", "m", "t")       # running indices stay italic, labels (p, s, sys) upright
        script_depth <<- script_depth + 1L
        atom <- paste0(atom, "<sub>", read_argument(!index), "</sub>")
        script_depth <<- script_depth - 1L
      }
      else if (peek() == "^") { advance(); script_depth <<- script_depth + 1L
        atom <- paste0(atom, "<sup>", read_argument(FALSE), "</sup>")
        script_depth <<- script_depth - 1L }
      else break
    }
    atom
  }

  # A run of letters that is a known function name (ln, SD, ...) stays upright
  protected <- math
  for (fn in MATH_FUNCTIONS) protected <- gsub(sprintf("(?<![A-Za-z\\\\])%s(?![A-Za-z])", fn), sprintf("\\\\mathrm{%s}", fn), protected, perl = TRUE)
  chars <- strsplit(protected, "")[[1]]
  n <- length(chars)
  position <- 1L
  html <- ""
  while (position <= n) html <- paste0(html, read_item(FALSE))
  gsub("</i><i>", "", html, fixed = TRUE)            # "LL" is one italic run
}

#' LaTeX-style text -> HTML for ggtext. Math is written between dollar signs; "\n" becomes a line break;
#' **text** is bold. (ggtext turns the characters "(c)" into a copyright sign: write "panel c".)
#' Vectorized; NA stays NA.
#' @examples
#' tex("Sympathetic-like state $x_s$ vs $\\tau_p$ (s)")
tex <- function(x) {
  one <- function(s) {
    if (is.na(s)) return(NA_character_)
    pieces <- regmatches(s, gregexpr("(?<!\\\\)\\$[^$]*\\$", s, perl = TRUE))[[1]]
    outside <- regmatches(s, gregexpr("(?<!\\\\)\\$[^$]*\\$", s, perl = TRUE), invert = TRUE)[[1]]
    out <- html_escape(outside[1])
    for (i in seq_along(pieces)) {
      out <- paste0(out, tex_math(substring(pieces[i], 2, nchar(pieces[i]) - 1)), html_escape(outside[i + 1]))
    }
    out <- gsub("&#42;&#42;(.+?)&#42;&#42;", "<b>\\1</b>", out, perl = TRUE)     # **bold**
    gsub("\n", "<br>", gsub("\\$", "$", out, fixed = TRUE), fixed = TRUE)
  }
  vapply(x, one, character(1), USE.NAMES = FALSE)
}

#' A named vector of tex() labels (names kept), e.g. for scale_*_manual(labels = ).
tex_named <- function(x) setNames(tex(x), names(x))


# ---- Variables of the result tables ------------------------------------------------------------------------------
#' LaTeX of every variable that reaches a figure or a table (column names of the summary tables).
SYMBOL_TEX <- c(
  # latent summaries and parameters
  Nu_0 = "$\\nu_0$", nu0 = "$\\nu_0$", Tau_P = "$\\tau_p$", Tau_S = "$\\tau_s$", tau_p = "$\\tau_p$", tau_s = "$\\tau_s$",
  Kappa_P = "$\\kappa_p$", Kappa_S = "$\\kappa_s$", kp = "$\\kappa_p$", ks = "$\\kappa_s$",
  lambda_P = "$\\lambda_P$", Lambda_P = "$\\lambda_P$", lp = "$\\lambda_P$", Lambda_R = "$\\lambda_R$", lR = "$\\lambda_R$",
  Sigma_sys2 = "$\\hat\\sigma^2_{sys}$", sigma2 = "$\\sigma^2_{sys}$",
  rr_jitter_ms = "RR timing jitter (ms)", sd_x_p = "SD of $x_p$", sd_x_s = "SD of $x_s$", sd_ratio = "SD ratio $x_p/x_s$",
  log_kappa_ratio = "$\\ln(\\kappa_p/\\kappa_s)$", share_p = "Share of $x_p$ variance",
  x_p = "$x_p$", x_s = "$x_s$", X_p = "$x_p$", X_s = "$x_s$",
  # conventional HRV and physiology
  mean_hr = "Mean HR", SDNN = "SDNN", RMSSD = "RMSSD", LF_HF = "LF/HF", LF = "LF power", HF = "HF power",
  f_resp = "Respiratory frequency", resp_sd_x_p = "Respiratory-band SD of $x_p$", resp_sd_x_s = "Respiratory-band SD of $x_s$",
  resp_share_x_p = "Respiratory-band share of $x_p$", resp_share_x_s = "Respiratory-band share of $x_s$",
  alpha_LF = "Baroreflex $\\alpha$ (LF)", alpha_HF = "Baroreflex $\\alpha$ (HF)", sbp = "Systolic pressure", map = "MAP",
  pulse_pressure = "Pulse pressure", hr = "HR",
  # perturbation responses
  d_x_p = "$\\Delta x_p$", d_x_s = "$\\Delta x_s$", z_d_x_p = "$\\Delta x_p$ (baseline SD)", z_d_x_s = "$\\Delta x_s$ (baseline SD)",
  tau_p_resp = "$\\tau_p$ (response)", tau_s_resp = "$\\tau_s$ (response)"
)

#' LaTeX of variable names (unknown names are returned unchanged, with underscores made readable).
sym_tex <- function(names) {
  known <- names %in% names(SYMBOL_TEX)
  ifelse(known, unname(SYMBOL_TEX[match(names, names(SYMBOL_TEX))]), gsub("_", " ", names, fixed = TRUE))
}

#' HTML of variable names for ggtext (figures).
sym <- function(names) tex(sym_tex(names))

#' Row / column labels written in LaTeX for knitr::kable (Quarto renders $...$ in HTML and PDF).
#' @param labels character vector in which `$...$` is math; returned unchanged (kept as a function so that every
#'   table label goes through one place, and so that '|' and '%' inside math do not break pipe tables)
table_names <- function(labels) gsub("|", "\\vert ", labels, fixed = TRUE)


# ---- Names of specifications, scenarios and sensitivity settings -------------------------------------------------
SPEC_TEX <- c(M0 = "$M_0$", M1 = "$M_1$", M2 = "$M_2$", M3 = "$M_3$")
spec_tex <- function(x) ifelse(x %in% names(SPEC_TEX), unname(SPEC_TEX[match(x, names(SPEC_TEX))]), x)

# ---- Names of analysis choices in words (a code identifier never appears in the text or on an axis) --------------------
RECOVERY_DEFINITION_TEXT <- c(sustained = "sustained return to within 1 rest SD", first = "first return to within 1 rest SD",
                              frac90 = "90% recovery of the excursion", expfit = "exponential time constant")
END_RULE_TEXT <- c(peak = "heart-rate peak", cp = "second change point", profile = "offset time of the fitted exercise profile")
RESPONSE_TEXT <- c(z_d_x_p = "standardized change of the vagal-like state", z_d_x_s = "standardized change of the sympathetic-like state",
                   d_mean_hr = "change in mean heart rate", d_SDNN = "change in SDNN", d_RMSSD = "change in RMSSD")
words_for <- function(x, labels) ifelse(x %in% names(labels), unname(labels[match(x, names(labels))]), x)

SCENARIO_LABELS <- c(in_prior = "True time scales\nwithin the priors", off_prior = "True time scales\noutside the priors",
                     vagal_dom = "Vagal-dominant\nvariance", symp_dom = "Sympathetic-dominant\nvariance")

#' Readable, tex()-ready label of a sensitivity setting ("tau_p_bounds=narrow", "M3 without tau").
setting_axis_label <- function(setting) {
  base <- c(tau_p_bounds = "$\\tau_p$ interval", tau_s_bounds = "$\\tau_s$ interval", topo_target_ratio = "Separation target $\\kappa_p/\\kappa_s$",
            topo_prior_weight = "Separation weight", jitter_bounds_ms = "Timing-jitter interval", crossover_hz = "Spectral crossover (Hz)",
            sns_leakage_at_boundary = "Sympathetic leakage", lower_lf_vagal_fraction = "Vagal floor fraction")
  term <- c(tau = "time-scale priors", topo = "separation prior", jitter = "timing-jitter prior", anchor = "spectral anchor")
  vapply(setting, function(x) {
    if (grepl("^M3 without ", x)) return(paste0("$M_3$ without the ", term[[sub("^M3 without ", "", x)]]))
    parts <- strsplit(x, "=", fixed = TRUE)[[1]]
    if (parts[1] %in% names(base)) paste0(base[[parts[1]]], ": ", parts[2]) else x
  }, character(1), USE.NAMES = FALSE)
}


# ---- Fonts and device -----------------------------------------------------------------------------------------------
#' Serif family for equations (the first installed of STIX, Times, TeX Gyre Termes, Liberation Serif, DejaVu Serif).
math_font <- function(candidates = c("STIX Two Text", "STIXGeneral", "Times New Roman", "Times", "TeX Gyre Termes", "Liberation Serif", "DejaVu Serif")) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return("serif")
  available <- unique(systemfonts::system_fonts()$family)
  found <- candidates[candidates %in% available]
  if (length(found)) found[1] else "serif"
}
#' Characters that figure labels can contain: lower- and upper-case Greek letters, the true minus sign, the operators and
#' symbols of tex() and the superscript two. A font used for figures must have all of them.
figure_glyphs <- function() {
  unique(c(unname(GREEK), unname(GREEK_UPPER), "\u2212", "\u00b2",
           unname(MATH_SYMBOLS[c("times", "cdot", "pm", "approx", "leq", "geq", "neq", "to", "infty", "partial", "degree")])))
}

#' Characters of `chars` that the font family lacks (character(0): complete; NA: the font could not be inspected).
font_missing_glyphs <- function(family, chars = figure_glyphs()) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(NA_character_)
  info <- tryCatch(systemfonts::glyph_info(chars, family = family), error = function(e) NULL)
  if (is.null(info)) return(NA_character_)
  chars[info$index == 0]
}

#' The font family for figures. The first installed family of `candidates` (cfg$figures$font first) that has every symbol
#' used in the labels; a family that lacks some (for example the Helvetica of some systems has no Greek letters) is skipped,
#' because the missing symbols would be drawn as empty boxes or not at all, mainly in the PDF versions. When none is complete the
#' installed family missing the fewest symbols is used and the missing ones are reported. "" (the graphics default) when no font is found.
figure_font <- function(cfg, candidates = c("Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "DejaVu Sans", "FreeSans", "Noto Sans", "Arial Unicode MS")) {
  wanted <- unique(c(cfg$figures$font, candidates))
  wanted <- wanted[nzchar(wanted)]
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(if (length(wanted)) wanted[1] else "")
  key <- paste(wanted, collapse = "|")
  if (!is.null(FIGURE_FONT_CACHE[[key]])) return(FIGURE_FONT_CACHE[[key]])
  available <- unique(systemfonts::system_fonts()$family)
  aliases <- c(Helvetica = "Arial", Arial = "Liberation Sans")      # metric-compatible substitutes
  installed <- character(0)
  for (family in wanted) {
    use <- if (family %in% available) family else if (family %in% names(aliases) && aliases[[family]] %in% available) aliases[[family]] else NA_character_
    if (!is.na(use) && !use %in% installed) installed <- c(installed, use)
  }
  chosen <- ""
  if (length(installed)) {
    missing <- lapply(installed, font_missing_glyphs)
    complete <- vapply(missing, function(m) length(m) == 0L, logical(1))
    if (any(complete)) {
      chosen <- installed[complete][1]
      if (!identical(chosen, installed[1]))
        message(sprintf("Figure font: '%s' lacks %s; using '%s' instead.", installed[1], paste(missing[[1]], collapse = " "), chosen))
    } else {
      n_missing <- vapply(missing, function(m) if (length(m) == 1L && is.na(m)) Inf else length(m), numeric(1))
      chosen <- installed[which.min(n_missing)]
      message(sprintf("Figure font: no installed family has every symbol; using '%s', which lacks %s. Install DejaVu Sans or another font with Greek and mathematical symbols.", chosen, paste(missing[[which.min(n_missing)]], collapse = " ")))
    }
  }
  FIGURE_FONT_CACHE[[key]] <- chosen
  chosen
}
FIGURE_FONT_CACHE <- new.env()

#' Graphics device for ggsave(): ragg (anti-aliased, full Unicode, font fallback) when installed.
figure_device <- function(format = "png") {
  if (format == "png" && requireNamespace("ragg", quietly = TRUE)) return(ragg::agg_png)
  if (format == "svg" && requireNamespace("svglite", quietly = TRUE)) return(svglite::svglite)
  NULL
}
