# =============================================================================
# Shared figure styling
# Sourced in every figure script after loading ggplot2
# =============================================================================
suppressMessages(library(ggplot2))

STROKE_PT <- 0.5
STROKE    <- STROKE_PT / (72.27 / 25.4)   # ggplot linewidth that renders at 0.5 pt

# geom-level defaults (survive theme changes) -------------------------------
local({
  line_geoms <- c("line", "path", "step", "segment", "hline", "vline", "abline",
                  "col", "bar", "boxplot", "errorbar", "errorbarh", "linerange",
                  "pointrange", "area", "ribbon", "density", "smooth", "tile",
                  "rect", "crossbar", "contour", "function")
  for (g in line_geoms) try(update_geom_defaults(g, list(linewidth = STROKE)), silent = TRUE)
  try(update_geom_defaults("point",     list(stroke = STROKE)), silent = TRUE)
  try(update_geom_defaults("count",     list(stroke = STROKE)), silent = TRUE)
})

# theme addon for axis lines/ticks (add AFTER theme_pubr()/theme_*() in a plot)
theme_strokes <- theme(
  axis.line    = element_line(linewidth = STROKE),
  axis.ticks   = element_line(linewidth = STROKE),
  panel.border = element_blank()
)
