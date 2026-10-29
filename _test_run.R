
# Prepare workspace -------------------------------------------------------

## Load packages
library(data.table)

## Load function
source("engine/model_engine.R")

spark_control()
spark_priors()

# -------------------------------------------------------------------------

## Load data
rr_data <- (fread("data/tmst/133.csv"))
rr_data <- (fread("data/tmst/136.csv"))
rr_data <- (fread("data/tmst/007.csv"))
rr_data <- (fread("data/tmst/008.csv"))
rr_data <- (fread("data/prcp/12726.csv")[1:1000])
rr_data <- (fread("data/tmst/009.csv"))
rr_data <- (fread("data/prcp/12734.csv"))


(fit <- fit_model(rr_data$rri)) |>
  (\(x) {
    fig <- visualize_model(x)
    print(fig)
    #report_model(fit)
  })()


rr_data$X_s <- c(fit$states[,"X_s"])
rr_data$X_p <- c(fit$states[,"X_p"])

rr_data <- na.omit(rr_data)

# 1. Calculate the 2D kernel density grid from your states
dens <- MASS::kde2d(rr_data$X_p, rr_data$X_s, n = 100)

# 2. Extract the Z (elevation) for each state to lift the path onto the surface
ix <- findInterval(rr_data$X_p, dens$x)
iy <- findInterval(rr_data$X_s, dens$y)

# Prevent out-of-bounds indices at the extreme edges
ix <- pmax(1, pmin(ix, length(dens$x)))
iy <- pmax(1, pmin(iy, length(dens$y)))

# Extract corresponding density height for each state
vec_z <- dens$z[cbind(ix, iy)]

# Panel B: 3D View
fig_3d <- wrap_elements(full = ~ {

  old_par <- par(mar = c(0, 1, 0, 0))
  # Draw surface and store projection matrix
  pmat <- persp(dens$x, dens$y, -dens$z,
                theta = 25, phi = 40, box = T, expand = 0.5,
                col = "dodgerblue3", shade = 0.4, border = NA,
                xlab = "PNS", ylab = "SNS", zlab = "ANS Topology")

  # Project the 3D trajectory (x, y, z) onto the 2D plot window
  pts <- trans3d(rr_data$X_p, rr_data$X_s, -vec_z + min(vec_z), pmat)

  # Draw the path
  lines(pts, col = "#FFF2")
  points(pts, col = "white", pch = 16, cex = 0.15)

  # Highlight the Start (green) and End (gold) of the trajectory
  points(pts$x[1], pts$y[1], col = "green3", pch = 19, cex = 1.5)
  points(pts$x[length(rr_data$X_p)], pts$y[length(rr_data$X_p)], col = "gold", pch = 19, cex = 1.5)

  par(old_par)
})

# Combine panels
fig <- (visualize_model(fit) / fig_3d) +
  patchwork::plot_layout(heights = c(0.5, 0.5))

ggsave(filename = "manuscript/figures/3d_density_hq.png", plot = fig,
       height = 900, width = 600, scale = 7, units = "px", dpi = 400)
