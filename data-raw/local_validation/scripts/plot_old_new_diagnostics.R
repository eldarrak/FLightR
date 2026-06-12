output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

bench <- utils::read.csv(file.path(output_dir, "benchmark_combined.csv"), stringsAsFactors = FALSE)
gps_sum <- utils::read.csv(file.path(output_dir, "gps_error_summary_by_result.csv"), stringsAsFactors = FALSE)
comparison <- utils::read.csv(file.path(output_dir, "benchmark_old_new_comparison.csv"), stringsAsFactors = FALSE)
by_twilight <- utils::read.csv(file.path(output_dir, "gps_old_new_by_twilight.csv"), stringsAsFactors = FALSE)
by_period <- utils::read.csv(file.path(output_dir, "gps_old_new_by_period.csv"), stringsAsFactors = FALSE)

version_cols <- c(old_master = "steelblue", new_branch = "firebrick", unknown = "grey50")
col_for <- function(version) version_cols[ifelse(version %in% names(version_cols), version, "unknown")]

png_plot <- function(name, expr) {
  grDevices::png(file.path(plot_dir, name), width = 1100, height = 800, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

plot_runtime <- function(y, ylab, file) {
  png_plot(file, {
    plot(bench$nParticles, bench[[y]], log = "x", pch = 19,
         col = col_for(bench$version), xlab = "nParticles", ylab = ylab)
    legend("topleft", legend = names(version_cols), col = version_cols, pch = 19, bty = "n")
  })
}

plot_gps <- function(y, ylab, file) {
  png_plot(file, {
    plot(gps_sum$nParticles, gps_sum[[y]], log = "x", pch = 19,
         col = col_for(gps_sum$version), xlab = "nParticles", ylab = ylab)
    legend("topleft", legend = names(version_cols), col = version_cols, pch = 19, bty = "n")
  })
}

if (nrow(bench) > 0) {
  plot_runtime("particle_filter_seconds", "Particle filter runtime (seconds)", "runtime_particle_filter_by_nParticles.png")
  plot_runtime("total_seconds", "Total runtime (seconds)", "runtime_total_by_nParticles.png")
}

if (nrow(gps_sum) > 0) {
  plot_gps("median_error_km", "Median GPS error (km)", "gps_median_error_by_nParticles.png")
  plot_gps("q75_error_km", "Q75 GPS error (km)", "gps_q75_error_by_nParticles.png")
  plot_gps("q95_error_km", "Q95 GPS error (km)", "gps_q95_error_by_nParticles.png")
}

if (nrow(by_twilight) > 0) {
  by_twilight$twilight_time <- as.POSIXct(by_twilight$twilight_time, tz = "GMT")
  png_plot("old_new_delta_error_through_time.png", {
    plot(by_twilight$twilight_time, by_twilight$delta_error_km, type = "p", pch = 19,
         xlab = "Twilight time", ylab = "Delta error km (new - old)")
    abline(h = 0, lty = 2)
  })
  png_plot("old_vs_new_error_scatter.png", {
    plot(by_twilight$error_km_old, by_twilight$error_km_new, pch = 19,
         xlab = "Old GPS error (km)", ylab = "New GPS error (km)")
    abline(0, 1, lty = 2)
  })
}

if (nrow(by_period) > 0 && nrow(by_twilight) > 0) {
  png_plot("delta_error_by_period_boxplot.png", {
    boxplot(delta_error_km ~ period, data = by_twilight,
            ylab = "Delta error km (new - old)", xlab = "Period")
    abline(h = 0, lty = 2)
  })
}

message("Plots written under: ", plot_dir)

