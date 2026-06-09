# Summarize local-only particle-filter profiling outputs.
# Does not read private data.

output_root <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
summary_path <- file.path(output_root, "particle_filter_hotspot_summary.csv")
inventory_path <- file.path(output_root, "sorting_rle_inventory.csv")
recommendation_path <- file.path(output_root, "next_optimization_recommendation.txt")

if (!file.exists(summary_path)) {
  stop("Hotspot summary not found. Run profile_particle_filter_hotspots.R first: ", summary_path, call. = FALSE)
}
if (!file.exists(inventory_path)) {
  stop("Sorting/RLE inventory not found. Run inventory_sorting_rle_operations.R first: ", inventory_path, call. = FALSE)
}

hotspots <- read.csv(summary_path)
inventory <- read.csv(inventory_path)

top_hotspots <- utils::head(hotspots[order(hotspots$self.time, decreasing = TRUE), ], 15)
inventory_summary <- aggregate(line_number ~ classification, inventory, length)
names(inventory_summary) <- c("classification", "n_occurrences")

recommendation <- c(
  "Next low-risk optimization recommendation",
  paste0("Date/time: ", Sys.time()),
  "",
  "Particle-history constraint:",
  "Do not collapse live particles to counts by current state. Particle identity, current state, trajectory history, weights, transition history, and ancestor/resampling relationships must remain aligned.",
  "",
  "Recommended next target:",
  "Focus on profiling-confirmed post-processing or stored-output RLE work first, especially get.transition.rle(), Points.rle/Transitions.rle creation, inverse.rle(), and get.coordinates.PF()/dist.fun if they appear near the top of the hotspot profile.",
  "",
  "Low-risk rule:",
  "A tabulate/counting rewrite is only safe for marginal summaries or stored RLE outputs after confirming downstream code does not require live particle identity. Do not apply count compression inside generate.points.dirs() or live propagation.",
  "",
  "Current evidence:",
  paste(capture.output(print(top_hotspots)), collapse = "\n"),
  "",
  "Sorting/RLE inventory summary:",
  paste(capture.output(print(inventory_summary)), collapse = "\n")
)

writeLines(recommendation, recommendation_path)
print(top_hotspots)
print(inventory_summary)
message("Recommendation saved to: ", recommendation_path)
