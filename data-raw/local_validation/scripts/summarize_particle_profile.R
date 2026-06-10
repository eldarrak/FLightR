# Summarize local particle-filter phase profiles and Rprof output.
# Reads aggregate diagnostic outputs only; does not read private GPS/GLS data.

output_root <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
canonical_label <- Sys.getenv("FLIGHTR_PROFILE_SUMMARY_LABEL", unset = "full_year_1e4_threads4_seed123")

phase_file <- file.path(output_root, paste0("pf_phase_profile_", canonical_label, ".csv"))
top_file <- file.path(output_root, paste0("pf_top_level_timing_", canonical_label, ".csv"))
rprof_file <- file.path(output_root, paste0("pf_rprof_", canonical_label, ".out"))
summary_file <- file.path(output_root, paste0("pf_profile_diagnostic_summary_", canonical_label, ".csv"))
notes_file <- file.path(output_root, paste0("pf_profile_diagnostic_report_", canonical_label, ".txt"))

if (!file.exists(phase_file)) stop("Phase profile file not found: ", phase_file)
if (!file.exists(top_file)) stop("Top-level timing file not found: ", top_file)

phase <- utils::read.csv(phase_file, stringsAsFactors = FALSE)
top <- utils::read.csv(top_file, stringsAsFactors = FALSE)

phase_columns <- intersect(c(
  "proposal_lookup",
  "propagate_particles",
  "physical_weights",
  "directional_weights",
  "smart_filter",
  "outlier_check",
  "cumulative_weight_product",
  "ESS_calculation",
  "resampling",
  "accepted_particles_assignment",
  "results_stack_append",
  "weights_stack_append",
  "point_rle_creation",
  "transition_rle_creation",
  "stack_drop",
  "final_smoothing",
  "final_stack_flush"
), names(phase))

phase_totals <- data.frame(
  phase = phase_columns,
  elapsed_seconds = vapply(phase[phase_columns], function(x) sum(x, na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)
phase_totals <- phase_totals[order(phase_totals$elapsed_seconds, decreasing = TRUE), ]
total_loop_seconds <- sum(phase$total_iteration_time, na.rm = TRUE)
phase_totals$fraction_of_profiled_loop <- if (total_loop_seconds > 0) phase_totals$elapsed_seconds / total_loop_seconds else NA_real_

utils::write.csv(phase_totals, summary_file, row.names = FALSE)

resampling_rows <- phase[phase$did_resample %in% TRUE, , drop = FALSE]
non_resampling_rows <- phase[!(phase$did_resample %in% TRUE) & !is.na(phase$Time.Period), , drop = FALSE]

cor_line <- function(xname) {
  if (!xname %in% names(phase)) return(paste0(xname, ": not available"))
  ok <- is.finite(phase[[xname]]) & is.finite(phase$total_iteration_time)
  if (sum(ok) < 3) return(paste0(xname, ": insufficient data"))
  paste0(xname, ": cor=", signif(stats::cor(phase[[xname]][ok], phase$total_iteration_time[ok]), 4))
}

rprof_total <- NULL
rprof_self <- NULL
if (file.exists(rprof_file)) {
  rprof_summary <- summaryRprof(rprof_file, memory = "both")
  rprof_total <- head(rprof_summary$by.total, 15)
  rprof_self <- head(rprof_summary$by.self, 15)
}

notes <- c(
  paste0("Particle-filter diagnostic report: ", canonical_label),
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "Top-level timings:",
  paste(capture.output(print(top)), collapse = "\n"),
  "",
  "Top phase totals:",
  paste(capture.output(print(head(phase_totals, 15))), collapse = "\n"),
  "",
  paste0("Profile rows=", nrow(phase)),
  paste0("Total profiled loop/final seconds=", signif(total_loop_seconds, 6)),
  paste0("Resampling events=", nrow(resampling_rows)),
  paste0("Median iteration seconds with resampling=", if (nrow(resampling_rows)) signif(stats::median(resampling_rows$total_iteration_time, na.rm = TRUE), 6) else NA),
  paste0("Median iteration seconds without resampling=", if (nrow(non_resampling_rows)) signif(stats::median(non_resampling_rows$total_iteration_time, na.rm = TRUE), 6) else NA),
  "",
  "Iteration-time correlations:",
  cor_line("n_unique_last_particles"),
  cor_line("n_groups"),
  cor_line("n_moving_groups"),
  cor_line("n_unique_new_particles"),
  cor_line("ESS"),
  cor_line("Results.stack_ncol"),
  cor_line("Weights.stack_ncol"),
  cor_line("n_unique_point_rle_values"),
  cor_line("n_unique_transition_rle_values"),
  "",
  "Top Rprof functions by total time:",
  if (is.null(rprof_total)) "Rprof output not found." else paste(capture.output(print(rprof_total)), collapse = "\n"),
  "",
  "Top Rprof functions by self time:",
  if (is.null(rprof_self)) "Rprof output not found." else paste(capture.output(print(rprof_self)), collapse = "\n"),
  "",
  "Initial recommendation rule:",
  "Choose the largest phase that is model-neutral and particle-history-safe. Do not implement from this report without reviewing the classified phase evidence.",
  "Raw private GPS/GLS rows are not read or printed by this summarizer."
)

writeLines(notes, notes_file)

message("Phase summary saved to: ", summary_file)
message("Diagnostic report saved to: ", notes_file)
message("Top phases:")
print(head(phase_totals, 10))

invisible(list(
  phase_summary = phase_totals,
  top_level = top,
  summary_file = summary_file,
  notes_file = notes_file
))
