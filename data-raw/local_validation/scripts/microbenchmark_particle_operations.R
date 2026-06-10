# Synthetic local microbenchmarks for particle-count-scaled operations.
# Does not read private data and does not modify package code.

output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "microbenchmark_particle_operations.csv")

n_grid <- as.integer(Sys.getenv("FLIGHTR_MICRO_N_GRID", "2591"))
nParticles <- as.numeric(strsplit(Sys.getenv("FLIGHTR_MICRO_NPARTICLES", "1e5,1e6"), ",")[[1]])
nParticles <- as.integer(nParticles)
n_reps <- as.integer(Sys.getenv("FLIGHTR_MICRO_REPS", "3"))
n_stack_cols <- as.integer(Sys.getenv("FLIGHTR_MICRO_STACK_COLS", "50"))
random_seed <- as.integer(Sys.getenv("FLIGHTR_MICRO_SEED", "123"))

set.seed(random_seed)

transition_base <- function(n_grid) n_grid + 1L

encode_transition <- function(from, to, n_grid) {
  from * transition_base(n_grid) + to
}

point_rle_old <- function(points) {
  rle(sort.int(points, method = "quick"))
}

point_rle_tabulate <- function(points, n_grid) {
  counts <- tabulate(points, nbins = n_grid)
  values <- which(counts > 0L)
  structure(list(lengths = as.integer(counts[values]), values = as.integer(values)), class = "rle")
}

transition_rle_old <- function(keys) {
  rle(sort.int(keys, method = "quick"))
}

transition_rle_radix <- function(keys) {
  rle(sort.int(keys, method = "radix"))
}

transition_rle_tabulate_dense <- function(keys, n_grid) {
  nbins <- n_grid * transition_base(n_grid) + n_grid
  counts <- tabulate(keys, nbins = nbins)
  values <- which(counts > 0L)
  structure(list(lengths = as.integer(counts[values]), values = as.integer(values)), class = "rle")
}

build_stack_cbind <- function(columns) {
  out <- NULL
  for (col in columns) out <- cbind(out, col)
  out
}

build_stack_prealloc <- function(columns) {
  out <- matrix(NA_real_, nrow = length(columns[[1]]), ncol = length(columns))
  for (i in seq_along(columns)) out[, i] <- columns[[i]]
  out
}

build_stack_list_docbind <- function(columns) {
  do.call(cbind, columns)
}

weighted_mean_rle <- function(values, lengths) {
  sum(values * lengths) / sum(lengths)
}

weighted_quantile_type1 <- function(values, lengths, probs = c(0.25, 0.5, 0.75, 0.95)) {
  ord <- order(values)
  values <- values[ord]
  lengths <- lengths[ord]
  cum <- cumsum(lengths)
  n <- sum(lengths)
  out <- vapply(probs, function(p) values[which(cum >= ceiling(p * n))[1]], numeric(1))
  names(out) <- paste0(probs * 100, "%")
  out
}

summarise_expanded <- function(rle_obj) {
  x <- inverse.rle(rle_obj)
  c(mean = mean(x), stats::quantile(x, probs = c(0.25, 0.5, 0.75, 0.95), type = 1))
}

summarise_weighted <- function(rle_obj) {
  c(
    mean = weighted_mean_rle(rle_obj$values, rle_obj$lengths),
    weighted_quantile_type1(rle_obj$values, rle_obj$lengths)
  )
}

rle_counts_equal <- function(a, b) {
  identical(as.integer(a$values), as.integer(b$values)) &&
    identical(as.integer(a$lengths), as.integer(b$lengths))
}

result_rows <- list()

append_result <- function(operation, candidate, n_particles, rep, elapsed, correct, notes = "") {
  row <- data.frame(
    operation = operation,
    candidate = candidate,
    n_grid = n_grid,
    nParticles = n_particles,
    rep = rep,
    elapsed_seconds = elapsed,
    correct = correct,
    notes = notes,
    stringsAsFactors = FALSE
  )
  result_rows[[length(result_rows) + 1L]] <<- row
  utils::write.csv(do.call(rbind, result_rows), output_path, row.names = FALSE)
}

time_one <- function(operation, candidate, n_particles, rep, expr, correct_fun, notes = "") {
  gc()
  value <- NULL
  err <- NULL
  elapsed <- system.time({
    value <- tryCatch(expr, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    })
  })[["elapsed"]]
  correct <- FALSE
  if (is.null(err)) {
    correct <- tryCatch(isTRUE(correct_fun(value)), error = function(e) FALSE)
  }
  append_result(
    operation = operation,
    candidate = candidate,
    n_particles = n_particles,
    rep = rep,
    elapsed = as.numeric(elapsed),
    correct = correct,
    notes = if (is.null(err)) notes else paste0("error: ", err)
  )
  value
}

message("Writing microbenchmark results to: ", output_path)
message("Synthetic only: n_grid=", n_grid,
        " nParticles=", paste(nParticles, collapse = ","),
        " reps=", n_reps,
        " stack_cols=", n_stack_cols)

for (n in nParticles) {
  message("nParticles=", n)
  for (rep in seq_len(n_reps)) {
    message("  rep=", rep)

    points <- sample.int(n_grid, n, replace = TRUE)
    point_ref <- NULL
    point_ref <- time_one(
      "point_marginal_rle", "old_sort_quick_rle", n, rep,
      point_rle_old(points),
      correct_fun = function(value) {
        point_ref <<- value
        TRUE
      }
    )
    invisible(time_one(
      "point_marginal_rle", "tabulate_rle", n, rep,
      point_rle_tabulate(points, n_grid),
      correct_fun = function(value) rle_counts_equal(point_ref, value)
    ))

    from <- sample.int(n_grid, n, replace = TRUE)
    to <- sample.int(n_grid, n, replace = TRUE)
    keys <- encode_transition(from, to, n_grid)
    transition_ref <- NULL
    transition_ref <- time_one(
      "transition_rle", "old_sort_quick_rle", n, rep,
      transition_rle_old(keys),
      correct_fun = function(value) {
        transition_ref <<- value
        TRUE
      }
    )
    invisible(time_one(
      "transition_rle", "old_sort_radix_rle", n, rep,
      transition_rle_radix(keys),
      correct_fun = function(value) rle_counts_equal(transition_ref, value)
    ))
    invisible(time_one(
      "transition_rle", "tabulate_dense_rle", n, rep,
      transition_rle_tabulate_dense(keys, n_grid),
      correct_fun = function(value) rle_counts_equal(transition_ref, value),
      notes = "Dense tabulate over transition key range; safe only for stored marginal transition counts, not live particle histories."
    ))

    columns <- lapply(seq_len(n_stack_cols), function(i) runif(n))
    stack_ref <- NULL
    stack_ref <- time_one(
      "weight_stack_building", "repeated_cbind", n, rep,
      build_stack_cbind(columns),
      correct_fun = function(value) {
        stack_ref <<- value
        is.matrix(value) && nrow(value) == n && ncol(value) == n_stack_cols
      },
      notes = paste0("columns=", n_stack_cols)
    )
    invisible(time_one(
      "weight_stack_building", "preallocated_matrix_assignment", n, rep,
      build_stack_prealloc(columns),
      correct_fun = function(value) isTRUE(all.equal(unname(stack_ref), unname(value), tolerance = 0)),
      notes = paste0("columns=", n_stack_cols)
    ))
    invisible(time_one(
      "weight_stack_building", "list_do_call_cbind", n, rep,
      build_stack_list_docbind(columns),
      correct_fun = function(value) isTRUE(all.equal(unname(stack_ref), unname(value), tolerance = 0)),
      notes = paste0("columns=", n_stack_cols)
    ))
    rm(columns, stack_ref)
    gc()

    summary_values <- sort(stats::rgamma(1000, shape = 2, scale = 100))
    summary_lengths <- as.integer(rmultinom(1, n, rep(1, length(summary_values)))[, 1])
    summary_values <- summary_values[summary_lengths > 0L]
    summary_lengths <- summary_lengths[summary_lengths > 0L]
    summary_rle <- structure(list(lengths = summary_lengths, values = summary_values), class = "rle")
    expanded_ref <- NULL
    expanded_ref <- time_one(
      "inverse_rle_summary", "expand_inverse_rle_mean_quantile", n, rep,
      summarise_expanded(summary_rle),
      correct_fun = function(value) {
        expanded_ref <<- value
        TRUE
      },
      notes = "quantile type=1 for exact weighted comparison"
    )
    invisible(time_one(
      "inverse_rle_summary", "weighted_mean_quantile_from_rle", n, rep,
      summarise_weighted(summary_rle),
      correct_fun = function(value) isTRUE(all.equal(unname(expanded_ref), unname(value), tolerance = 1e-12)),
      notes = "Safe for marginal summaries only; does not preserve live particle identity."
    ))
  }
}

results <- do.call(rbind, result_rows)

median_times <- aggregate(
  elapsed_seconds ~ operation + candidate + nParticles,
  results,
  median
)

baseline_candidates <- c(
  point_marginal_rle = "old_sort_quick_rle",
  transition_rle = "old_sort_quick_rle",
  weight_stack_building = "repeated_cbind",
  inverse_rle_summary = "expand_inverse_rle_mean_quantile"
)

speedups <- do.call(rbind, lapply(split(median_times, median_times$operation), function(df) {
  baseline <- baseline_candidates[[unique(df$operation)]]
  if (is.null(baseline)) return(NULL)
  do.call(rbind, lapply(split(df, df$nParticles), function(chunk) {
    base_time <- chunk$elapsed_seconds[chunk$candidate == baseline]
    if (!length(base_time)) return(NULL)
    chunk$speedup_vs_baseline <- base_time / chunk$elapsed_seconds
    chunk
  }))
}))

message("Median timing and speedup summary:")
print(speedups[order(speedups$operation, speedups$nParticles, -speedups$speedup_vs_baseline), ])

message("Results saved to: ", output_path)
