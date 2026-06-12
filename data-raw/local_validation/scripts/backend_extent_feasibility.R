# Local-only backend/grid extent feasibility benchmark for FLightR propagation backends.
# Aggregate outputs only; no raw private rows printed.

.libPaths(c('C:/Users/eldar/Documents/Codex/Rlibs-4.6', .libPaths()))
repo <- 'D:/GitHub/FLightR'
out_dir <- 'D:/GitHub/FLightR/data-raw/local_validation/outputs'
script <- 'D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(repo)

seed <- Sys.getenv('FLIGHTR_EXTENT_SEED', unset = '123')
subset_start <- Sys.getenv('FLIGHTR_EXTENT_START_DATE', unset = 'NA')
subset_end <- Sys.getenv('FLIGHTR_EXTENT_END_DATE', unset = '2013-10-15')
nParticles <- as.integer(Sys.getenv('FLIGHTR_EXTENT_NPARTICLES', unset = '10000'))
full_cache_skip_mb <- as.numeric(Sys.getenv('FLIGHTR_FULL_CACHE_SKIP_MB', unset = '4096'))
full_cache_skip_grid_size <- as.integer(Sys.getenv('FLIGHTR_FULL_CACHE_SKIP_GRID_SIZE', unset = '10000'))

extents <- data.frame(
  extent_name = c('validation_current', 'larger_realistic', 'stress_large'),
  left = c(-14, -30, -90), bottom = c(30, 20, -20), right = c(13, 40, 80), top = c(57, 70, 75),
  stringsAsFactors = FALSE
)
backends <- c('legacy', 'cached')

estimate_cache_bytes <- function(grid_size) {
  # Full cached backend stores candidate integer indices, distances, bearings, plus list overhead.
  # Worst case all destinations candidates for each origin: int + two numeric vectors.
  as.numeric(grid_size) * as.numeric(grid_size) * (4 + 8 + 8)
}

run_one <- function(ext, backend) {
  run_label <- paste0('backend_extent_', ext$extent_name, '_', backend, '_Np', nParticles, '_seed', seed)
  old <- Sys.getenv(c('FLIGHTR_RUN_LABEL','FLIGHTR_NPARTICLES','FLIGHTR_PRERUN_THREADS','FLIGHTR_PF_THREADS','FLIGHTR_THREADS','FLIGHTR_SEED','FLIGHTR_KNOWN_LAST','FLIGHTR_PROC_START_DATE','FLIGHTR_PROC_END_DATE','FLIGHTR_PROPAGATION_BACKEND','FLIGHTR_PROFILE_TOP_LEVEL','FLIGHTR_PROFILE_PHASES','FLIGHTR_VERBOSE','FLIGHTR_GRID_LEFT','FLIGHTR_GRID_RIGHT','FLIGHTR_GRID_BOTTOM','FLIGHTR_GRID_TOP'))
  on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)
  Sys.setenv(
    FLIGHTR_RUN_LABEL = run_label,
    FLIGHTR_NPARTICLES = as.character(nParticles),
    FLIGHTR_PRERUN_THREADS = '4',
    FLIGHTR_PF_THREADS = '1',
    FLIGHTR_SEED = seed,
    FLIGHTR_KNOWN_LAST = 'true',
    FLIGHTR_PROC_START_DATE = subset_start,
    FLIGHTR_PROC_END_DATE = subset_end,
    FLIGHTR_PROPAGATION_BACKEND = backend,
    FLIGHTR_PROFILE_TOP_LEVEL = 'true',
    FLIGHTR_PROFILE_PHASES = 'false',
    FLIGHTR_VERBOSE = Sys.getenv('FLIGHTR_VERBOSE', unset='quiet'),
    FLIGHTR_GRID_LEFT = as.character(ext$left),
    FLIGHTR_GRID_RIGHT = as.character(ext$right),
    FLIGHTR_GRID_BOTTOM = as.character(ext$bottom),
    FLIGHTR_GRID_TOP = as.character(ext$top)
  )
  value <- source(script, local = new.env())$value
  b <- value$benchmark
  top <- value$Result$Results$top_level_profile
  pf_core <- top$elapsed_seconds[top$phase == 'pf.run.parallel.SO.resample'][1]
  run_pf <- top$elapsed_seconds[top$phase == 'total_run.particle.filter'][1]
  data.frame(
    extent_name = ext$extent_name, left = ext$left, bottom = ext$bottom, right = ext$right, top = ext$top,
    backend = backend, cache_type = if (backend == 'cached') 'full_cached' else 'none',
    nParticles = nParticles, threads = 1L, subset_start = subset_start, subset_end = subset_end,
    grid_size = b$grid_size, n_twilights = b$n_twilights,
    movement_a = 45, movement_b = 1500,
    estimated_full_cache_bytes_worst_case = estimate_cache_bytes(b$grid_size),
    estimated_full_cache_mb_worst_case = estimate_cache_bytes(b$grid_size) / 1024^2,
    actual_cache_object_bytes = NA_real_,
    make_prerun_seconds = b$prerun_seconds,
    pf_core_seconds = pf_core,
    run_particle_filter_seconds = run_pf,
    total_benchmark_seconds = b$total_seconds,
    LL = b$LL,
    completed = TRUE,
    error_or_warning = NA_character_,
    benchmark_path = value$benchmark_path,
    result_path = value$result_path,
    stringsAsFactors = FALSE
  )
}

results <- list()
for (ei in seq_len(nrow(extents))) {
  ext <- extents[ei, ]
  for (backend in backends) {
    results[[length(results)+1]] <- tryCatch(run_one(ext, backend), error = function(e) {
      data.frame(
        extent_name = ext$extent_name, left = ext$left, bottom = ext$bottom, right = ext$right, top = ext$top,
        backend = backend, cache_type = if (backend == 'cached') 'full_cached' else 'none',
        nParticles = nParticles, threads = 1L, subset_start = subset_start, subset_end = subset_end,
        grid_size = NA_integer_, n_twilights = NA_integer_, movement_a = 45, movement_b = 1500,
        estimated_full_cache_bytes_worst_case = NA_real_, estimated_full_cache_mb_worst_case = NA_real_, actual_cache_object_bytes = NA_real_,
        make_prerun_seconds = NA_real_, pf_core_seconds = NA_real_, run_particle_filter_seconds = NA_real_, total_benchmark_seconds = NA_real_, LL = NA_real_,
        completed = FALSE, error_or_warning = conditionMessage(e), benchmark_path = NA_character_, result_path = NA_character_, stringsAsFactors = FALSE
      )
    })
    utils::write.csv(do.call(rbind, results), file.path(out_dir, 'backend_extent_feasibility_runtime.csv'), row.names = FALSE)
  }
}
runtime <- do.call(rbind, results)
utils::write.csv(runtime, file.path(out_dir, 'backend_extent_feasibility_runtime.csv'), row.names = FALSE)
memory <- runtime[, c('extent_name','backend','cache_type','grid_size','estimated_full_cache_bytes_worst_case','estimated_full_cache_mb_worst_case','actual_cache_object_bytes','completed','error_or_warning')]
utils::write.csv(memory, file.path(out_dir, 'backend_extent_feasibility_memory.csv'), row.names = FALSE)

notes <- c(
  paste('Created:', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')),
  paste('nParticles:', nParticles),
  paste('Subset start:', subset_start),
  paste('Subset end:', subset_end),
  'Available backends in current code: legacy and full cached. Partial cached backend was not found.',
  'Full cached candidates are built inside run.particle.filter()/pf.run.parallel.SO.resample(), not make.prerun.object().',
  'Worst-case full-cache memory estimate assumes all destinations are candidates for all origins; actual candidate list can be smaller due movement bounds.',
  '',
  'Runtime table:',
  paste(capture.output(print(runtime)), collapse = '\n'),
  '',
  'Raw private GPS/GLS rows were not printed. Nothing was staged or committed.'
)
writeLines(notes, file.path(out_dir, 'backend_extent_feasibility_notes.txt'))
print(runtime)
print(memory)
