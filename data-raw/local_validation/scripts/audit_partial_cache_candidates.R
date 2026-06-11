# Local-only diagnostic: compare full cached and partial_cached candidate sets.
# Aggregate grid/candidate diagnostics only; no private GPS/GLS rows are read or printed.

.libPaths(c('C:/Users/eldar/Documents/Codex/Rlibs-4.6', .libPaths()))
repo <- 'D:/GitHub/FLightR'
out_dir <- 'D:/GitHub/FLightR/data-raw/local_validation/outputs'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(repo)

if (requireNamespace('pkgload', quietly = TRUE)) {
  pkgload::load_all(repo, quiet = TRUE, export_all = FALSE)
} else if (requireNamespace('devtools', quietly = TRUE)) {
  devtools::load_all(repo, quiet = TRUE, export_all = FALSE)
} else {
  library(FLightR)
}

env <- new.env(parent = globalenv())
sys.source('R/run_particle_filter.R', envir = env)

extent_name <- Sys.getenv('FLIGHTR_AUDIT_EXTENT_NAME', unset = 'validation_current')
grid_left <- as.numeric(Sys.getenv('FLIGHTR_GRID_LEFT', unset = '-14'))
grid_right <- as.numeric(Sys.getenv('FLIGHTR_GRID_RIGHT', unset = '13'))
grid_bottom <- as.numeric(Sys.getenv('FLIGHTR_GRID_BOTTOM', unset = '30'))
grid_top <- as.numeric(Sys.getenv('FLIGHTR_GRID_TOP', unset = '57'))
a <- as.numeric(Sys.getenv('FLIGHTR_MOVEMENT_A', unset = '45'))
b <- as.numeric(Sys.getenv('FLIGHTR_MOVEMENT_B', unset = '1500'))
all_origins <- tolower(Sys.getenv('FLIGHTR_AUDIT_ALL_ORIGINS', unset = 'true')) %in% c('1','true','yes')

Grid <- make.grid(
  left = grid_left, bottom = grid_bottom, right = grid_right, top = grid_top,
  distance.from.land.allowed.to.use = c(-Inf, Inf),
  distance.from.land.allowed.to.stay = c(-Inf, Inf),
  plot = FALSE
)

message('Building full cached candidates for ', nrow(Grid), ' grid cells')
full <- env$build.grid.movement.candidates(Grid, a = a, b = b)
cache <- env$build.partial.movement.cache(Grid, a = a, b = b)

if (all_origins || nrow(Grid) <= 5000) {
  origins <- seq_len(nrow(Grid))
} else {
  lon <- Grid[,1]
  lat <- Grid[,2]
  edge <- unique(c(
    order(lon)[seq_len(min(50, length(lon)))],
    order(lon, decreasing = TRUE)[seq_len(min(50, length(lon)))],
    order(lat)[seq_len(min(50, length(lat)))],
    order(lat, decreasing = TRUE)[seq_len(min(50, length(lat)))]
  ))
  full_counts <- vapply(full$to, length, integer(1))
  large <- order(full_counts, decreasing = TRUE)[seq_len(min(200, length(full_counts)))]
  set.seed(123)
  random <- sample(seq_len(nrow(Grid)), min(500, nrow(Grid)))
  origins <- sort(unique(c(edge, large, random)))
}

compare_one <- function(origin) {
  full_to <- full$to[[origin]]
  full_dist <- full$distance[[origin]]
  partial <- env$get.partial.movement.entry(cache, origin)
  partial_to <- partial$to
  partial_dist <- partial$distance
  missing <- setdiff(full_to, partial_to)
  extra <- setdiff(partial_to, full_to)
  miss_dist <- if (length(missing)>0) full_dist[match(missing, full_to)] else numeric()
  extra_dist <- if (length(extra)>0) partial_dist[match(extra, partial_to)] else numeric()
  data.frame(
    extent_name = extent_name,
    origin = origin,
    origin_lon = Grid[origin,1],
    origin_lat = Grid[origin,2],
    n_full = length(full_to),
    n_partial = length(partial_to),
    n_missing_from_partial = length(missing),
    n_extra_in_partial = length(extra),
    min_missing_distance_km = if (length(miss_dist)>0) min(miss_dist) else NA_real_,
    max_missing_distance_km = if (length(miss_dist)>0) max(miss_dist) else NA_real_,
    min_extra_distance_km = if (length(extra_dist)>0) min(extra_dist) else NA_real_,
    max_extra_distance_km = if (length(extra_dist)>0) max(extra_dist) else NA_real_,
    missing_near_a = if (length(miss_dist)>0) sum(abs(miss_dist-a) <= 1) else 0L,
    missing_near_b = if (length(miss_dist)>0) sum(abs(miss_dist-b) <= 1) else 0L,
    extra_near_a = if (length(extra_dist)>0) sum(abs(extra_dist-a) <= 1) else 0L,
    extra_near_b = if (length(extra_dist)>0) sum(abs(extra_dist-b) <= 1) else 0L,
    west_edge = Grid[origin,1] <= quantile(Grid[,1], 0.05),
    east_edge = Grid[origin,1] >= quantile(Grid[,1], 0.95),
    south_edge = Grid[origin,2] <= quantile(Grid[,2], 0.05),
    north_edge = Grid[origin,2] >= quantile(Grid[,2], 0.95),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(origins, compare_one))
runtime_stamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
out_csv <- file.path(out_dir, 'partial_cache_candidate_audit.csv')
out_csv_labeled <- file.path(out_dir, paste0('partial_cache_candidate_audit_', extent_name, '_', runtime_stamp, '.csv'))
utils::write.csv(audit, out_csv, row.names = FALSE)
utils::write.csv(audit, out_csv_labeled, row.names = FALSE)

mismatches <- audit[audit$n_missing_from_partial > 0 | audit$n_extra_in_partial > 0,]
summary <- data.frame(
  extent_name = extent_name,
  grid_size = nrow(Grid),
  origins_checked = length(origins),
  exact_matches = sum(audit$n_missing_from_partial == 0 & audit$n_extra_in_partial == 0),
  mismatches = nrow(mismatches),
  total_missing_candidate_destinations = sum(audit$n_missing_from_partial),
  total_extra_candidate_destinations = sum(audit$n_extra_in_partial),
  missing_near_a = sum(audit$missing_near_a),
  missing_near_b = sum(audit$missing_near_b),
  extra_near_a = sum(audit$extra_near_a),
  extra_near_b = sum(audit$extra_near_b),
  cached_origins = env$partial.movement.cache.diagnostics(cache)$cached_origins,
  cache_hits = env$partial.movement.cache.diagnostics(cache)$cache_hits,
  cache_misses = env$partial.movement.cache.diagnostics(cache)$cache_misses,
  stringsAsFactors = FALSE
)
summary_csv <- file.path(out_dir, 'partial_cache_candidate_audit_summary.csv')
utils::write.csv(summary, summary_csv, row.names = FALSE)

top_mismatch <- if (nrow(mismatches)>0) {
  mismatches[order(mismatches$n_missing_from_partial, decreasing = TRUE),
             c('origin','origin_lon','origin_lat','n_full','n_partial','n_missing_from_partial','n_extra_in_partial','min_missing_distance_km','max_missing_distance_km','missing_near_a','missing_near_b','west_edge','east_edge','south_edge','north_edge')]
} else mismatches
if (nrow(top_mismatch)>20) top_mismatch <- top_mismatch[seq_len(20),]

notes <- c(
  paste('Created:', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')),
  paste('Extent:', extent_name),
  paste('Grid bounds:', paste(c(left=grid_left, bottom=grid_bottom, right=grid_right, top=grid_top), collapse=', ')),
  paste('Movement bounds a/b km:', a, b),
  paste('Grid size:', nrow(Grid)),
  paste('Origins checked:', length(origins)),
  paste('Exact matches:', summary$exact_matches),
  paste('Mismatches:', summary$mismatches),
  paste('Total missing from partial:', summary$total_missing_candidate_destinations),
  paste('Total extra in partial:', summary$total_extra_candidate_destinations),
  paste('Missing near a:', summary$missing_near_a),
  paste('Missing near b:', summary$missing_near_b),
  '',
  'Top mismatching origins by missing count:',
  paste(capture.output(print(top_mismatch, row.names = FALSE)), collapse = '\n'),
  '',
  'Raw private GPS/GLS rows were not read or printed.'
)
notes_path <- file.path(out_dir, 'partial_cache_candidate_audit_notes.txt')
writeLines(notes, notes_path)

print(summary)
print(top_mismatch, row.names = FALSE)
message('Audit saved to: ', out_csv)
message('Audit notes saved to: ', notes_path)