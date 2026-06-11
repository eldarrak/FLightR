# Partial/Lazy Candidate Cache Design Notes

Date: 2026-06-10
Branch context: speed-distance-kernel

## Problem

The current cached propagation backend is excellent on the local validation grid with 2591 cells, but benchmark evidence shows that its advantage disappears on a larger realistic grid with 11821 cells. The dense/full candidate construction has an n_grid^2 memory and setup-time component, so it is not a safe universal default.

Observed evidence:

- 2591-cell validation grid: full cached backend is much faster than legacy propagation.
- 11821-cell larger grid: full cached backend is about equal to legacy.
- Estimated dense full-cache memory: about 128 MB for 2591 cells and about 2.7 GB for 11821 cells.
- Particle-filter PSOCK threads are slow in the current benchmark setup; preparation/prerun parallelism should remain separate and can still be useful.
- Local prerun-object caching is available and should be used for repeated validation benchmarks.

## Current Full-Cached Backend

The current backend builds movement candidates eagerly in `build.grid.movement.candidates()` before particle propagation starts. For every origin grid cell it stores:

- integer destination indices within movement bounds `[a, b]`
- distance in km for those destinations
- bearing in degrees from origin/current cell to destination cell
- cache metadata: `a`, `b`, and `grid_n`

Although the stored candidate lists are sparse after filtering by `[a, b]`, construction still computes `sf::st_distance()` from every origin to the entire grid. This is the source of the n_grid^2 setup cost. Bearings are computed only for retained candidates.

Backend selection currently accepts `propagation.backend = c("auto", "cached", "legacy")`. `auto` attempts the same eager cached construction and falls back to legacy only if construction fails. It does not yet estimate memory, choose a partial cache, or choose legacy based on grid size.

## Proposed Backend: Partial/Lazy Candidate Cache

Add an explicit experimental backend, for example `propagation.backend = "partial_cached"`, without changing defaults initially.

The core idea is to cache candidates only for origin cells actually visited by particles. Origins that are never reached do not pay distance, bearing, or memory cost.

### Candidate Structure

Use a run-local environment or list stored under `in.Data$Spatial$tmp`, for example:

```r
in.Data$Spatial$tmp$partial_movement_cache <- list(
  grid_n = nrow(Grid),
  a = a,
  b = b,
  entries = new.env(parent = emptyenv()),
  order = integer(),
  hits = 0L,
  misses = 0L,
  builds = 0L,
  evictions = 0L,
  max_cache_mb = 512,
  max_cached_origins = NA_integer_
)
```

Each origin entry should contain:

```r
list(
  from = from_index,
  to = integer_destination_indices,
  distance = numeric_distance_km,
  bearing = numeric_bearing_degrees,
  n_candidates = length(to),
  object_size_bytes = as.numeric(object.size(entry))
)
```

Store bearings because directional proposals may be used at any twilight. If memory becomes an issue, a later refinement could store bearings only after the first `Kappa > 0` request, but that adds complexity and should be benchmarked separately.

## Candidate Construction

For one requested origin:

1. Read source longitude/latitude.
2. Apply a fast lon/lat bounding-box prefilter using movement upper bound `b`:

```r
lat_min <- lat0 - b / 111
lat_max <- lat0 + b / 111
lon_delta <- b / (111 * cos(lat0 * pi / 180))
lon_min <- lon0 - lon_delta
lon_max <- lon0 + lon_delta
```

3. Handle polar cases where `cos(lat0)` is near zero by disabling longitude filtering.
4. Handle dateline crossing by accepting `lon >= lon_min | lon <= lon_max` after wrapping longitudes to a common `[-180, 180]` convention.
5. Compute exact geodesic distances only for bounding-box candidates.
6. Keep destinations with exact distance between `a` and `b`.
7. Compute bearings only for retained destinations, using the existing corrected direction convention: origin/current cell -> destination cell.
8. If no candidate exists, keep the existing fallback behavior by returning the origin cell as the only destination.

The initial implementation should use the same distance/bearing libraries as the existing code where practical, then benchmark alternatives later. A numeric `geosphere::distGeo()` candidate build may avoid `sf` object construction for each lazy miss, but the first goal is correctness and bounded memory.

## Lazy Behavior

`generate.points.dirs()` or a helper called by it should request candidates for `from_index` from the partial cache.

- Cache hit: return existing candidate vectors.
- Cache miss: build candidates for that origin, store them if within memory limits, and return them.
- Never-visited origins are never built.
- Use FIFO or LRU eviction if the estimated cache size exceeds the configured limit.

Suggested counters:

- cached origins
- visited origins
- cache hits
- cache misses
- builds
- evictions
- median/max candidate count
- candidate fraction of full grid
- approximate cache object size

These diagnostics should be stored in the result object or local benchmark output when profiling is enabled, not printed per iteration.

## Memory Controls

Recommended controls for the experimental backend:

- `partial_cache_max_mb_per_run`, default 512 for sequential local validation.
- `partial_cache_max_origins`, optional hard cap.
- `partial_cache_eviction = c("lru", "fifo")`, default `lru` if easy, otherwise `fifo`.
- Disable or strongly warn for partial cached backend with PSOCK `threads > 1`, because each worker would build its own cache and serialization overhead may dominate.

Memory estimate per entry:

```r
entry_bytes ~= length(to) * (4 + 8 + 8) + R/list overhead
```

Use a conservative estimate such as 32 to 48 bytes per candidate when deciding whether to store an entry.

## Backend Policy

Do not change defaults until benchmarks confirm the policy. Proposed future `auto` policy:

1. Estimate dense/full cache memory before building it.
2. Use full cached when estimated dense memory is below a safe threshold and grid size is moderate.
3. Use partial/lazy cached for larger grids where full cached is too large but candidate fractions are expected to be small.
4. Use legacy if partial cache is unavailable, candidate fractions are high, or memory limits are exceeded.
5. Keep particle-filter `threads = 1` as the recommended setting unless that exact dataset/backend is benchmarked.
6. Keep prerun/preparation thread controls separate because preparation parallelism can still help.

A possible starting threshold:

- full cached if estimated dense memory <= 512 MB and `n_grid <= 5000`
- partial cached if estimated dense memory > 512 MB or `n_grid > 5000`
- legacy if partial candidate diagnostics show median candidate fraction is high enough that lazy caching gives little benefit

These numbers are placeholders for benchmarking, not final defaults.

## Correctness Tests

Focused tests should use synthetic grids only:

1. Candidate construction includes east/north nearby points and excludes points outside `[a, b]`.
2. Dateline-crossing bounding-box prefilter includes valid nearby destinations.
3. Partial cache candidate list equals a full exact scan for small grids.
4. Repeated calls to the same origin are cache hits and return identical candidates.
5. Cache eviction works with a tiny memory/origin cap.
6. `generate.points.dirs()` with `partial_cached` returns valid indices.
7. With `Kappa <= 0`, partial cached and full cached produce equivalent distance-only distributions.
8. With `Kappa > 0`, partial cached favors the intended direction.
9. Compare partial vs legacy and partial vs full cached with the same seed where the RNG path is intentionally preserved, or use distributional checks if the sampling path differs.

## Benchmark Plan

Do not start with long runs.

1. Current validation grid, 2591 cells:
   - movement-containing subset
   - nParticles = 1e4
   - threads = 1
   - compare legacy, full cached, partial cached

2. Larger realistic grid, 11821 cells:
   - same subset length/settings
   - nParticles = 1e4
   - threads = 1
   - compare legacy, full cached, partial cached

3. If partial cached is promising:
   - nParticles = 1e5 on movement-containing subset
   - then 1e6 only if memory diagnostics remain healthy

Outputs should include runtime, LL, posterior/movement comparison, cache hit/miss counters, cached-origin count, candidate counts, candidate fraction, and object size. Store only aggregate outputs under `data-raw/local_validation/outputs/`.

## Risks

- If particles eventually visit most origins and candidate fractions are high, partial caching may approach full-cache memory while adding lookup overhead.
- Different sampling implementations can alter RNG streams even with equivalent probabilities; exact equality requires preserving the full-grid `sample.int()` probability vector path.
- PSOCK parallel execution would create separate worker-local caches and may remain slower due serialization and memory duplication.
- Bounding-box filtering must handle dateline and high-latitude cases carefully.
- Backend-specific cached structures should not be stored in prerun objects unless backend/cache mode becomes part of the prerun cache key.

## Minimal Prototype Recommendation

Implement `partial_cached` only after the current staged cleanup commits are finalized. The minimal prototype should be explicit opt-in, sequential-first, and diagnostics-heavy. It should not change `auto` until it wins on the 11821-cell benchmark.
