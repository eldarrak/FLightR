# Local Prerun Cache

This cache is only for local validation and benchmarking scripts. It is not package user-facing behavior.

`make.prerun.object()` can dominate repeated benchmark runtime, especially for full-year or large-grid runs. The cache stores the resulting prerun object under:

`data-raw/local_validation/outputs/prerun_cache/`

These files can contain derived information from private validation data and must not be committed.

## Enable reuse

Set:

```r
Sys.setenv(FLIGHTR_REUSE_PRERUN = "1")
```

or in PowerShell:

```powershell
$env:FLIGHTR_REUSE_PRERUN = "1"
```

When enabled, scripts load a matching prerun object if the cache key and basic dimensions match.

## Force rebuild

Set:

```r
Sys.setenv(FLIGHTR_FORCE_REBUILD_PRERUN = "1")
```

This rebuilds and overwrites the matching cache file.

## Save without reuse

Set:

```r
Sys.setenv(FLIGHTR_SAVE_PRERUN = "1")
```

This builds normally but saves the prerun object for a later reuse-enabled run.

## Cache key contents

The key includes public/aggregate setup information such as input file path, size, mtime and MD5, subset dates, calibration periods/settings, likelihood correction, grid bounds and grid size, land-distance settings, start/stop coordinates, make.prerun arguments, FLightR commit, dirty working tree flag, and package version.

Raw GPS/GLS rows are not written into metadata.

## Backend note

In the current speed-distance-kernel branch, cached propagation candidates are built inside `run.particle.filter()` / `pf.run.parallel.SO.resample()`, not inside `make.prerun.object()`. Therefore the prerun cache is backend-independent. If future code stores backend-specific structures in the prerun object, backend/cache mode must be added to the cache key.
