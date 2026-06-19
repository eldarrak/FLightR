# ver 0.5.7 Jun 2026

* Improved `get.tags.data()` import controls for newer or custom loggers:
  users can now specify logger type, raw/log light scale, custom raw and
  log-light bounds, and combined
  `tag.settings`. Automatic tag/light-scale detection now records structured
  diagnostics and stops when detection is ambiguous or fails; use
  explicit settings for custom or newly supported trackers. Import
  and calibration settings are also carried in structured `Metadata` blocks
  through `Proc.data`, `Calibration`, prerun objects, and particle-filter
  results while preserving existing top-level fields.

# ver 0.5.6 Jun 2026
Release-candidate updates for the particle-filter implementation and local
validation workflow.

* Added fast cached and partial-cached movement proposal backends for
  `run.particle.filter()`, while keeping the legacy backend available for
  fallback and reproducibility.
* In local double-tag validation benchmarks, the new `partial_cached`
  propagation backend reduced particle-filter runtime by about 20-25x for a
  full-length 1e6-particle run on the tested dataset and machine. Actual
  speedups depend on grid size, backend, particle count, hardware, and thread
  settings.
* Added focused fixes for sequential `threads = 1`, start/stop point handling,
  directional bearing orientation, and transition encoding on large grids.
* Added quiet-by-default particle-filter logging via `verbose = FALSE`; use
  `verbose = TRUE` for normal progress and `verbose = "debug"` for detailed
  loop diagnostics.
* Added warnings that particle-filter PSOCK parallel execution can be slower
  than sequential execution for tested cached backends. Preparation-stage
  parallelism in `make.prerun.object()` remains available and can still be
  beneficial.
* Added local-only validation and benchmarking utilities for private
  double-tagged data. These files are excluded from CRAN builds.

# ver 0.5.5 Jun 2024
Small (and not affecting the results) bug was corrected in the summary function; removed dependcy on ggsn, as the package was archived (so scalebar cannot be plotted anymore); checked for the new ggmap; imporved unix cluster handling

# ver 0.5.4 Sept 2023
Removed dependencies from sp, maptools, rgeos and rgdal

# ver 0.5.3 Sept 2023 
Minor change in testhat part to comply with new CRAN checks.

# ver 0.5.2 Jan 2022 
Moved GeoLight from Imports to suggests. Corrected significant bug in get_utilisation_points that affected size of the estimated error distributions in plotting and summary functions!!! 

# ver 0.5.1 July 2021
Fixed minor bug in find times distribution, set time zones to gmt everywhere, replaced cat() with message(), improved documentation, conserved par values throughout the code.

# ver 0.5.0 
Small fixes in test for R 4.1.0 compartibility.

# ver 0.4.9
Small bug fixes in calibration function, FLightR2Movebank function added.

# ver 0.4.8
Made FLightR compartible with ggmap >= 3.0.0.

# ver 0.4.7 
Many small stability improvements, typo corrections and new error mesages.

# ver 0.4.6 
Added support for old swiss loggers, improved find stationary location function small bug fixes, citation update, added vignette, fixed many small bugs.

# ver 0.4.5 
Corrected mistakes in check submitted to CRAN.

# ver 0.4.4 
Many changes, all tests, helpfiles, ready for CRAN submission.

# ver. 0.4.0 
Solved dateline isssues for plotting.

# ver. 0.3.8 
added bigMemory and new workflow.

# ver. 0.3.9 
No preestimated distances and angles - attempt to save memory on large spatial extents, new workflow.

# ver 0.3.6 03 Dec 2015
New boundaries for Intigeo tags, changes in black-tailed godwit examples many bug fixes.

# ver 0.3.5 Nov 2015
Correction in time - positions match.

# ver 0.3.4 08 Oct 2014 
Added fast lm implementation.

# ver 0.2.4 June 16 2015  
pf filter became MUCH faster! the rest is more or less unchanged.

# ver 0.2.3 June 10 2015 
Added the double calibration option reorganized the workflow.

# ver 0.2.2 
Functions reorganized and the whole package is cleaned all output is in the results section now.

# ver 0.2.1 May 2015 
Many changes with reorganizing the whole thing, putting functions in thematic files etc.

# ver 0.1.6 
Brute force calibration and the whole bunch of small changes.

# ver 0.1.5 
Functions 5.0 added.

# ver 0.1.4 June 25 
`correct.hours` corrected for summer time.., also added some stability to the movies.

# ver 0.1.3 June 7 
`get.irradiance` added.

# ver 0.0.2. May 15 2013 
added median in the plotting function.

