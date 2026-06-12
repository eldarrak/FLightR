## Test environments

* local Windows 11 x64, R 4.6.0 ucrt
* GitHub Actions: macOS, Windows, Ubuntu, R-release/devel/oldrel
* AppVeyor: Windows, R-release
* R-hub v2 checks: pending/submitted/passed

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

This release includes a faster particle-filter propagation backend and quieter
default particle-filter output. Local validation benchmarks showed around
20x faster particle-filter runtime on the validation dataset, with unchanged
checked posterior summaries.

`RcppArmadillo` is listed in Imports because FLightR calls
`RcppArmadillo::fastLmPure()` from R code, not only through compiled code.