## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
This update solves r check error on numeric input in the package version.

## Test environments

* ubuntu 
* macOS
* win-builder (R devel, R release, R old release) r-cmd-check and local
---