## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
This is an update with a few bug fixes and a new function to save results in the Movebank repository format. 

## Test environments
* ubuntu 14.04 (R devel, R release, R old release) on Travis
* macOS 10.13 (R release, R old release) on Travis
* win-builder (R devel, R release, R old release) Appveyor and local
---