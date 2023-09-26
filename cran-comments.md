## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
This update removes dependencies from sp, maptools, rgdal and rgeos.

## Test environments

* ubuntu 
* macOS
* win-builder (R devel, R release, R old release) r-hub r-cmd-check on GitHub and locally
---