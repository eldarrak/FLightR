## RESUBMISSION
The package was archived as it had dependencies on archived package GeoLight. These dependencies were not essential and are removed in the new verison. 

## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
This update has some bug fixes and also does not depend on the archived GeoLight anymore.

## Test environments
* rhub
* ubuntu 
* macOS
* win-builder (R devel, R release, R old release) Appveyor and local
---