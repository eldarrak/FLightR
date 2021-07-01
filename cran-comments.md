## RESUBMISSION
 
## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
The package was erroring in linux environments due to new r warning on time zone mismatch. I now directly define all the time zones within the package and the tests.
## Test environments
* ubuntu 14.04 (R devel, R release, R old release) on Travis
* macOS 10.13 (R release, R old release) on Travis
* win-builder (R devel, R release, R old release) Appveyor and local
---