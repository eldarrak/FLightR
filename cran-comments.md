## RESUBMISSION

## R CMD check results
There were no ERRORs or WARNINGs AND 2 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## 2 NOTE there was a note that ggmap does not have function has_google_key. It does and I this note happened because the ggmap was not update on the machine that ran the CRAN check

## Update
This is an update with the only changes related to make it compartible with ggmap>=3.0.0 

## Test environments
* ubuntu 14.04 (R devel, R release, R old release) on Travis
* macOS 10.13 (R release) on Travis
* win-builder (R devel, R release, R old release) Appveyor and local
---