## resubmission
in the submission I had two problems:

1. My old email address eldar@nioz.nl is no longer active, and I cannot send a confimration message from there.

2. there was a comment: 'Found the following (possibly) invalid file URI:
   URI: www.dx.doi.org/10.1111/2041-210X.12765
     From: 
inst/doc/FLightR_with_black-tailed_godwit_vignette_from_MEE_2017.html'

This is false positive and address to my Methods in Ecology and Evolution paper is correct. 


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