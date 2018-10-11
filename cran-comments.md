## Resubmission
This is a resubmission. In this version I have:
 
 * Solved the Google maps api key issue. The testing environemnt now will check the functions relying on Google maps only if the api key is present as environemntal variable. I have set up the hidden environemntal api keys at travis and appveyor, to check that everything works there. CRAN check will not test these functions.
 
 * I have also added an explanation for users that they need to obtain the Google maps api key and and an error message for the case when functions are ran with old version of ggmap.

 *Corrected package title to title case
 
 *Extended description
 
 *Clarified roles of authors

## Note that there were 3 other notes - 
 * One on RcppArmadillo, that I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo
 
 * There is note that has_goog_key() is not a function in ggmap. This is rignt, in the old ggmap that is on CRAN at the moment, there was no such a function but it exists in the new ggmap (available on GitHub). FLightR will check whether the user has the new version installed and will not run the has_goog_key with old ggmap.
 
 * one false postivie on potential misspelling.
 

## Update
This is an update. 
In 0.4.7 version I have 
* fixed many small bugs

## Test environments
* ubuntu 12.04 (R devel, R release, R old release)
* OSX (release)
* win-builder (devel and release)
---