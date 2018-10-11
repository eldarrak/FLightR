## Resubmission
This is a resubmission. In this version I have:
 
 * Solved the Google maps api key issue. The testing environemnt now will check the functions relying on Google maps only if the api key is present as environemntal variable. I have set up the hidden environemntal api keys at travis and appveyor, to check that everything works there. CRAN check will not test these functions.
 
 * I have also added an explanation for users that they need to obtain the Google maps api key.

 *Corrected package title to title case
 
 *Extended description
 
 *Clarified roles of authors

## Note that there were 4 other notes - 
 * One on RcppArmadillo, that I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo
 
 * two more were on tests and examples that take more then 10 seconds, but they took less than 11 seconds, so I would argue this is not critical.
 
 * one false postivie on potential misspelling.
 

## Update
This is an update. 
In 0.4.7 version I have 
* fixed many small bugs

## Test environments
* ubuntu 12.04 (on travis-ci), R devel, R release, R oldrelease
* OSX
* win-builder (devel and release)
---