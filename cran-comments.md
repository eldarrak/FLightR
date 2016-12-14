## Resubmission
This is a resubmission. In this version I have:

1. Possibly mis-spelled words in DESCRIPTION:
  postions (9:59)
changed description to avoid use of word package. Hope this is what was needed.
  
2. * checking package dependencies ... NOTE
Package in Depends/Imports which should probably only be in LinkingTo: ‘RcppArmadillo’

 I use RcppArmadillo::fastLmPure in my package, so RcppArmadillo has to be in Imports but not LinkingTo

3. and the running the examples fails for me with...

Fixed, sorry about that!

## Test environments
* ubuntu 12.04 (on travis-ci), R 3.3.2
* win-builder (devel and release)
---