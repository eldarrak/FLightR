## resubmission
this is the resubmission in which I adressed all the comments/suggestions listed below:

> Please write references in the description of the DESCRIPTION file in the form authors (year) <doi:...>. -> please write the year in parentheses.

Done


> Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. 
> \value{No return value, called for side effects} or similar) -> Missing
> Rd-tags:
      plot_lon_lat.Rd: \value

Done


> \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. 
> Please replace \dontrun with \donttest.

Done

> Please unwrap the examples if they are executable in < 5 sec, or replace dontrun{} with \donttest{}.

Done

> You write information messages to the console that cannot be easily suppressed.
> It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. 
> Instead of cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. 
> (except for print, summary, interactive functions) ->R/spatial_likelihood_estimation.R

Done - cat() is no more present in the package

> Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies. 
> Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir(). -> R/run_particle_filter.R

everything is now done with the use of the tempfiles. Also warnings are added with the file paths to ease finding these files.


## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
This update removes dependency on the archived package ggsn.

## Test environments

* ubuntu 
* macOS
* win-builder (R devel, R release, R old release) r-hub r-cmd-check on GitHub and locally
---