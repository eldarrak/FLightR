## RESUBMISSION
I have corrected all the suggestions by Julia Haider. See below.

Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. (If a
function does not return a value, please document that too, e.g.
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      find.stationary.location.Rd: \value
      make.calibration.Rd: \value
      map.FLightR.ggmap.Rd: \value
      plot_likelihood.Rd: \value
      plot_lon_lat.Rd: \value
      plot_slopes_by_location.Rd: \value

DONE
	  
You write information messages to the console that cannot be easily
suppressed. It is more R like to generate objects that can be used to
extract the information a user is interested in, and then print() that
object. Instead of print()/cat() rather use message()/warning() or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console. (except for print, summary, interactive functions)

DONE

Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. e.g.:
...
oldpar <- par(no.readonly = TRUE)    # code line i
on.exit(par(oldpar))            # code line i + 1
...
par(mfrow=c(2,2))            # somewhere after
...
e.g.: R/data_preparation.R,  R/new_plotting_functions.R,
R/simulations_functions.R, R/spatial_likelihood_estimation.R 

DONE

## R CMD check results
There were no ERRORs or WARNINGs AND 1 NOTES

## 1 NOTE "Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'"
This note on RcppArmadillo, I get every time: I use RcppArmadillo::fastLmPure in FLightR, so RcppArmadillo has to be in Imports but not LinkingTo

## Update
The package was erroring in linux environments due to new r warning on time zone mismatch. I now directly define all the time zones within the package and the tests.
## Test environments
* rhub
* ubuntu 14.04 (R devel, R release, R old release) on Travis
* macOS 10.13 (R release, R old release) on Travis
* win-builder (R devel, R release, R old release) Appveyor and local
---