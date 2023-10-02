[![R-CMD-check](https://github.com/eldarrak/FLightR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eldarrak/FLightR/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/eldarrak/FLightR/master.svg)](https://app.codecov.io/github/eldarrak/FLightR?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/FLightR)](https://CRAN.R-project.org/package=FLightR)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/FLightR?color=brightgreen)](https://CRAN.R-project.org/package=FLightR)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/FLightR?color=brightgreen)](https://CRAN.R-project.org/package=FLightR)


FLightR
=======

Welcome to FLightR, an R package that deals with solar geolocation data. 
The scope of the package is to position animal using data from solar geolocation archival tags
the package is in the active development.

---------------	
```{r}
    install.packages("FLightR")
```
I update CRAN version several times a year, but to get the latest version from GitHub try:
```{r}
    devtools::install_github("eldarrak/FLightR")
```

```{r}
   library(FLightR)
```

---------------

## NB
Version 0.4.4 All function names started from `plot.` are replaced wit `plot_`, e.g. `plot.lon.lat()` became `plot_lon_lat()`
Version 0.3.9 has got two major changes:

1. Workflow was completely rewritten and simplified. Updated workflow is [here](https://github.com/eldarrak/FLightR/blob/master/examples/Black-Tailed_Godwit_FLightR_vignette/FLightR_analysis_workflow.Rmd)
2. New version does not require lots of RAM, but it became slower.

Very important changes were made for version 0.3.6
Just contact me if results come out strange.

## Examples of packages use:
1. [new (>=0.4.7) workflow for red-backed shrike example with Intigeo tag](https://github.com/eldarrak/FLightR_shrike_vignette)
2.  [new (>=0.3.9) workflow for black-tailed godwit example with Intigeo tag](https://github.com/eldarrak/FLightR/blob/master/examples/Black-Tailed_Godwit_FLightR_vignette/FLightR_analysis_workflow.Rmd)
3.  [old (<0.3.9) workflow for tree swallow example with BAS tag](https://github.com/eldarrak/FLightR/blob/master/examples/tree_swallow_BAS_tag_example/tree_swallow_analysis.Rmd)
4.  [old (<0.3.9) workflow for black-tailed godwit example with Intigeo tag](https://github.com/eldarrak/FLightR/blob/master/examples/Black-Tailed_Godwit_JAB_example/A6_FLightR_analysis.Rmd)

Do not know what the difference between BAS and Intigeo is? In short intigeo are being currently produced by Migrate Technology Ltd, and measure data up to very high sun elevation angles, BAS tags are the old ones produced initially by British Antarctic Survey, then by Migrate Technology (till ~2013) and Lotek (still available) and measured data at the low sun angles (with maximum written at 64). More on the tag specific differences can be found [here](https://github.com/eldarrak/FLightR/wiki/setting-up-tag-specific-boundaries).

-------------

Vignette for the package is available [here](https://github.com/eldarrak/FLightR/blob/master/vignettes/FLightR_with_black-tailed_godwit_vignette_from_MEE_2017.Rmd)

-------------

## References
1. Rakhimberdiev, E., Winkler, D.W., Bridge, E., Seavy, N.E., Sheldon, D., Piersma, T. & Saveliev, A. (2015). A hidden Markov model for reconstructing animal paths from solar geolocation loggers using templates for light intensity. Movement Ecology, 3, 25. [Check it](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-015-0062-5)

2. Rakhimberdiev, E., Senner, N.R., Verhoeven, M.A., Winkler, D.W., Bouten, W. & Piersma, T. (2016). Comparing inferences of solar geolocation data against high-precision GPS data: annual movements of a double-tagged black-tailed godwit. Journal of Avian Biology, 47, 589–596. [Check it](https://onlinelibrary.wiley.com/doi/10.1111/jav.00891)

3. Rakhimberdiev, E. (2016). Ornithology by light levels today: dealing with a developing teenager. Wader Study, 123, 1–3. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.166000.svg)](https://doi.org/10.5281/zenodo.166000)

4. Rakhimberdiev, E.Saveliev A., Piersma, T. & Karagicheva, J. (2017). FLightR: An R package for reconstructing animal paths from solar geolocation loggers. Methods in Ecology and evolution. [Check it](https://doi.org/10.1111/2041-210X.12765)

-------------
Discussion web forum for solar geolocation is available at [ornithologyexchange](https://ornithologyexchange.org/forums/forum/259-geolocator-discussion-support/). Ask there if you need help with the FLightR _per se_ or with solar geolocation in general.

