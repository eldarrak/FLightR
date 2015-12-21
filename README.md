FLightR
=======

Welcome to FLightR, an R package that deals with solar geolocation data. 
The scope of the package is to position animal using data from solar geolocation archival tags
the package is in the active development.

---------------	
to install the package try:
```{r}
    require("devtools")
    install_github("eldarrak/FLightR")
	library(FLightR)
```
---------------
## NB
Very important changes were made for version 0.3.6
Just contact me if results come out strange.

##Examples of packages use:

1.  [tree swallow example with BAS tag](https://github.com/eldarrak/FLightR/blob/master/examples/tree_swallow_BAS_tag_example/tree_swallow_analysis.Rmd)
2.  [black-tailed godwit example with Intigeo tag](https://github.com/eldarrak/FLightR/blob/master/examples/Black-Tailed_Godwit_JAB_example/A6_FLightR_analysis.Rmd)

Do not know what the difference between BAS and Intigeo is? In short intigeo are being currently produced by Migrate Technology Ltd, and measure data up to very high sun elevation angles, BAS tags are the old ones produced inititally by British Antarctic Survey, then by Migrate Technology (till ~2013) and Lotek (still available) and measured data at the low sun angles (with maximum written at 64). More on the tag specific differences can be found [here](https://github.com/eldarrak/FLightR/wiki/setting-up-tag-specific-boundaries).

-------------

Help files and vignette for the package are still to come :(

-------------

I am also slowly developing [wiki pages for the package] (https://github.com/eldarrak/FLightR/wiki)

