FLightR
=======

here is a source code for R package FLightR
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
Examples of packages use:

1.  [tree swallow example with BAS tag](https://github.com/eldarrak/FLightR/blob/master/examples/tree_swallow_BAS_tag_example/tree_swallow_analysis.Rmd)
2.  [black-tailed godwit example with Intigeo tag](https://github.com/eldarrak/FLightR/blob/master/examples/black-tailed%20godwit_Intigeo_tag_example/godwit_intigeo_analysis.Rmd)

Do not know what the difference between BAS and Intigeo is? In short intigeo are being currently produced by James Fox, and measure data up to very high sun elevation angles, BAS tags were produced till ~2013 (?) and measured data at the low sun angles (with maximum written at 64). More on the tag specific differences can be found [here](https://github.com/eldarrak/FLightR/wiki/log.light-and-log.irrad-boundaries).

-------------

Help files and vignette for the package are still to come :(

-------------

 I am also slowly developing [wiki pages for the package] (https://github.com/eldarrak/FLightR/wiki)

