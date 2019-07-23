
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Combining serological and contact data to derive target immunity levels for achieving and maintaining measles elimination

This repository contains the data and code for our paper:

> Sebastian Funk, Jennifer K. Knapp, Emmaculate Lebo, Susan E. Reef,
> Alya J. Dabbagh, Katrina Kretsinger, Mark Jit, W. John Edmunds and
> Peter M. Strebel (2019). *Combining serological and contact data to
> derive target immunity levels for achieving and maintaining measles
> elimination*. bioRxiv, online at <https://doi.org/10.1101/201574>

### How to download or install

You can download the compendium as a zip from from this URL:
<https://github.com/sbfnk/immunity.thresholds/archive/master.zip>

Or you can install this compendium as an R package,
`immunity.thresholds`, from GitHub with:

``` r
# install.packages("devtools")
remotes::install_github("sbfnk/immunity.thresholds")
```

### Included data sets

The package includes a data set of national seroprevalence levels found
in the ESEN2 study as reported by [Andrews et
al. (2008)](https://www.who.int/bulletin/volumes/86/3/07-041129/en/).

``` r
data(esen2_serology)
```

The `esen2_serology` data set is a data frame of four columns,
`lower.age.limit` (the lower limit of each represented age group),
`country` (the country in which the study was conducted), `survey.yr`
(the year(s)) in which the study was conducted), `eqi_pos` (the
proportion of samples equivocal or positive) and `neg` (the proportion
of samples negative). For the manuscript, we sampled from a more
comprehensive, individualised data set which is available upon request
(see “Data availability” statement in the manuscript for contacts). A
toy data set of this type can be created with the `toy_sero_ds` function
included in the package:

``` r
ms.sero <- toy_sero_ds(survey.yr=2002, pos_eqi=list(Germany=c(0.8, 0.1), Finland=c(0.9, 0)), lower.age.limits=seq(0, 70, by=5), n=100)
```

Using the toy data set or the original data set, once can generate a
results on adjusted immunity levels from the different models used in
the manuscript, by running the `serology_sample.r` script:

``` r
source("analysis/serology_sample.r")
```

This generates a variable called `sero_adjImm` containing the results.
The results form the paper using the original data set are contained in
the `esen2_adjImm` data set in the package:

``` r
data(esen2_adjImm)
sero_adjImm <- esen2_adjImm
```

To generate scenarios of adjusted immunity, one can run

``` r
source("analysis/scenarios_sample.r")
```

This creates a variable called `scenarios_adjImm` containing the results
for all countries in the paper except the South-East Asian countries
from the SMILI study, for which contact data are available upon request
(see “Data availability” statement in the manuscript for contacts). The
full results including all contact data are included in the package and
can be loaded with

``` r
data(scenarios_adjImm)
```

### Table and figures

To generate Figures 1 and 2 and Tables 1 and 2 from the paper from a
`sero_adjImm` variable (see above), one can run

``` r
source("analysis/serology_analysis.r")
```

To generate Figure 3 from the paper from a `scenarios_adjImm` variable
(see above), one can run

``` r
source("analysis/scenarios_analysis.r")
```
