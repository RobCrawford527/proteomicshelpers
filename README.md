
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicshelpers

<!-- badges: start -->
<!-- badges: end -->

A package containing useful helper functions for various aspects of
proteomics data analysis, including quality control, gene ontology
analysis and fuzzy clustering. Many of the functions are wrappers for
functions from other packages, including `clusterProfiler` ([Yu et al,
2012](https://doi.org/10.1089/omi.2011.0118)) and `Mfuzz` ([Kumar &
Futschik, 2007](https://doi.org/10.6026%2F97320630002005)).

The package includes functions for: + Converting names:
`names_conversion()` + Combining replicates: `combine_reps()` + Quality
control: `protein_count()`, `replicate_correlations()`,
`pca_enhanced()`, `assign_missingness()` and `data_completeness()` +
Fuzzy clustering: `mfuzz_prep()`, `plot_fuzzy_clusters()` and
`alpha_core_modified()` + Gene Ontology enrichment analysis:
`enrichGO_enhanced()`, `enrichGO_plot()`, `clusters_enrichGO_enhanced()`
and `clusters_enrichGO_plot()` + Calculating ribosome engagement from
multiple fractions: `ribosome_engagement()`, `deltaR()` and `zscore()`

## Installation

You can install the development version of proteomicshelpers from
[GitHub](https://github.com/) with:

Note: if it is not already installed, remove the comment sign (#) at the
start of the first line to install `devtools`.

``` r
# install.packages("devtools")
devtools::install_github("RobCrawford527/proteomicshelpers", dependencies = TRUE)
```

The dependencies = TRUE argument ensures that other packages that are
required for proteomicshelpers to work correctly are installed too.

## Quality control

This is a basic example which shows you how to solve a common problem:

``` r
library(proteomicshelpers)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
