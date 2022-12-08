
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

The package includes functions for:

- **Converting names** from one format to another: `names_conversion()`
- **Calculating the mean across replicates**: `combine_reps()`
- **Quality control**:
  - Counting the number of proteins per sample: `protein_count()`
  - Calculating pairwise linear correlations between replicates:
    `replicate_correlations()`
  - Assessing the similarity between samples using PCA: `pca_enhanced()`
  - Assigning missingness for samples: `assign_missingness()`
  - Calculating the proportion of values that are imputed for each
    protein: `data_completeness()`
- **Fuzzy clustering**:
  - Preparing data for fuzzy clustering: `mfuzz_prep()`
  - Plotting the output of fuzzy clustering: `plot_fuzzy_clusters()`
  - Extracting cluster cores from fuzzy clustering:
    `alpha_core_modified()`
- **Gene Ontology enrichment analysis**:
  - Identifying enriched GO terms from a gene list:
    `enrichGO_enhanced()`
  - Identifying enriched GO terms in clusters:
    `clusters_enrichGO_enhanced()`
  - Plotting the output from GO analysis: `enrichGO_plot()`
  - Plotting the output from GO analysis of clusters:
    `clusters_enrichGO_plot()`
- **Calculating ribosome engagement from multiple fractions**:
  - Calculating ribosome engagement, by summing multiple ribosomal
    fractions and normalising to totals: `ribosome_engagement()`
  - Calculating differences in ribosome engagement between samples:
    `deltaR()`
  - Calculating z-scores and p-values for differences in ribosome
    engagement, relative to a reference Population: `zscore()`

## Installation

You can install the development version of proteomicshelpers from
[GitHub](https://github.com/) with:

Note: if it is not already installed, remove the comment sign (#) at the
start of the first line to install `devtools`.

``` r
# install.packages("devtools")
devtools::install_github("RobCrawford527/proteomicshelpers", dependencies = TRUE)
```

The `dependencies = TRUE` argument ensures that other packages that are
required for proteomicshelpers to work correctly are installed too.

## Converting names

The `names_conversion()` function takes an input list of gene/protein
names and converts them to the format(s) of your choice. It uses the
`bitr()` function from `clusterProfiler`, so fromType and toType must be
compatible with that. The output from the function is a data frame with
the names in different formats. Write this into an object to store for
reference.

## Calculating means across replicates

The `combine_reps()` function calculates mean intensities across the
replicates from a proteomics experiment. By setting min_reps, proteins
with missing values can be excluded. The default is `min_reps = 1`, so
means are calculated for proteins with at least one value. The function
returns a data frame in the same format as the original, but containing
only the mean values.

## Quality control

You can count the number of proteins per sample and make a plot using
`protein_count()`.

Pairwise linear correlations between replicates of the same sample can
be calculated using `replicate_correlations()`. The function returns a
plot showing the R-squared value and number of shared proteins for each.

The `pca_enhanced()` function performs principal component analysis of a
data frame and returns a labelled plot of the first two PCs, coloured as
appropriate. The plot can be saved directly.

You can assign missingness to a set of samples using
`assign_missingness()`, which is an adaptation of
\[protti::assign_missingness()\]. This compares the number of values
that are present for a given sample and the reference sample with two
thresholds, to assign samples as “complete”, “MAR” or “MNAR”. Proteins
that do not fit in any of these categories (i.e. have too few values in
the sample and reference) are removed from the analysis. The output from
`assign_missingness()` is compatible with `protti`: missing values can
be imputed for MAR and MNAR comparisons using \[protti::impute()\].

Following imputation of missing values, `data_completeness()` assesses
how what proportion of the values for each protein are real (i.e. not
imputed) and assigns each to one of five categories, depending on the
thresholds that are set. This output can be written into a new object
and merged with the imputed data frame.
