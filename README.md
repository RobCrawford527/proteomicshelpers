
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicshelpers

<!-- badges: start -->
<!-- badges: end -->

A package containing useful helper functions for various aspects of
proteomics data analysis, including quality control, gene ontology
analysis and fuzzy clustering. Many of the functions are wrappers for
functions from other packages, including `clusterProfiler` ([Yu et al,
2012](https://doi.org/10.1089/omi.2011.0118)), `Mfuzz` ([Kumar &
Futschik, 2007](https://doi.org/10.6026%2F97320630002005)) and `protti`
([Quast et al, 2021](https://doi.org/10.1093/bioadv/vbab041)).

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
`protti::assign_missingness()`. This compares the number of values that
are present for a given sample and the reference sample with two
thresholds, to assign samples as “complete”, “MAR” or “MNAR”. Proteins
that do not fit in any of these categories (i.e. have too few values in
the sample and reference) are removed from the analysis. The output from
`assign_missingness()` is compatible with `protti`. Missing values can
then be imputed for MAR and MNAR comparisons using `protti::impute()`.

Following imputation of missing values, `data_completeness()` assesses
how what proportion of the values for each protein are real (i.e. not
imputed) and assigns each to one of five categories, depending on the
thresholds that are set. This output can be written into a new object
and merged with the imputed data frame.

## Fuzzy clustering

The package can be used for fuzzy clustering along with `Mfuzz`. The
`mfuzz_prep()` function takes a data frame (which should already be in
“wide” format, i.e. with proteins as rows and samples as columns) and
converts it to an expression set, which is required for `Mfuzz`. Three
other preparation steps are bundled in to this function: filtering for
NAs, replacing NAs and standardising the data.

Following `mfuzz_prep()`, fuzzy clustering itself is performed using the
functions in `Mfuzz`. The output can then be visualised using
`plot_fuzzy_clusters()`. This function plots profiles for proteins above
the membership threshold for each cluster, optionally showing the
background (i.e. all other proteins) and cluster centres.

Finally, cluster alpha cores can be extracted and combined into a single
data frame using `alpha_core_modified()`. Again this is really just a
wrapper for `Mfuzz::acore()` that provides the output in a different
format.

## Gene Ontology enrichment analysis

These functions utilise the package `clusterProfiler` to perform GO
analysis on lists of genes (`enrichGO_enhanced()`) or clusters
(`clusters_enrichGO_enhanced()`). There are options to simplify the
output, evaluate the terms for plotting and convert the names to a
different format. Converting the names relies on a data frame produced
using `names_conversion()`, which contains the protein names in multiple
formats.

The output from both of the `enrichGO_enhanced` functions can be plotted
using `enrichGO_plot()` or `clusters_enrichGO_plot()`, as appropriate.
`enrichGO_plot()` shows gene ratios for enriched GO terms, coloured and
sized according to adjusted p-value, while `clusters_enrichGO_plot()`
shows enriched terms in the different clusters, coloured and sized
according to adjusted p-value.

## Calculating ribosome engagement

Finally, the package contains functions to sum intensities across
multiple fractions. Specifically this is intended to calculate overall
ribosome engagement from multiple ribosome fractions, but would be
applicable to other uses too. `ribosome_engagement()` takes a data frame
in which the “fraction” aspect of the sample names has been separated
from the other parts, and sums across the defined “ribosome_fractions”.
The summed values can be normalised to the “total” fraction that is
defined.

Differences in ribosome engagement between samples can then be
calculated using `deltaR()`. The sample names must again be rearranged
to put them into the correct format for making the comparisons. For
example, if the sample names are in the format
strain_condition_fraction, they should be arranged to “strain_condition”
and “fraction” columns for `ribosome_engagement()`, then rearranged to
“strain_fraction” and “condition” columns for `deltaR()` (to compare the
different conditions). By default, all pairwise comparisons will be
performed, but you can set this to only a subset if desired.

The `zscore()` function calculates z-scores and p-values for the
differences calculated by `deltaR()`, based on a defined reference
population (e.g. the ribosomal proteins themselves). Multiple testing
correction is performed by default using the BH procedure for each
comparison.
