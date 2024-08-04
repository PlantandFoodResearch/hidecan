# hidecan 1.1.0.9000

* Added the `manhattan_plot()` function to generate Manhattan plots from a table of GWAS results.

* Fixed bug where specifying chromosome limits for `hidecan_plot_from_gwaspoly()` would mess up the ordering of the chromosomes in the plot.

# hidecan 1.1.0

* Removed `get_gwaspoly_example_data()` function so that the package doesn't depend on GWASpoly (for CRAN submission)

* Saved the example object previously returned by `get_gwaspoly_example_data()` in as an `rda` file in the `extdata` folder
