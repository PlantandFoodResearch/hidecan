# hidecan 1.2.0.0

### Added option to display QTL regions

* New S3 classes to store QTL regions: `QTL_data` and `QTL_data_thr`.

* `hidecan_plot()` now has arguments `qtl_list` (to pass a list of QTL mapping results) and `score_thr_qtl` (to set the significance threshold on the score of QTL regions). The significant QTL regions are represented as rectangles spanning the length of the region on the tracks.

### Changes to enable custom tracks

* New S3 classes to store custom genomic features: `CUSTOM_data` and `CUSTOM_data_thr`, which can be represented either as genomic locations (points) or regions (rectangles).

* `create_hidecan_plot()` now accepts `CUSTOM_data_thr` objects as input, which allows users to add custom tracks to their plot. It has an additional argument `custom_aes` to handle the aesthetics for custom tracks.

* Added the `hidecan_aes()` function which returns the default aesthetics for the different types of data in the plot.

### Other changes

* Added the `manhattan_plot()` function to generate Manhattan plots from a table of GWAS results.

* Fixed bug where specifying chromosome limits for `hidecan_plot_from_gwaspoly()` would mess up the ordering of the chromosomes in the plot.

# hidecan 1.1.0

* Removed `get_gwaspoly_example_data()` function so that the package doesn't depend on GWASpoly (for CRAN submission)

* Saved the example object previously returned by `get_gwaspoly_example_data()` in as an `rda` file in the `extdata` folder
