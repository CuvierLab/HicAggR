# HicAggR 1.1.2

* minor fix in doc for preparePlotgardener

# HicAggR 1.1.1

* HicAggR is now on bioconductor release 3_19
* added a minor fix to ExtractSubmatrix when using GRanges
object to extract, retain only GRanges with seqnames that
are present in hicList
* added preparePlotgardener function to generate a data.frame
that can be directly used in plotgardener's plotHicTriangle,
plotHicRectangle and plotHicSquare directly

# HicAggR 0.99.7

* caching in vignette

# HicAggR 0.99.6

* default set column names in OrientateMatrix are now removed,
this is for ggAPA
* ggAPA: "rf" mode has it's own customized axes labels now
* added introduction to the package's vignette

# HicAggR 0.99.5

* added option to remove duplicated submatrices in SearchPairs

# HicAggR 0.99.4

* removed CITATION file

# HicAggR 0.99.3

* NAMESPACE was generated with roxygen2 to define exportable functions.
* CompareToBackground: to correct the skewedness of o/e values towards long distances,
    computation of z.scores is now calculated using residuals from a polynomial
    model that fits the background couples (log(counts)~distance).
* SearchPairs: added option to remove self interacting bins.
* all internal functions are now in utilities.R.
* Docs were reviewed.

# HicAggR 0.99.2

* Implemented import of corrected matrices for data in .hic and cool/mcool format.
* Implemented import of O/E matrix for data in .hic format.
* Implemented import of raw data in .h5 format.
* Added GetInfo function to get info on a hic data (.hic, cool/mcool/h5 formats).
* Removed dependency to BSDA::z.test in CompareToBackground.
* Removed dependency to InteractionSet and added it as package to import in NAMESPACE to remove the PackageStartUpMessages.
* Removed chatty package start up message and replaced it with nicer message.
* Some BiocCheck NOTES were also takedn into consideration: changing sapply to vapply etc.

# HicAggR 0.99.1

* Corrected with Bioconductor's reviews
* Added PrepareMtxList, ImportLoops, plotMultiAPA & CompareToBackground functions
* Corrected bugs on seqlevels consistancy and name column for GRanges objects
* Encapsulated small and internal functions in utilities.R
* ExtractSubMatrix has option to remove duplicated submatrices
* Corrected quantilization operations in Aggregation
* Corrected over all code with suggestions from BiocCheck

# HicAggR 0.99.0

* Submitted to Bioconductor