# receptoR

This is the R code for performing receptor gene expression analysis on transcriptomics data. The principle behind **receptoR** is that by pooling large numbers of existing transcriptomics datasets, we can identify molecular receptors expressed by specific cell types, and use this information to generate hypotheses about signaling modalities that influence cell behaviour. A live implementation of is available at https://wcm.ucalgary.ca/ungrinlab/receptoR, using publicly available microarray datasets from the [GEO database](https://www.ncbi.nlm.nih.gov/geo). Using the code provided here, our analysis pipeline can be implemented on any type of data.

## Installation
receptoR is a package for the R statisical computing language, and was built and tested on versions >= 3.5, and certain dependencies did not work during testing on 3.3. To install receptoR, type either of these commands into the console:

``` r
devtools::install_github("derektoms/receptoR")
```

...or...

``` r
source("https://install-github.me/derektoms/receptoR")
```

Packages imported by receptoR come from CRAN ('dplyr', 'dbplyr', 'tidyr', 'ggplot2', 'RColorBrewer', 'readr', 'stringr', 'shiny', 'shinythemes', 'shinyjs', 'DT', 'pool', 'writexl', 'BiocManager') and from Bioconductor ('affy', 'limma', 'annotate', 'pheatmap', 'mixOmics', 'cowplot'). Occasionally, there are problems installing Bioconductor packages on installation. If this is the case, install the following packages via BiocManager:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('affy', 'limma', 'annotate', 'mixOmics')) ## these are usually the offending packages, but others should be installed as needed
````

## Run the Shiny app

There's only one exported function in the package and it runs the Shiny app:

``` r
receptoR::launchApp()
```
This will open a browser window where you can interact with the **receptoR** application. Compiled read tables can be annotated in the "Upload Transcriptome Data" tab and ultimately processed to R data files (.RDA) that can be used by the app for analysis. This is performed on the "Load Expression Datasets" tab. Currently the package comes with an RNA-seq example dataset comparing expression between mouse retina and retinal pigment epithelium using control animal data from three publicly available GEO datasets (accession: [GSE114945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114945), [GSE121858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121858), and [GSE131954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131954)).
### Uploading read tables
To search and characterize publicly available microarray data, have a look at the online version of [**receptoR**](https://wcm.ucalgary.ca/ungrinlab/receptoR). For RNA-seq or unpublished data, you are able to upload a read count table and perform the same analysis. Processing raw transcriptome data is beyong the scope of this document, but ultimately the format is a matrix of samples (columns) by features (rows). A properly formatted count table should be free of missing values, and not be normalised (i.e. raw counts). Gene symbols are the required identifier for features and should be based on human (HGNC) and mouse (MGI) nomenclature.

---

# License
This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3

