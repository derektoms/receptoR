# receptoR

This is the R code for performing receptor gene expression analysis on transcriptomics data. The principle behind **receptoR** is that by pooling large numbers of existing transcriptomics datasets, we can identify molecular receptors expressed by specific cell types, and use this information to generate hypotheses about signaling modalities that influence cell behaviour. A live implementation of is available at https://wcm.ucalgary.ca/ungrinlab/receptoR, using publicly available microarray datasets from the [GEO database](https://www.ncbi.nlm.nih.gov/geo). Using the code provided here, our analysis pipeline can be implemented on any type of data.

## Installation

``` r
devtools::install_github("derektoms/receptoR")
```

...or...

``` r
source("https://install-github.me/derektoms/receptoR")
```

---

receptoR working version of R (>= 3.5) and packages from CRAN ('dplyr', 'dbplyr', 'tidyr', 'ggplot2', 'RColorBrewer', 'readr', 'stringr', 'shiny', 'shinythemes', 'shinyjs', 'DT', 'pool', 'writexl') and from Bioconductor ('GEOmetadb', 'GEOquery', 'affy', 'limma', 'annotate', 'pheatmap', 'mixOmics', 'cowplot'). 

## Run the Shiny app

There's only one exported function in the package and it runs the Shiny app:

``` r
shinyAppDemo::launchApp()
```

---

# License
This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3

