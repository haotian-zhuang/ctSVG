# ctSVG

[10x Visium HD](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression) provides expression profiles for the whole transcriptome at single cell-scale resolution. The R package **ctSVG** implements a computational method to extract single-cell gene expression profiles from Visium HD data and identify cell-type-specific SVGs.

## Installation

To install the latest version of the R package from GitHub, please run following commands in R:

```
Sys.setenv("LIBARROW_MINIMAL" = FALSE)
```
Set the environment variable to install the R package **arrow** on Linux. More information is available at this [website](https://arrow.apache.org/docs/r/articles/install.html).

```
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("SparseArray")

if (!require("devtools"))
install.packages("devtools")
devtools::install_github("haotian-zhuang/findPC")
devtools::install_github("haotian-zhuang/ctSVG")
```

## Vignettes

Detailed vignettes are available at this [website](https://haotian-zhuang.github.io/ctSVG/).

## Citation

Please cite the following paper: [Submitted](https://www.biorxiv.org)

## Contact

Authors: Haotian Zhuang, Zhicheng Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Haotian Zhuang (haotian.zhuang@duke.edu)

Or open a new issue on this GitHub page
