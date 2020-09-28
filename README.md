# flowGraph

[![DOI](https://zenodo.org/badge/DOI/10.1101/837765.svg)](https://doi.org/10.1101/837765)

This repository contains and R flowGraph package.

flowGraph takes flowType Phenotypes cell populations as input and outputs SpecEnr values for each cell populations nodes (immunophenotypes), within a flow cytometry sample, based on their expected proportion. In this way, SpecEnr accounts for relations between cell populations to produce cell population quantification whose changes are not induced by neighbouring cell populations but represent real differential abundance behaviour.

## Citation

The theory, proof, and algorithm behind the SpecEnr statistic used in the flowGraph package can be found in the following [paper](https://www.biorxiv.org/content/10.1101/837765v3.abstract). Please consider citing if you found it helpful.

bibtex:
```
@article{yue2019identifying,
  title={Identifying differential cell populations in flow cytometry data accounting for marker frequency},
  author={Yue, Alice and Chauve, Cedric and Libbrecht, Maxwell and Brinkman, Ryan},
  journal={BioRxiv},
  pages={837765},
  year={2019},
  publisher={Cold Spring Harbor Laboratory}
}
```

The scripts and data from the paper can be downloaded on [Zenodo](https://zenodo.org/record/3991166).


## Installation

In R, install the package via BiocManager

```{r}
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("aya49/flowGraph")
```

## Usage

See our [vignette](vignettes/flowGraph.Rmd) for different use cases of the package: generating features, calculating summary statistics, visualizing results, and data specific use cases.

flowGraph takes as input either a (sample x cell population) raw cell count matrix OR an output from the [flowType](https://doi.org/doi:10.18129/B9.bioc.flowType) package.

flowGraph gives as output a flowGraph object that contains the SpecEnr abundance statistic and optionally its q-values. The flowGraph object can then be used to calculate more statistics, plot and visualize q-values on a cell hierarchy, etc.

```{r}
# load the package
library("flowGraph)
```
