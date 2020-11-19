# flowGraph

[![DOI](https://zenodo.org/badge/DOI/10.1101/837765.svg)](https://doi.org/10.1101/837765)

[Summary blog post](https://aya49.github.io/2020/09/30/flowGraph/)

flowGraph is an R package used to identify candidate biomarkers for disease diagnosis in flow cytometry data. It does so by identifying driver cell populations whose abundance changes significantly and independently given a disease.

flowGraph takes cell counts as input and outputs SpecEnr values for each cell population within a flow cytometry sample, based on their expected proportion. SpecEnr accounts for dependencies between cell populations such that we can use it to flag only cell populations whose abundance change is incurred wholly or in part because of its association with a sample class (e.g. healthy vs sick).

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

**flowGraph** can be installed via Bioconductor.

You can also intall the development version directly from Github using BiocManager:

```{r}
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("aya49/flowGraph")
```

## Usage

See our [vignette](vignettes/flowGraph.Rmd) for different use cases of the package: generating features, calculating summary statistics, visualizing results, and data specific use cases.

flowGraph takes as input a (sample x cell population) raw cell count matrix.

flowGraph gives as output a flowGraph object that contains the SpecEnr abundance statistic and optionally its q-values. The flowGraph object can then be used to calculate more statistics, plot and visualize q-values on a cell hierarchy, etc.

```{r}
# load the package
library("flowGraph)
```
