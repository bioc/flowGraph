PLEASE DO NOT MAKE CHANGES TO THIS FILE.
Last saved: 2021-02-17 12:17:06

This folder contains files in a flowGraph object created by the flowGraph package.

For a given flow cytometry sample and its cell populations, flowGraph generates SpecEnr values for each cell population or immunophenotype of based their expected proportion.
By doing so, flowGraph accounts for the relation between cell populations and flags truly enriched cell populations rather than those with induced changes.

The Summary_statistics folder contains p-values obtained from features associated with each immunophenotype or/ relation between immunophenotypes. The plots folder contains plots to help interpret those p-values.

The Features folder contains features for:
- 20 flow cytometry samples.
- 4 markers: A, B, C, D
- 81 cell population nodes (nodes for short) and 216 edges on 4 layers.

Cell hierarchy plots: By default, for SpecEnr features, we generate two cell_hierarchy plots. The original one is where the colours represent difference between mean SpecEnr values across sample classes. We also add one where the colours represent difference betwee mean original values across sample classes. Original here is usually proportion: the feature used to create SpecEnr. Note that if a node is coloured lightly on the second plot, then the difference is very small, meaning the SpecEnr value may become sporadic. Therefore, when analyzing the plots, for most cases we recommend looking at most important cell populations as the ones with large difference in both plots.

