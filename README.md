# :book: Spatial Fingerprinting - research compendium
[![R-CMD-check](https://github.com/McCannLab/spatial_fingerprints/workflows/R-CMD-check/badge.svg)](https://github.com/McCannLab/spatial_fingerprints/actions)

## Installation

```R
install.packages("remotes")
remotes::install_github("McCannLab/spatial_fingerprints")
```

## Data 

Raw data are included in `data/data.rda` and will be available on Dryad upon manuscript acceptance. See `R.get_data_ready.R` to see how to access the data within the package. Note that data can easy be written in `.csv` format using the following:


```R
write.csv(get_data_ready(), file = "data_f_cs.csv",
row.names = FALSE)
# transformed data (PCA)
write.csv(get_data_ready(pca = TRUE), file = "inst/extdata/data_f_pca.csv",
row.names = FALSE)
```


## Reproducibility 

Note that this project involved relatively long simulation (roughly .5core year) that were run on [Compute Canada](https://www.computecanada.ca/?lang=fr)'servers. The code to run these simulations are included in this repository: 

- for the simulations performed with [R](https://www.r-project.org/), see `R/scr_simulations.R`;
- for the simulations performed with [Julia](https://julialang.org/), see `inst/julia/scripts.jl` and `inst/julia/README.md`.

The results of those simulations were save in `output/` but not included in this repository due to size constrains (~10Go total).

Once the package is installed, the following can be run:


```R
library(spatialfingerprints)
pipeline()
```

Note that all steps to reproduce figure that require data that are not included in this repository have been commented out.