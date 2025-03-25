# High-resolution model of coastal wave fetch

Calculates wave fetch for coastal areas.

To install, use:

``` {r}
devtools::install_github("iobis/mpaeu_wavefetch/rwavefetch")
```

All wave fetch codes were developed by **Mike Burrows (SAMS)**.

## MPA Europe application

As part of MPA Europe project, a global map of wave fetch was generated. The codes to generate the global map are available on the folder `mpaeu-application`. Before running the main code (`calc-wavefetch.R`), run `get-coastline.R` to get the landmass raster.