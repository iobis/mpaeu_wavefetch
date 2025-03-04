## code to prepare `europe_wave` dataset goes here
library(terra)
fs::dir_create("data")

land1 <- rast("data-sources/eu5dg13_latlon.tif")
land1 <- wrap(land1)
usethis::use_data(land1, overwrite = TRUE)

land2 <- rast("data-sources/eu5dg5km_latlon.tif")
land2 <- wrap(land2)
usethis::use_data(land2, overwrite = TRUE)