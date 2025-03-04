# Load packages -----
library(rwavefetch)
library(terra)


# Settings -----
outfile_1 <- "wavefetch_fetch.tif"
outfile_2 <- "wavefetch_orientation.tif"
outfile_3 <- "wavefetch_directionality.tif"


# Load example data -----
data(land1)
data(land2)
land1 <- unwrap(land1)
land2 <- unwrap(land2)

land1
land2

plot(land1)
plot(land2)


# Prepare data ------
land1[is.na(land1)] <- 0
land1_matrix <- as.matrix(land1, wide = TRUE)
dim(land1_matrix)

land2_matrix <- as.matrix(land2, wide = TRUE)
dim(land2_matrix)

land1llx <- ext(land1)[1]
land1uly <- ext(land1)[4]
land1cellszx <- res(land1)[1]
land1cellszy <- res(land1)[2]

land2llx <- ext(land2)[1]
land2uly <- ext(land2)[4]
land2cellszx <- res(land2)[1]
land2cellszy <- res(land2)[2]

land1_ext <- ext(land1)
land1_crs <- crs(land1)
rm(land1)


# Computes the wave fetch metric ------
wave_fetch <- coastal_wave_fetch(
    land1 = land1_matrix,
    land2 = land2_matrix, 
    dx = 50, dwx = 2000, cl = 2, jit = 1,
    land1llx, land1uly, land1cellszx, land1cellszy,
    land2llx, land2uly, land2cellszx, land2cellszy
)
wave_fetch_r <- to_raster(wave_fetch, land1_ext, land1_crs)
plot(wave_fetch_r, main = "Fetch")
writeRaster(wave_fetch_r, outfile_1, overwrite = T)


# Computes the wave fetch orientation ------
wave_orientation <- coastal_wave_orientation(
    land1 = land1_matrix,
    land2 = land2_matrix, 
    dx = 50, dwx = 2000, cl = 2, jit = 1,
    land1llx, land1uly, land1cellszx, land1cellszy,
    land2llx, land2uly, land2cellszx, land2cellszy
)
wave_orientation_r <- to_raster(wave_orientation, land1_ext, land1_crs)
plot(wave_orientation_r, main = "Orientation")
writeRaster(wave_orientation_r, outfile_2, overwrite = T)


# Computes the wave fetch directionality ------
wave_direction <- coastal_wave_direction(
    land1 = land1_matrix,
    land2 = land2_matrix, 
    dx = 50, dwx = 2000, cl = 2, jit = 1,
    land1llx, land1uly, land1cellszx, land1cellszy,
    land2llx, land2uly, land2cellszx, land2cellszy
)
wave_direction_r <- to_raster(wave_direction, land1_ext, land1_crs)
plot(wave_direction_r, main = "Directionality")
writeRaster(wave_direction_r, outfile_3, overwrite = T)