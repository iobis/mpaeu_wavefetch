# Get the landmass raster

# Load packages
library(terra)
library(fasterize)
options(timeout = 9999)
sf::sf_use_s2(FALSE)
fs::dir_create("mpaeu-data/landmass")

# Download landmass shapefile
# https://www.soest.hawaii.edu/pwessel/gshhg/
gshhg_url <- "http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip"

outfile <- file.path("mpaeu-data", basename(gshhg_url))

download.file(gshhg_url, outfile)
unzip(outfile, exdir = "mpaeu-data")
file.remove(outfile)

# Load data
coast <- vect("mpaeu-data/GSHHS_shp/f/GSHHS_f_L1.shp")
coast_ant <- vect("mpaeu-data/GSHHS_shp/f/GSHHS_f_L6.shp")
coast <- rbind(coast, coast_ant)
coast <- sf::st_as_sf(coast)

# Create high resolution rasters
deg_target <- 1/1113.25 # ~100m

target_rast <- rast(res = deg_target)

tiles <- getTileExtents(target_rast, rast(ncol = 20, nrow = 20))

rf <- lapply(seq_len(nrow(tiles)), \(idx){
    message("Processing ", idx)
    tile_rast <- raster::raster(res = deg_target)
    tile_rast <- crop(tile_rast, raster::extent(tiles[idx,]))
    coast_rast <- fasterize(coast, tile_rast)
    outfile <- paste0("mpaeu-data/landmass/landmass_100m_gshhg_", 
                      sprintf("%04d", idx), ".tif")
    terra::writeRaster(rast(coast_rast),
                       outfile,
                       datatype = "INT1U", overwrite = TRUE)
    return(outfile)
})

# Plot data
rf_vrt <- vrt(unlist(rf, use.names = F))
plot(rf_vrt)

# Get coarser resolution landmass
target_rast_c <- raster::raster(res = 0.05)
coast_rast <- fasterize(coast, target_rast_c)
coast_rast <- rast(coast_rast)
plot(coast_rast)

terra::writeRaster(coast_rast, "mpaeu-data/landmass_5km_gshhg.tif", 
                   datatype = "INT1U")
