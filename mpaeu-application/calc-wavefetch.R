# Calculate a global map of wave fetch
# Before running this code, execute `get-coastline.R`

# Load packages
library(terra)
library(sf)
library(rwavefetch)
library(furrr)

# Settings -----
outfolder <- "results"
fs::dir_create(outfolder)

# Parallel settings
n_cores <- 20 # set here the number of cores
plan(multisession, workers = n_cores)

# terra maximum memory usage by core
max_mem <- 0.8/n_cores



# Load corse layer and make grid -----
landmass_coarse <- rast("mpaeu-data/landmass_5km_gshhg.tif")

grid <- st_make_grid(landmass_coarse, cellsize = 5)
grid <- vect(grid)
grid$id <- seq_len(length(grid))
writeVector(grid, "mpaeu-data/grid_5dg.gpkg")

# Check which are completely empty or completely full (thus, no coast)
has_value <- future_map(seq_len(length(grid)), \(id){
    fine_vrt <- vrt(list.files("mpaeu-data/landmass", full.names = T))
    grid <- vect("mpaeu-data/grid_5dg.gpkg")

    sel_grid <- grid[grid$id == id,]
    sel_land <- crop(fine_vrt, sel_grid)
    vals_land <- values(sel_land)[,1]

    if (all(is.na(vals_land))) {
        return(0)
    } else if (all(!is.na(vals_land))) {
        return(1)
    } else {
        return(2)
    }
}, .progress = TRUE)
has_value <- unlist(has_value, use.names = T)

target_grid <- grid$id[has_value == 2]

# Get adjacent cells (to use with the coarse)
adjacent_grid <- adjacent(grid, type = "queen")


# Create a function to get the maps in parallel -----
process_grid <- function(id, grid_path, coarse_path,
                         fine_path, adjacent_grid, outfolder, max_mem) {
    
    terra::terraOptions(memfrac = max_mem)

    grid <- vect(grid_path)
    sel_grid <- grid[grid$id == id,]

    coarse_layer <- rast(coarse_path)

    fine_vrt <- vrt(list.files(fine_path, full.names = T))

    land1 <- crop(fine_vrt, sel_grid)

    square <- adjacent_grid[adjacent_grid[,1] == id, 2]
    square <- grid[grid$id %in% c(square, id),]

    land2 <- crop(coarse_layer, square)

    land1[is.na(land1)] <- 0
    land1_matrix <- as.matrix(land1, wide = TRUE)

    land2_matrix <- as.matrix(land2, wide = TRUE)

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
    rm(land1, land2)    

    for (metric in c("fetch", "orientation", "direction")) {
        result <- coastal_wave(
            type = metric,
            land1 = land1_matrix,
            land2 = land2_matrix, 
            dx = 50, dwx = 2000, cl = 2, jit = 1,
            land1llx, land1uly, land1cellszx, land1cellszy,
            land2llx, land2uly, land2cellszx, land2cellszy
        )
        result_r <- to_raster(result, land1_ext, land1_crs)
        
        outfile <- file.path(
            outfolder,
            paste0(
                "wavefetch_", metric, "_id",
                sprintf("%04d", id), ".tif"
            )
        )
        #plot(result_r)
        writeRaster(result_r, outfile, overwrite = T)
        rm(result)
    }

    return(paste("Grid id", id, "done"))

}



# Run in parallel ------
tictoc::tic()
result <- furrr::future_map(target_grid,
                            process_grid,
                            grid_path = "mpaeu-data/grid_5dg.gpkg",
                            coarse_path = "mpaeu-data/landmass_5km_gshhg.tif",
                            fine_path = "mpaeu-data/landmass",
                            adjacent_grid = adjacent_grid,
                            outfolder = outfolder,
                            max_mem = max_mem,
                            .progress = TRUE)
tictoc::toc() #34939.25



# Plot results ------
all_result <- list.files(outfolder, full.names = T)

fetch <- all_result[grepl("_fetch", all_result)]
fetch <- vrt(fetch)
plot(fetch)
plot(fetch, xlim = c(-4.3, -4), ylim = c(50.2, 50.5), 
     col = RColorBrewer::brewer.pal(n=11, "Spectral"))



# Optional: plot and save grid IDs. -------
library(ggplot2)
centroids <- sf::st_as_sf(centroids, coords = c("x", "y"), crs = "EPSG:4326")
ggplot() +
    geom_sf(data = sf::st_as_sf(grid)) +
    geom_text(data = centroids, aes(x = x, y=y, label = id), size = 2)
ggsave("grid_id.png", width = 20, height = 15)



# Optional: zip files -------
all_result <- list.files(outfolder, full.names = T)
fetch <- all_result[grepl("_fetch", all_result)]
direction <- all_result[grepl("_direction", all_result)]
orientation <- all_result[grepl("_orientation", all_result)]
zip::zip("wave_fetch.zip", fetch, mode = "cherry-pick")
zip::zip("wave_direction.zip", direction, mode = "cherry-pick")
zip::zip("wave_orientation.zip", orientation, mode = "cherry-pick")
