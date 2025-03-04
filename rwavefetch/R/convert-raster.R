#' Convert output of wave fetch functions to raster
#'
#' @param metric_matrix any matrix returned by one of the functions of this package
#' @param extent an optional extent object (generated with [terra::ext()])
#' @param crs an optional CRS declaration
#' @param fill_na if TRUE, it will fill all 0s with `NA`
#'
#' @return SpatRaster
#' @export
#' 
#' @details
#' Convert results to a SpatRaster.
#'
#' @examples
#' \dontrun{
#' to_raster(wave_fetch)
#' }
#'
to_raster <- function(metric_matrix, extent = NULL, crs = NULL, fill_na = TRUE) {
    if (fill_na) {
        metric_matrix[metric_matrix == 0] <- NA
    }
    conv_rast <- terra::rast(metric_matrix)
    if (!is.null(extent)) {
        if (!inherits(extent, "SpatExtent")) stop("`extent` should be of class SpatExtent")
        terra::ext(conv_rast) <- extent
    }
    if (!is.null(crs)) {
        terra::crs(conv_rast) <- crs
    }
    return(conv_rast)
}
