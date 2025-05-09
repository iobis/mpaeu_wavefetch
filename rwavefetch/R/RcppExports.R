# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

meanC <- function(x) {
    .Call(`_rwavefetch_meanC`, x)
}

isitcoast <- function(land) {
    .Call(`_rwavefetch_isitcoast`, land)
}

isitnearcoast <- function(land, dx) {
    .Call(`_rwavefetch_isitnearcoast`, land, dx)
}

#' Compute the wave fetch metric using log-transformed distances.
#'
#' This function calculates a wave fetch metric (wx32) for coastal cells,
#' and uses a log-transformed approach for the random distance search.
#' An alternative mapping search is performed when the primary search is out-of-bound.
#'
#' The original name of this function was `coastwx32twomaplatlonf`
#'
#' @param land1 IntegerMatrix representing the primary land area.
#' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
#' @param dx Integer specifying the search distance for near-coast detection.
#' @param dwx Integer controlling the magnitude of the random search distance.
#' @param cl Integer indicating the coastal cell class to process.
#' @param nsearches Integer specifying the number of iterations for the random search.
#' @param land1llx Double: lower left x-coordinate for land1.
#' @param land1uly Double: upper left y-coordinate for land1.
#' @param land1cellszx Double: cell size in the x-direction for land1.
#' @param land1cellszy Double: cell size in the y-direction for land1.
#' @param land2llx Double: lower left x-coordinate for land2.
#' @param land2uly Double: upper left y-coordinate for land2.
#' @param land2cellszx Double: cell size in the x-direction for land2.
#' @param land2cellszy Double: cell size in the y-direction for land2.
#' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
#' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
#' @export
coastal_wave_fetch <- function(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose = 1L) {
    .Call(`_rwavefetch_coastal_wave_fetch`, land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose)
}

#' Compute the wave orientation metric using log-transformed distances.
#'
#' This function calculates a wave orientation metric (wx32) for coastal cells,
#' and uses a log-transformed approach for the random distance search.
#' An alternative mapping search is performed when the primary search is out-of-bound.
#'
#' The original name of this function was `coastwx32twomaplatlono`
#'
#' @param land1 IntegerMatrix representing the primary land area.
#' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
#' @param dx Integer specifying the search distance for near-coast detection.
#' @param dwx Integer controlling the magnitude of the random search distance.
#' @param cl Integer indicating the coastal cell class to process.
#' @param nsearches Integer specifying the number of iterations for the random search.
#' @param land1llx Double: lower left x-coordinate for land1.
#' @param land1uly Double: upper left y-coordinate for land1.
#' @param land1cellszx Double: cell size in the x-direction for land1.
#' @param land1cellszy Double: cell size in the y-direction for land1.
#' @param land2llx Double: lower left x-coordinate for land2.
#' @param land2uly Double: upper left y-coordinate for land2.
#' @param land2cellszx Double: cell size in the x-direction for land2.
#' @param land2cellszy Double: cell size in the y-direction for land2.
#' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
#' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
#' @export
coastal_wave_orientation <- function(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose = 1L) {
    .Call(`_rwavefetch_coastal_wave_orientation`, land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose)
}

#' Compute the wave directionality metric using log-transformed distances.
#'
#' This function calculates a wave directionality metric (wx32) for coastal cells,
#' and uses a log-transformed approach for the random distance search.
#' An alternative mapping search is performed when the primary search is out-of-bound.
#'
#' The original name of this function was `coastwx32twomaplatlonv`
#'
#' @param land1 IntegerMatrix representing the primary land area.
#' @param land2 IntegerMatrix representing the secondary land area used for alternative mapping.
#' @param dx Integer specifying the search distance for near-coast detection.
#' @param dwx Integer controlling the magnitude of the random search distance.
#' @param cl Integer indicating the coastal cell class to process.
#' @param nsearches Integer specifying the number of iterations for the random search.
#' @param land1llx Double: lower left x-coordinate for land1.
#' @param land1uly Double: upper left y-coordinate for land1.
#' @param land1cellszx Double: cell size in the x-direction for land1.
#' @param land1cellszy Double: cell size in the y-direction for land1.
#' @param land2llx Double: lower left x-coordinate for land2.
#' @param land2uly Double: upper left y-coordinate for land2.
#' @param land2cellszx Double: cell size in the x-direction for land2.
#' @param land2cellszy Double: cell size in the y-direction for land2.
#' @param verbose Integer print messages to track progress. Use 0 for no messages, 1 for progress bar or 2 for text messages.
#' @return A NumericMatrix with the computed wave fetch metric based on the log-transformed search approach.
#' @export
coastal_wave_direction <- function(land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose = 1L) {
    .Call(`_rwavefetch_coastal_wave_direction`, land1, land2, dx, dwx, cl, jit, land1llx, land1uly, land1cellszx, land1cellszy, land2llx, land2uly, land2cellszx, land2cellszy, verbose)
}

