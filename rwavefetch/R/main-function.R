#' Calculates coastal wave metrics
#'
#' @param type one of "fetch", "direction", or "orientation"
#' @param ... additional parameters for the wave metrics functions
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
#'
#' @return A matrix with the computed wave fetch metric based on the log-transformed search approach.
#' @export
#' 
#' @details
#' This function is a wrapper around the three individual functions used to calculate
#' the metrics.
#'
#' @examples
#' \dontrun{
#' cw_result <- coastal_wave(type = "fetch", 
#'  land1 = land1_matrix,
#'  land2 = land2_matrix, 
#'  dx = 50, dwx = 2000, cl = 2, jit = 1,
#'  land1llx, land1uly, land1cellszx, land1cellszy,
#'  land2llx, land2uly, land2cellszx, land2cellszy)
#' }
#'
coastal_wave <- function(type = "fetch", ...) {
    if (length(type) > 1 || !type %in% c("fetch", "direction", "orientation")) {
        stop('Type should be one of "fetch", "direction", or "orientation"')
    }
    switch(type,
           fetch = coastal_wave_fetch(...),
           direction = coastal_wave_direction(...),
           orientation = coastal_wave_orientation(...))
}