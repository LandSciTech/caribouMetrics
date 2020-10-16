#' lstools: A package for calculating caribou resource selection probability.
#'
#' The lstools package provides three important functions: loadStaticInputs,
#' loadSimInputs, and getCaribouHabitatUse. These can be used to calculate the
#' probability of caribou habitat use for ranges in Ontario.
#'
#' @docType package
#' @name lstools
#'
#' @import sf
#' @import dplyr
#' @import tidyr
#' @import rsyncrosim
#' @import purrr
#' @importFrom tidyselect all_of
#' @importFrom raster compareRaster cover crs focal focalWeight layerize mask ncell
#'   nlayers projectRaster raster reclassify res addLayer
#'   
#'
#'   
NULL