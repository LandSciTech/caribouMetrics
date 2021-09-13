#' caribouMetrics: Models and metrics of boreal caribou responses to forest
#' landscapes
#'
#' The caribouMetrics package provides implementations of several different
#' models for boreal caribou. These include the resource selection probability
#' functions (RSPF) described in [Hornseth and Rempel
#' (2016)](https://doi.org/10.1139/cjz-2015-0101) and the disturbance metrics
#' described in Table 52 of Environment Canada Scientific Assessment to Inform
#' the Identification of Critical Habitat for Woodland Caribou (\emph{Rangifer
#' tarandus caribou}), Boreal Population, in Canada 2011 Update.
#'
#' @docType package
#' @name caribouMetrics
#'
#' @import sf
#' @import dplyr
#' @import tidyr
#' @import rsyncrosim
#' @import purrr
#' @importFrom tidyselect all_of
#' @importFrom raster compareRaster cover crs focal focalWeight layerize mask
#'   ncell nlayers projectRaster raster reclassify res addLayer
#'
#'
#'   
NULL