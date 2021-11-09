#' caribouMetrics: Models and metrics of boreal caribou responses to forest
#' landscapes
#' 
#' caribouMetrics provides implementations of several different models.
#' Firstly, it implements the caribou resource selection probability functions 
#' described in Hornseth and Rempel (2016) Seasonal resource selection of woodland 
#' caribou (Rangifer tarandus caribou) across a gradient of anthropogenic disturbance.
#' This allows for a spatial prediction of caribou habitat use across 13 caribou ranges
#' in Ontario. The package also includes an implementation of the population and
#' demographic models in Science to inform policy: linking population dynamics 
#' to habitat for a threatened species in Canada by Johnson et. al. (2020) and 
#' the "Environment Canada Scientific Assessment to Inform the Identification of
#'     Critical Habitat for Woodland Caribou (Rangifer tarandus caribou), Boreal 
#'     Population, in Canada 2011 Update" report. These functions allow users to 
#' calculate metrics of disturbance, predict demographic rates for a given level 
#' of disturbance and simulate population growth over time.
#'
#' @docType package
#' @name caribouMetrics
#'
#' @import sf
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom raster compareRaster cover crs focal focalWeight layerize mask
#'   ncell nlayers projectRaster raster reclassify res addLayer
#' @importFrom methods as is new slot slot<- slotNames
#' @importFrom stats qbeta qlnorm qnorm rbeta rbinom rnorm
#' @importFrom utils read.csv
#' 
#'   
NULL