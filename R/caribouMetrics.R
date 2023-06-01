#' caribouMetrics: Models and metrics of boreal caribou responses to forest
#' landscapes
#'
#' The caribouMetrics R package provides reproducible and open source
#' implementations of several models of Boreal woodland caribou (_Rangifer
#' tarandus caribou_) demography and habitat use. These include a population and
#' demographic model that allows users predict to demographic rates for a given
#' level of disturbance and project population growth over time. Demographic
#' rates are predicted using model coefficients published in [Johnson et. al.
#' (2020)](doi:10.1111/1365-2664.13637).
#' Population growth is projected using a two-stage demographic model with
#' density dependence and interannual variability based on Johnson et. al.
#' (2020) but with some modifications as described in [Dyson et al.
#' (2022)](https://doi.org/10.1101/2022.06.01.494350). In addition to these
#' national scale models, we provide a simple Bayesian integrated population
#' model that integrates prior information from national analysis of
#' demographic-disturbance relationships with available local demographic data
#' to reduce uncertainty in population viability projections. Our model is an
#' extension of work by [Eacker et al. (2019)](https://doi.org/10.1002/wsb.950)
#' with some modifications and an added ability to simulate observation data
#' given parameters that define a common caribou monitoring program.   Finally,
#' caribouMetrics contains a set of functions which implement the caribou
#' resource selection probability functions (RSPF) for Ontario boreal caribou
#' ranges described in [Hornseth and Rempel
#' (2016)](https://doi.org/10.1139/cjz-2015-0101).
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
#' @importFrom stats qbeta qlnorm qnorm rbeta rbinom rnorm runif as.formula
#'   ks.test quantile time
#' @importFrom utils read.csv write.csv
#'
#' 
NULL