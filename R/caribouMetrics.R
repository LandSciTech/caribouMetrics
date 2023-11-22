#' caribouMetrics: Models and metrics of boreal caribou responses to forest
#' landscapes
#'
#' The caribouMetrics R package provides reproducible open source
#' implementations of several models of Boreal woodland caribou (_Rangifer
#' tarandus caribou_) demography and habitat use. A national two-stage demographic model with
#' density dependence and interannual variability follows [Johnson et. al.
#' (2020)](doi:10.1111/1365-2664.13637) with modifications described in [Dyson et al.
#' (2022)](https://doi.org/10.1101/2022.06.01.494350). Demographic
#' rates vary with disturbance as estimated by [Johnson et. al.
#' (2020)](doi:10.1111/1365-2664.13637). The package also includes a Bayesian population model 
#' designed to integrate prior information from Johnson et al's national analysis of demographic-disturbance 
#' relationships with available local demographic data to reduce uncertainty in population viability projections.
#' The Bayesian population model is an extension of work by 
#' [Eacker et al. (2019)](https://doi.org/10.1002/wsb.950). The national model can be used to simulate 
#' example population trajectories, and combined with a simple observation model and the Bayesian population model
#' to show how monitoring requirements depend on landscape condition. Finally,
#' caribouMetrics contains an implementation of [Hornseth and Rempel's
#' (2016)](https://doi.org/10.1139/cjz-2015-0101) Ontario boreal caribou resource selection model
#' as described in [Dyson et al. (2022)](https://doi.org/10.1101/2022.06.01.494350). Model implementation
#' is intended to be modular and flexible, allowing reuse of components in a variety of contexts including
#' projections of the cumulative effects of disturbance and climate change [(e.g. Stewart et al. 2023)](https://doi.org/10.1002/eap.2816)
#' and a [Shiny app](https://landscitech.github.io/BayesianCaribouDemographicProjection/) designed to allow allow exploration of user-specified monitoring and disturbance scenarios.
#'
#' @docType package
#' @name caribouMetrics
#'
#' @import sf
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom methods as is new slot slot<- slotNames
#' @importFrom stats qbeta qlnorm qnorm rbeta rbinom rnorm runif as.formula
#'   ks.test quantile time setNames weighted.mean var
#' @importFrom utils read.csv write.csv packageVersion
#' @importFrom terra plot
#' @importFrom rlang .data
#'
#' 
NULL