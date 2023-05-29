#' Model coefficients for caribou ranges
#'
#' A table containing the coefficients modelled by Hornseth and Rempel for the
#' resource selection probability functions for 13 caribou ranges in Northern
#' Ontario.
#'
#' @format A data frame with 547 rows and 5 variables: \describe{
#'   \item{Season}{season, one of Spring, Summer, Fall, Winter}
#'   \item{Variable}{the shortform variable name for this resource type}
#'   \item{Coefficient}{modelled coefficients} 
#'   \item{WinArea}{window area for moving window average} 
#'   \item{Range}{the caribou range name}}
#'
#' @source Rempel, R.S., Carlson, M., Rodgers, A.R., Shuter, J.L., Farrell,
#'   C.E., Cairns, D., Stelfox, B., Hunt, L.M., Mackereth, R.W. and Jackson,
#'   J.M., 2021. Modeling Cumulative Effects of Climate and Development on
#'   Moose, Wolf, and Caribou Populations. The Journal of Wildlife Management.
"coefTableHR"

#' Model coefficients from standardized data for caribou ranges
#'
#' A table containing the standardized coefficients modelled by Hornseth and
#' Rempel for the resource selection probability functions for 6 caribou ranges
#' in Northern Ontario.
#'
#' @format A data frame with 224 rows and 5 variables: \describe{
#'   \item{Season}{season, one of Spring, Summer, Fall, Winter}
#'   \item{Variable}{the shortform variable name for this resource type}
#'   \item{Coefficient}{modelled coefficients} 
#'   \item{WinArea}{window area for moving window average} 
#'   \item{Range}{the caribou range name} }
#' @source Table 3 of: Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource
#'   selection of woodland caribou (Rangifer tarandus caribou) across a gradient
#'   of anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'   <https://doi.org/10.1139/cjz-2015-0101>
"coefTableStd"

#' Category 2 thresholds table
#'
#' A table containing the thresholds used to assign habitat as Category 2 (ie in
#' use) based on probability of use. Developed by Rempel and Hornseth using
#' Youden's J and a false negative cost of 5. 
#'
#' @format A data frame with 48 rows and 6 variables: \describe{
#'   \item{Range}{the caribou range name} 
#'   \item{Season}{season, one of Spring, Summer, Fall, Winter}
#'   \item{Sensitivity}{the proportion of presences correctly predicted as presence}
#'   \item{Specificity}{the proportion of absences correctly predicted as absence} 
#'   \item{Threshold}{the value for the lowest probability of use that is considered 
#'   a predicted presence}
#'   \item{AUC}{the area under the resource operating characteristic (ROC) curve; 
#'   it ranges from 0-1 and indicates the overall success of the model for predicting
#'    both presence and absence} }
#' @source Rempel, R.S., Carlson, M., Rodgers, A.R., Shuter, J.L., Farrell,
#'   C.E., Cairns, D., Stelfox, B., Hunt, L.M., Mackereth, R.W. and Jackson,
#'   J.M., 2021. Modeling Cumulative Effects of Climate and Development on
#'   Moose, Wolf, and Caribou Populations. The Journal of Wildlife Management.
"threshTable"

#' Lookup table for PLC to resource type
#'
#' A table to convert Provincial Land Cover classes to resource types used in
#' caribou resource selection model.
#' 
#' @format A data frame with 30 rows and 2 variables: \describe{
#'   \item{PLCCode}{Provincial landcover number code}
#'   \item{ResourceType}{Letter code indicating the resource type} 
#' } 
#' @source LSL script for work published in:
#'  Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
#'  woodland caribou (Rangifer tarandus caribou) across a gradient of
#'  anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'  <https://doi.org/10.1139/cjz-2015-0101>
"plcToResType"

#' Lookup table for FNLC to resource type
#'
#' A table to convert Far North Land Cover classes to resource types used in
#' caribou resource selection model.
#'
#' The FNLC classes were linked to resource types by comparing the descriptions
#' of FNLC classes to the PLC classes and then linking them to the corresponding
#' resource types
#'
#' @format A data frame with 26 rows and 2 variables: \describe{
#'   \item{PLCCode}{FNLC number code}
#'   \item{ResourceType}{Letter code indicating the resource type} }
"fnlcToResType"

#' Lookup table for RFU to resource type
#'
#' A table to convert Regional Forest Units to resource types used in
#' caribou resource selection model
#' @format A data frame with 17 rows and 2 variables: \describe{
#'   \item{RegionalForestUnit}{Regional forest unit 5 letter or number code}
#'   \item{ResourceType}{Letter code indicating the resource type} 
#' } 
#' @source LSL script for work published in:
#'  Hornseth, M.L. and Rempel, R.S., 2016. Seasonal resource selection of
#'  woodland caribou (Rangifer tarandus caribou) across a gradient of
#'  anthropogenic disturbance. Canadian Journal of Zoology, 94(2), pp.79-93.
#'  <https://doi.org/10.1139/cjz-2015-0101>
"rfuToResType"

#' Lookup table for resource type codes
#' 
#' A table to get numeric value associated with each resource type in rasters.
#' @format A data frame with 9 rows and 2 columns: \describe{
#'   \item{ResourceType}{Letter code indicating the resource type}
#'   \item{code}{the raster value corresponding to the resource type}
#'   }
"resTypeCode"
