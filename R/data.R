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
#'   \item{Range}{the caribou range name} 
#' }
#' @source LSL script for work published in:
#'  \url{https://doi.org/10.1139/cjz-2015-0101}
"coefTableHR"

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
#' @source Rempel, R. and M. Hornseth. 2018. Range-specific seasonal resource 
#'  selection probability functions for 13 caribou ranges in Northern Ontario.
#'  Ontario Ministry of Natural Resources and Forestry, Science and Research Branch, 
#'  Peterborough, ON. Science and Research Internal File Report IFR-01. 
#'  42 p. + appends.
"threshTable"

#' Lookup table for PLC to resource type
#'
#' A table to convert Provincial Landcover Classes to resource types used in
#' caribou resource selection model.
#' 
#' @format A data frame with 27 rows and 2 variables: \describe{
#'   \item{PLCCode}{Provincial landcover number code}
#'   \item{ResourceType}{Resource type code} 
#' } 
#' @source LSL script for work published in:
#'  \url{https://doi.org/10.1139/cjz-2015-0101}
"plcToResType"

#' Lookup table for RFU to resource type
#'
#' A table to convert Regional Forest Units to resource types used in
#' caribou resource selection model
#' @format A data frame with 17 rows and 2 variables: \describe{
#'   \item{RegionalForestUnit}{Regional forest unit 5 letter or number code}
#'   \item{ResourceType}{Resource type code} 
#' } 
#' @source LSL script for work published in:
#'  \url{https://doi.org/10.1139/cjz-2015-0101}
"rfuToResType"