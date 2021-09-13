#' Sample demographic regression model coefficients
#'
#' A wrapper around \code{\link{getCoefs}} to select coefficients for the
#' appropriate model version and \code{\link{sampleCoefs}} to sample
#' coefficients for each replicate population, for both the survival and 
#' recruitment models. Also optionally, generates quantiles with
#' \code{\link{getQuantiles}}.
#'
#' @param replicates integer. Number of replicate populations.
#' @param modelVersion character. Which model version to use. Options are "ECCC"
#'   for the model used in the ECCC Report (2011) and "Johnson" for the model
#'   used in Johnson et. al. (2020)
#' @param survivalModelNumber,recruitmentModelNumber character. Which model
#'   number to use see \code{\link{popGrowthTableJohnsonECCC}} for options.
#' @param randomQuantiles logical. Should each replicate population be assigned a random quantile to be used for
#'   sampling the from the distribution of demographic parameters around the means?
#' @param populationGrowthTable data.frame. By default
#'   \code{\link{popGrowthTableJohnsonECCC}} is used. A custom table of model parameters
#'   can be provided but it must match the column names of
#'   \code{\link{popGrowthTableJohnsonECCC}}.
#'
#' @return A list with elements:
#' \describe{
#'     \item{"modelVersion"}{The name of the model version}
#'     \item{"coefSamples_Survival" and "coefSamples_Recruitment"}{
#'       lists with elements:
#'        \describe{
#'           \item{"coefSamples"}{Bootstrapped coefficients with \code{replicates} 
#'             rows}
#'           \item{"coefValues"}{Coefficient values taken from 
#'             \code{populationGrowthTable}}
#'           \item{"quantiles"}{A vector of randomly selected quantiles between
#'              0.025 and 0.975 with length \code{replicates}}
#'         }
#'       }
#'   } 
#'
#' @references ECCC. 2011. Scientific assessment to inform the identification of
#'   critical habitat for woodland caribou (Rangifer tarandus caribou), boreal
#'   population, in Canada. Canadian Wildlife Service, Ottawa.
#'   \url{http://epe.lac-bac.gc.ca/100/200/301/environment_can/2011/scientific_assessment_inform-ef/CW66-296-2011-eng.pdf}.
#'    Accessed 26 Mar 2021.
#'
#'   Johnson, C.A., Sutherland, G.D., Neave, E., Leblond, M., Kirby, P.,
#'   Superbie, C. and McLoughlin, P.D., 2020. Science to inform policy: linking
#'   population dynamics to habitat for a threatened species in Canada. Journal
#'   of Applied Ecology, 57(7), pp.1314-1327.
#'   \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.13637}
#'
#' @export

demographicCoefficients <- function(replicates,
                           modelVersion = "Johnson",
                           survivalModelNumber = "M1",
                           recruitmentModelNumber = "M4",
                           randomQuantiles = TRUE,
                           populationGrowthTable = popGrowthTableJohnsonECCC){
  
  if(!all(colnames(popGrowthTableJohnsonECCC) %in% 
          colnames(populationGrowthTable))) {
    stop("populationGrowthTable must contain all colnames in popGrowthTableJohnsonECCC")
  }
  
  populationGrowthTable <- data.table::data.table(populationGrowthTable)
  
  DT_S <- getCoefs(populationGrowthTable, 
                    resVar = "femaleSurvival", 
                    modVer = modelVersion, 
                    modNum = survivalModelNumber)[[1]]
  
  coefSamples_S <- sampleCoefs(DT_S, replicates)
  
  DT_R <- getCoefs(populationGrowthTable, 
                    resVar = "recruitment",
                    modVer = modelVersion, 
                    modNum = recruitmentModelNumber)[[1]]
  
  coefSamples_R <- sampleCoefs(DT_R, replicates)
  
  if(randomQuantiles){
    coefSamples_S$quantiles=sample(getQuantiles(nrow(coefSamples_S$coefSamples)),replace=F)
    coefSamples_R$quantiles=sample(getQuantiles(nrow(coefSamples_R$coefSamples)),replace=F)
  }  
  
  return(list(modelVersion = modelVersion,
              coefSamples_Survival = coefSamples_S,
              coefSamples_Recruitment = coefSamples_R))  
}

getQuantiles<-function(x,low=0.025,high=0.975){
  return(low+(seq(0,x-1)/(x-1))*(high-low))
}