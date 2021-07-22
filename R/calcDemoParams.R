#' Wrapper function to generate the population growth model parameters under
#' default conditions that match those needed by most users
#'
#' @param replicates
#' @param modelVersion
#' @param survivalModelNumber
#' @param recruitmentModelNumber
#' @param Type
#' @param populationGrowthTable
#' 
#' @description ...
#' 
#' @return ...
#' 
#' @export

calcDemoParams <- function(replicates,
                           modelVersion = "Johnson",
                           survivalModelNumber = "M1",
                           recruitmentModelNumber = "M4",
                           populationGrowthTable = 
                             read.csv("./data/populationGrowthTable.csv")){
  if (is.null(populationGrowthTable)) {
    stop("Please supply a directory containing a populationGrowthTable")
  }
  
  populationGrowthTable <- data.table(populationGrowthTable)
  
  DT_S <- makePopDT(populationGrowthTable, 
                             responseVariable = "femaleSurvival", 
                             modelVersion = modelVersion, 
                             modelNumber = survivalModelNumber)[[1]]
  
  coeffTable_S <- buildCoefficientsTable(DT_S, replicates)
  
  DT_R <- makePopDT(populationGrowthTable, 
                    responseVariable = "recruitment",
                    modelVersion = modelVersion, 
                    modelNumber = recruitmentModelNumber)[[1]]
  
  coeffTable_R <- buildCoefficientsTable(DT_R, replicates)
  
  return(list(modelVersion = modelVersion,
              coeffTable_Survival = coeffTable_S,
              coeffTable_Recruitment = coeffTable_R))  
}