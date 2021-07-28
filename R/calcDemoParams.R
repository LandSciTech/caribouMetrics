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
                    resVar = "femaleSurvival", 
                    modVer = modelVersion, 
                    modNum = survivalModelNumber)[[1]]
  
  coeffTable_S <- buildCoefTable(DT_S, replicates)
  
  DT_R <- makePopDT(populationGrowthTable, 
                    resVar = "recruitment",
                    modVer = modelVersion, 
                    modNum = recruitmentModelNumber)[[1]]
  
  coeffTable_R <- buildCoefTable(DT_R, replicates)
  
  return(list(modelVersion = modelVersion,
              coeffTable_Survival = coeffTable_S,
              coeffTable_Recruitment = coeffTable_R))  
}