#' Wrapper function to generate the population growth model parameters under
#' default conditions that match those needed by most users
#'
#' @param replicates
#' @param modelVersion
#' @param survivalModelNumber
#' @param recruitmentModelNumber
#' @param randomQuantiles
#' @param populationGrowthTable
#' 
#' @export

calcDemoParams <- function(replicates,
                           modelVersion = "Johnson",
                           survivalModelNumber = "M1",
                           recruitmentModelNumber = "M4",
                           randomQuantiles=T,
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
  
  if(randomQuantiles){
    coeffTable_S$quantiles=sample(getQuantiles(nrow(coeffTable_S$coeffTable)),replace=F)
    coeffTable_R$quantiles=sample(getQuantiles(nrow(coeffTable_R$coeffTable)),replace=F)
  }  
  
  return(list(modelVersion = modelVersion,
              coeffTable_Survival = coeffTable_S,
              coeffTable_Recruitment = coeffTable_R))  
}

getQuantiles<-function(x,low=0.025,high=0.95){
  return(low+(seq(0,x-1)/(x-1))*high)
}