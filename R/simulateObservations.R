#' Simulate observations
#'
#' Parameters specify a monitoring program that is applied to simulate observations from the example trajectories. 
#' Parameters for the caribou monitoring program, disturbance scenario and the true population
#' trajectory can be specified with `getScenarioDefaults()`.
#' 
#' For a detailed description of the process for simulating data see the
#' [vignette](https://landscitech.github.io/caribouMetrics/articles/bayesian-model-outputs.html#simulation-of-local-population-dynamics-and-monitoring)
#' (`vignette("bayesian-model-outputs", package = "caribouMetrics")`) and 
#' [Hughes et al. 2025](https://doi.org/10.1016/j.ecoinf.2025.103095).
#'
#' @param trajectories data.frame.
#' @param paramTable data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param printPlot logical. print a plot of the true population trajectory?
#' @param cowCounts data.frame. Optional. Number of cows counted in aerial
#'   surveys each year. If NULL, and `paramTable` contains `cowMult` the number
#'   of cows that survive calving based on the collar data is multiplied by
#'   `cowMult` to determine the number of cows counted in aerial surveys. If
#'   `paramTable` does not contain `cowMult` `paramTable$cowCount` is used to
#'   set the number of cows counted in aerial surveys each year. If a data.frame
#'   is provided it must have columns "Year"  and "Cows".
#' @param freqStartsByYear data.frame. Optional. Number of collars deployed in
#'   each year. If NULL `paramTable$collarCount` is used as the target number of
#'   collars and each year that collars are deployed they will be topped up to
#'   this number. If a data.frame is provided it must have 2 columns "Year" and
#'   "numStarts" and the "numStarts" is the absolute number of collars deployed
#'   in that year.
#' @param collarNumYears integer. Number of years until collar falls off
#' @param collarOffTime integer. Month that collars fall off. A number from 1
#'   (January) to 12 (December)
#' @param collarOnTime integer. Month that collars are deployed. A number from 1
#'   (January) to 12 (December)
#' @param caribouYearStart integer. The first month of the year for caribou.
#' @param distScen data.frame. Disturbance scenario. Must have columns "Year",
#'   "Anthro", and "fire_excl_anthro" containing the year, percentage of the
#'   landscape covered by anthropogenic disturbance buffered by 500 m, and the
#'   percentage covered by fire that does not overlap anthropogenic disturbance.
#'   See [disturbanceMetrics()]. If NULL the disturbance scenario is simulated
#'   based on `paramTable`
#' @inheritParams demographicCoefficients
#' @param writeFilesDir character. If not NULL `simSurvObs` and `simRecruitObs`
#'   results will be saved to csv files in the directory provided
#' @param surv_data data.frame. Optional existing survival data in bboudata format. Will be combined with simulated data if ... Otherwise ignored.
#' @param recruit_data data.frame. Optional existing recruitment data in bboudata format. Will be combined with simulated data if ... Otherwise ignored.
#'
#' @return a list with elements:
#'   * minYr: first year in the simulations,
#'   * maxYr: last year in the simulations,
#'   * simDisturbance: a data frame with columns Anthro, fire_excl_anthro, Total_dist, and  Year,
#'   * simSurvObs: a data frame of survival data in bboutools format,
#'   * simRecruitObs: a data frame of recruitment data in bboutools format,
#'   * exData: a tibble of expected population metrics based on the initial model,
#'   * paramTable: a data frame recording the input parameters for the simulation.
#'
#' @family demography
#' @export
#' 
#' @references    
#'   Hughes, J., Endicott, S., Calvert, A.M. and Johnson, C.A., 2025.
#'   Integration of national demographic-disturbance relationships and local
#'   data can improve caribou population viability projections and inform
#'   monitoring decisions. Ecological Informatics, 87, p.103095.
#'   <https://doi.org/10.1016/j.ecoinf.2025.103095>
#'   
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#' 
simulateObservations <- function(trajectories, paramTable,
         cowCounts = NULL,
         freqStartsByYear = NULL,
         collarNumYears = 4, collarOffTime = 4,
         collarOnTime = 4,
         caribouYearStart = 4,
         recSurveyMonth = 3,
         recSurveyDay = 15,
         distScen = NULL,
         writeFilesDir = NULL,
         surv_data=NULL,
         recruit_data=NULL) {
  #paramTable=cs;cowCounts=NULL;freqStartByYear=NULL; collarNumYears = ePars$collarNumYears
  #collarOffTime = ePars$collarOffTime; collarOnTime = ePars$collarOnTime;caribouYearStart=4; writeFileDir=NULL

  if(nrow(paramTable) > 1){
    stop("paramTable cannot have multiple rows", call. = FALSE)
  }
  
  includeTimes = seq((paramTable$startYear+paramTable$preYears),
                       (paramTable$startYear+paramTable$preYears+paramTable$obsYears))
  
  
  if(!is.null(surv_data)){
    if(!is.element("Annual",names(surv_data))){
      surv_data <- bboutools::bb_fit_survival(surv_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)$data
    }
    
    surv_data <- subset(surv_data,as.numeric(as.character(Annual))<=max(includeTimes))
    surv_data$Month <- as.numeric(as.character(surv_data$Month))
  }
  if(!is.null(recruit_data)){
    if(!is.element("Annual",names(recruit_data))){
      recruit_data <- bboutools::bb_fit_recruitment(recruit_data, multi_pop = TRUE, allow_missing = TRUE, quiet = TRUE, do_fit=FALSE)$data
    }
    
    recruit_data <- subset(recruit_data,as.numeric(as.character(Annual))<=max(includeTimes))
  }
  
  # Simulate covariate table
  if (is.null(distScen)) {
    covariates <- simCovariates(paramTable$iAnthro, paramTable$iFire, 
                                paramTable$preYears+paramTable$obsYears + paramTable$projYears, 
                                paramTable$obsAnthroSlope, paramTable$projAnthroSlope, 
                                paramTable$obsYears + paramTable$preYears + 1)
    simDisturbance <- covariates
    simDisturbance$Year <- paramTable$startYear + simDisturbance$time - 1
    
    if (!is.null(writeFilesDir)) {
      write.csv(simDisturbance,
                file.path(writeFilesDir,
                          paste0("simDisturbance", paramTable$label, ".csv")),
                row.names = FALSE)
    }
  } else {
    simDisturbance <- distScen
    simDisturbance$time <- simDisturbance$Year - paramTable$startYear + 1
    simDisturbance <- filter(simDisturbance, .data$Year <= (paramTable$startYear + paramTable$preYears + paramTable$obsYears - 1 + paramTable$projYears) &
                               .data$Year >= paramTable$startYear)
  }
  
  #remove irrelevant disturbance combinations from the example trajectory.
  distMerge <- subset(simDisturbance, select=c(Anthro,fire_excl_anthro,Year))
  distMerge$fire_excl_anthro=round(distMerge$fire_excl_anthro);distMerge$Anthro=round(distMerge$Anthro)
  distMerge=unique(distMerge)
  names(distMerge) <- c("AnthroID","fire_excl_anthroID","Year")
  tt<- merge(trajectories,distMerge)
  check <- unique(subset(tt,select=names(distMerge)))
  if(nrow(check)!=nrow(distMerge)){
    simDisturbance <- unique(subset(trajectories,select=c(AnthroID,fire_excl_anthroID,Year)))
    if(max(table(simDisturbance$Year))>1){
      stop("The example trajectories do not include the disturbance scenario specified, and they include more than one disturbance scenario. Either provide a trajectory that does not include multiple disturbance scnenario, or specify a disturbance scenario that is included in the trajectories.")
    }else{
      warning("The example trajectories do not include the disturbance scenario. Ignoring the disturbance scenario.")
    }
  }else{
    trajectories<-tt
  }
  trajectories$AnthroID=NULL;trajectories$fire_excl_anthroID=NULL
  
  #table(subset(trajectories,Replicate=="xV1")$Year)
  #subset(trajectories,Replicate=="xV1"&Year==2064&MetricTypeID=="Anthro")
  if(length(unique(table(subset(trajectories,Replicate=="xV1")$Year)))>1){stop("Deal with this case")}

  # if cowCounts and freqStartsByYear not provided build from cowCount and collarCount
  cowCountsIn <- cowCounts
  freqStartsByYearIn <- freqStartsByYear
  
  includeYears <- sort(intersect(trajectories$Year,includeTimes))
  
  iYrs <- c()
  if(!is.null(freqStartsByYear)){
    iYrs <- c(iYrs,freqStartsByYear$Year)
    if(!is.null(surv_data)){iYrs <- c(iYrs,as.numeric(as.character(surv_data$Annual)))}
  }
  if(!is.null(cowCounts)){
    iYrs <- c(iYrs,cowCounts$Year)
    if(!is.null(recruit_data)){iYrs <- c(iYrs,as.numeric(as.character(recruit_data$Annual)))}
  }
  
  if(length(iYrs)>0){
    includeYears <- sort(intersect(includeYears,iYrs))
  }

  popMetrics = subset(trajectories,is.element(Year,includeYears))
  
  includeYears = unique(popMetrics$Year)
  
  exData <- tidyr::pivot_wider(popMetrics, id_cols = c("Replicate", "Year","Timestep","PopulationName"),
                               names_from = "MetricTypeID",
                               values_from = "Amount")
  
  if(!is.null(recruit_data)){
    recruitYrs <- sort(setdiff(includeYears,subset(recruit_data,!is.na(Calves))$Annual))
  }else{
    recruitYrs <- includeYears
  }
  
  if(!is.null(cowCounts)){
    recruitYrs <-intersect(recruitYrs,cowCounts$Year)
    cowCounts <- subset(cowCounts,is.element(Year,recruitYrs))
    
    testTable(cowCounts, c("Year", "Cows"),
              req_vals = list(Year = recruitYrs))
  } else if(hasName(paramTable, "cowCount")){
    cowCounts <- expand.grid(Year = recruitYrs,
                             Cows = paramTable$cowCount)
  } else if(!hasName(paramTable, "cowCount") & !hasName(paramTable, "cowMult")){
    stop("One of cowCounts or paramTable$cowMult must be provided",
         call. = FALSE)
  }

                       
  if(!is.null(surv_data)){
    survYrs <- sort(setdiff(includeYears,subset(surv_data,!is.na(MortalitiesCertain))$Annual))
  }else{
    survYrs <- includeYears
  }

  inclSurvYrs = survYrs
  inclSurvYrs[survYrs==max(survYrs)]=survYrs[survYrs==max(survYrs)]-1
  inclSurvYrs=unique(inclSurvYrs)
  
  if(!is.null(freqStartsByYear)){
    survYrs <- intersect(survYrs, freqStartsByYear$Year) 
                         
    freqStartsByYear <- subset(freqStartsByYear,is.element(Year,survYrs))
    testTable(freqStartsByYear, c("Year","numStarts"),
              acc_vals = list(Year = survYrs))
  } else if(!is.null(paramTable$collarCount)){
    freqStartsByYear <- expand.grid(Year = survYrs,numStarts = paramTable$collarCount)
  }else {
    stop("One of freqStartsByYear or paramTable$collarCount must be provided",
         call. = FALSE)
  }
  
  if(!is.element("PopulationName",names(freqStartsByYear))){
    freqStartsByYear=merge(freqStartsByYear,data.frame(PopulationName=unique(exData$PopulationName)))
    if(!is.null(cowCounts)){
      cowCounts=merge(cowCounts,data.frame(PopulationName=unique(exData$PopulationName)))
    }
  }
  
  # simulate survival data from survival probability.
  if (is.element("collarInterval", names(paramTable)) & is.null(freqStartsByYearIn)) {
    renewYrs <- intersect(min(exData$Year) + seq(0, 100) * paramTable$collarInterval,
                          unique(exData$Year))
    freqStartsByYear$numStarts[!is.element(freqStartsByYear$Year, renewYrs)] <- 0
  }
  
  if(!is.null(surv_data)&&(length(unique(surv_data$Month))>1)){
    forceMonths = T
  }else{forceMonths=F}

  if(nrow(freqStartsByYear)>0){

    if (is.null(freqStartsByYearIn)) {
      #freqStartsByYear$numStarts=0
      simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears, collarOffTime,
                                    collarOnTime, caribouYearStart,topUp = T,forceMonths=forceMonths)
    } else {
      simSurvObs <- simSurvivalData(freqStartsByYear, exData, collarNumYears,
                                    collarOffTime, collarOnTime,caribouYearStart,forceMonths=forceMonths)
    }
    simSurvObs$survival=NULL
    simSurvObs$Annual <- simSurvObs$Year
    simSurvObs$Annual[simSurvObs$Month<caribouYearStart] <- simSurvObs$Annual[simSurvObs$Month<caribouYearStart]-1
    simSurvObs=subset(simSurvObs,is.element(Annual,survYrs))
    
    # if cowMult is provided, set cows as a function of number of surviving cows at
    # year start
    if (is.element("cowMult", names(paramTable)) & all(paramTable$cowMult != 1) & is.null(cowCountsIn)) {
      
      # if multiple months filter to the start of caribou year
      # if only one subtract all the mortalities for the year
      if(n_distinct(simSurvObs$Month) > 1){
        surveyMonth <- caribouYearStart - 1
        if(surveyMonth == 0){
          surveyMonth <- 12
        }
        survsCalving <- subset(simSurvObs, simSurvObs$Month == surveyMonth) 
      } else {
        survsCalving <- simSurvObs 
      }
      survsCalving <- survsCalving %>% 
        mutate(surviving = StartTotal - MortalitiesCertain)
      
      if (nrow(survsCalving) > 0) {
        cowCounts <- subset(survsCalving, select=c("PopulationName","Replicate","Year","StartTotal", "surviving"))
        cowCounts$Cows <- paramTable$cowMult * cowCounts$surviving
        cowCounts$surviving <- NULL
      } else {
        cowCounts <- unique(subset(trajectories,select=c("PopulationName","Replicate","Year")))
        cowCounts$Cows <- 0
      }
      cowCounts$Bulls=0;cowCounts$UnknownAdults=0;cowCounts$Yearlings=0
      cowCounts= subset(cowCounts,is.element(Year,recruitYrs))
    }
    
    # given observed total animals & proportion calfs/cows from simulation - get
    # calf/cow ratio
    simRecruitObs <- simCalfCowRatios(cowCounts, exData)
    if(nrow(simRecruitObs)>0){simRecruitObs$Day = recSurveyDay;simRecruitObs$Month = recSurveyMonth}
    if (!is.null(writeFilesDir)) {
      write.csv(simRecruitObs,
                file.path(writeFilesDir, paste0("simRecruitData", paramTable$label, ".csv")),
                row.names = FALSE)
      write.csv(simSurvObs,
                file.path(writeFilesDir, paste0("simSurvData", paramTable$label, ".csv")),
                row.names = FALSE)
    }
    if(!is.null(recruit_data)){
      if(nrow(simRecruitObs)>0){
        recruit_data <- merge(recruit_data,data.frame(Replicate=unique(simRecruitObs$Replicate)))
        missing = setdiff(names(simRecruitObs),names(recruit_data))
        add= unique(subset(simRecruitObs,select=missing))
        if(nrow(add)>1){
          stop("Error in simulateObservations: not clear how to combine simulated and existing recruitment data.")
        }
        
        if(is.element("Month",missing)){
          if(recSurveyMonth<caribouYearStart){
            recruit_data$Year = recruit_data$Year+1
          }
        }
        
        recruit_data = merge(subset(recruit_data,!is.na(Calves),select=intersect(names(recruit_data),names(simRecruitObs))),add)
        simRecruitObs = subset(simRecruitObs,!is.na(Calves))
        simRecruitObs = rbind(recruit_data,simRecruitObs)
        simRecruitObs = simRecruitObs[order(simRecruitObs$Year),]
        
        dups = table(subset(simRecruitObs,select=c(Year,Month,PopulationName,Replicate)))
        if(max(dups)>1){
          stop("Error in simulateObservations: duplication in simulated and existing recruitment data.")
        }
      }else{
        simRecruitObs <- recruit_data
        if(!is.element("Month",names(simRecruitObs))){
          if(recSurveyMonth<caribouYearStart){
            simRecruitObs$Year = simRecruitObs$Year+1
            simRecruitObs$Month = recSurveyMonth
            simRecruitObs$Day = recSurveyDay
          }
        }
      }  
    }
    
    simSurvObs = subset(simSurvObs,is.element(Annual,inclSurvYrs))
    simSurvObs$Annual <- NULL
    
    if(!is.null(surv_data)){
      if(nrow(simSurvObs>0)){
        surv_data <- merge(surv_data,data.frame(Replicate=unique(simSurvObs$Replicate)))
        surv_data$Month=as.numeric(as.character(surv_data$Month))
        missing = setdiff(names(simSurvObs),names(surv_data))
        if(length(missing)>0){
          add= unique(subset(simSurvObs,select=missing))
          if(nrow(add)>1){
            stop("Error in simulateObservations: not clear how to combine simulated and existing survival data.")
          }
          #id month/year combos that are in the existing data and remove from 
          #dups = merge(simSurvObs,subset(surv_data,!is.na(MortalitiesCertain),select=c(Year,Month)))
          
          surv_data = merge(subset(surv_data,select=intersect(names(surv_data),names(simSurvObs))),add)
        }else{
          surv_data = subset(surv_data,select=intersect(names(surv_data),names(simSurvObs)))
        }
        simSurvObs = subset(simSurvObs,!is.na(MortalitiesCertain))
        simSurvObs = rbind(subset(surv_data,!is.na(MortalitiesCertain)),simSurvObs)
        simSurvObs = simSurvObs[order(simSurvObs$Year),]
        
        dups = table(subset(simSurvObs,select=c(Year,Month,PopulationName,Replicate)))
        if(max(dups)>1){
          stop("Error in simulateObservations: duplication in simulated and existing survival data.")
        }
      }
    }else{
      surv_data <- simSurvObs 
    }
  }else{
    simSurvObs <- surv_data
    simRecruitObs <- recruit_data
    if(!is.element("Month",names(simRecruitObs))){
      if(recSurveyMonth<caribouYearStart){
        simRecruitObs$Year = simRecruitObs$Year+1
        simRecruitObs$Month = recSurveyMonth
        simRecruitObs$Day = recSurveyDay
      }
    }
  }
  if(paramTable$obsYears==0){
    simSurvObs$Mortalities[!is.na(simSurvObs$MortalitiesCertain)] = NA
    simRecruitObs$Calves[!is.na(simRecruitObs$Calves)] = NA
  }
  
  if(!is.element("UnknownAdults",names(simRecruitObs))){simRecruitObs$UnknownAdults <- 0}
  if(!is.element("Yearlings",names(simRecruitObs))){simRecruitObs$Yearlings <- 0}
  
  if(!is.element("CowsBulls",names(simRecruitObs))){
    simRecruitObs$CowsBulls <- simRecruitObs$Cows
  }
  if(!is.element("Bulls",names(simRecruitObs))){
    simRecruitObs$Bulls <- simRecruitObs$CowsBulls-simRecruitObs$Cows
  }
  
  retList = list(minYr=min(includeYears),maxYr = max(simDisturbance$Year),
                simSurvObs = simSurvObs, simRecruitObs = simRecruitObs,
                 exData = trajectories, paramTable = paramTable)
  if(nrow(simDisturbance)>0){retList$simDisturbance=simDisturbance}
  
  return(retList)
}
