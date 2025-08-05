#' Summarize results of Bayesian demographic model in tables
#'
#' Produces summary tables for Bayesian caribou population model
#' results. 
#'
#' @param caribouBayesDemogMod caribou Bayesian demographic model results
#'   produced by calling [caribouBayesianPM()]
#' @inheritParams caribouBayesianPM
#' @param paramTable data.frame. Optional. Scenario parameters see
#'   [simulateObservations()]
#' @param exData data.frame. Optional. Output of [simulateObservations()] that
#'   records the true population metrics of the population that observations
#'   were simulated from.
#' @param simInitial Initial simulation results, produced by calling
#'   [getSimsInitial()]
#'
#' @return a list of tables:
#' * rr.summary.all: Mean parameter values for each year and standard deviation,
#'   upper and lower credible intervals projected by the Bayesian model, as well
#'   as scenario input parameters.
#' * sim.all: Mean parameter values and upper and lower credible intervals from
#'   the initial model for each year, as well as scenario input parameters.
#' * obs.all: Observed parameter values with column "Type" identifying if it is
#'   the "true" value of the simulated population or the "observed" value
#'   simulated based on the collaring program parameters.
#' 
#' @family demography
#' @export
#'
#' @examples
#' scns <- getScenarioDefaults(projYears = 10, obsYears = 10,
#'                             obsAnthroSlope = 1, projAnthroSlope = 5,
#'                             collarCount = 20, cowMult = 5)
#'
#' simO <- simulateObservations(scns)
#'
#' out <- caribouBayesianPM(surv_data = simO$simSurvObs, recruit_data = simO$simRecruitObs,
#'                           disturbance = simO$simDisturbance,
#'                           niters=10)
#'
#' outTables <- getOutputTables(out, exData = simO$exData, paramTable = simO$paramTable,
#'                              simInitial = getSimsInitial())
#'                              
#' str(outTables, max.level = 2, give.attr = FALSE)
                             
getOutputTables <- function(caribouBayesDemogMod, 
                            startYear = min(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            endYear = max(caribouBayesDemogMod$inData$disturbanceIn$Year), 
                            paramTable = data.frame(param = "observed"),
                            exData = NULL,
                            simInitial = NULL) {
  # caribouBayesDemogMod = out; startYear = oo$minYr;endYear = oo$maxYr; simInitial = simInitial
  # exData = oo$exData; paramTable = oo$paramTable
  
  result <- caribouBayesDemogMod$result
  survInput <- caribouBayesDemogMod$result$surv_data
  recInput <- caribouBayesDemogMod$result$recruit_data
  distInput <- caribouBayesDemogMod$inData$disturbanceIn
  
  # get summary info for plots
  rr.summary <- caribouBayesDemogMod$result$summary

  survInput$Year=as.numeric(as.character(survInput$Annual))
  recInput$Year=as.numeric(as.character(recInput$Annual))
  
  obsSurv <- survInput %>% group_by(PopulationName,Year)%>% summarize(Mortalities = sum(MortalitiesCertain,na.rm=T),StartTotal = max(StartTotal,na.rm=T)) 
  obsSurv$Mean <- 1-obsSurv$Mortalities/obsSurv$StartTotal
  obsSurv$Parameter <- "Adult female survival"
  obsSurv$MetricTypeID <- "S"
  obsSurv$Type <- "observed"

  obsRec <- subset(recInput, 
                select = intersect(names(recInput), c("PopulationName","Year", "Cows", "Calves","UnknownAdults","Yearlings")))
  adult_female_proportion = 0.65; sex_ratio=0.5 #TO DO: get from model object
  obsRec$Mean <- obsRec$Calves / (obsRec$Cows + obsRec$UnknownAdults*adult_female_proportion+obsRec$Yearlings*sex_ratio)
  obsRec$Parameter <- "Recruitment"
  obsRec$MetricTypeID <- "R"
  obsRec$Type <- "observed"

  #hist(subset(rr.summary,Parameter=="Recruitment")$Mean)
  
  obsAll <- rbind(subset(obsRec, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                            "Type")),
                  subset(obsSurv, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                             "Type")))
  
  if(!is.null(exData)){
    
    exData <- merge(exData,unique(subset(rr.summary,select=c(MetricTypeID,Parameter))),all.x=T)
    names(exData)[names(exData)=="Amount"] = "Mean"
    exData$Type = "true"
    
    obsAll <- rbind(obsAll,
                    subset(exData, select = c("PopulationName","Year", "Mean", "Parameter", "MetricTypeID",
                                              "Type")))
  } 

  
  if(!is.null(distInput)){
    # combine paramTable and simDisturbance and add to all output tables, nest params in a list
    dist_params <- merge(distInput, paramTable)
  }else{
    dist_params <- paramTable
  }
  
  if(!is.null(simInitial)){
    summaries <- simInitial$summary

    if(is.element("AnthroID",names(summaries))&&any(!is.na(summaries$AnthroID))){
      
      if(!is.element("Year",names(summaries))&!is.element("Anthro",names(dist_params))){
        stop("Set disturbance in caribouBayesianPM function call in order to compare to national model simulations.", call. = FALSE)
      }
      if(!all(unique(distInput$Anthro) %in% summaries$AnthroID)){
        message("recalculating initial sims to match anthropogenic distubance scenario")
        simInitial <- getSimsInitial(cPars=paramTable)
        summaries <- simInitial$summary
      }
      #remove irrelevant disturbance combinations from the summaries, and add Year if missing.
      distMerge <- subset(dist_params, select=c(Anthro,fire_excl_anthro,Year))
      distMerge$fire_excl_anthro=round(distMerge$fire_excl_anthro);distMerge$Anthro=round(distMerge$Anthro)
      distMerge=unique(distMerge)
      names(distMerge) <- c("AnthroID","fire_excl_anthroID","Year")
      tt<- merge(summaries,distMerge)
      check <- unique(subset(tt,select=names(distMerge)))
      if(nrow(check)!=nrow(distMerge)){
        stop("Handle this case")
      }
      simBigO<-tt
      simBigO$AnthroID=NULL;simBigO$fire_excl_anthroID=NULL
    }else{
      summaries$AnthroID=NULL;summaries$fire_excl_anthroID=NULL
      by_col <- intersect(names(summaries), names(dist_params))
      if(length(by_col) == 0){
        if(any(table(summaries$Year[summaries$MetricTypeID=="recruitment"])>length(unique(summaries$PopulationName)))){
          stop("Cannot merge caribouBayesDemogMod$inData$disturbanceIn and simInitial$summary because there are no columns shared between them", call. = FALSE)
        }
        simBigO <- summaries
      }else{
        matches <- intersect(summaries[[by_col]], dist_params[[by_col]])
        if(length(matches) == 0){
          stop("Cannot merge caribouBayesDemogMod$inData$disturbanceIn and simInitial$summary because there are no overlapping values in ", by_col, call. = FALSE)
        }
        simBigO <- merge(summaries, dist_params)
      }
    }
  } else {
    simBigO <- NULL
  }
  
  rr.summary <- merge(rr.summary, dist_params)
  obsAll <- merge(obsAll, dist_params)

  return(list(rr.summary.all = rr.summary, sim.all = simBigO,
              obs.all = obsAll))
}
