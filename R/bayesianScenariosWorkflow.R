#' Run the Bayesian population model for multiple parameter sets
#'
#' Define scenarios in a table and [simulateObservations()], run the
#' [bayesianTrajectoryWorkflow()] model and [compareTrajectories()] for each scenario.
#'
#' @param scns data.frame. Parameters for the simulations. See
#'   [getScenarioDefaults()] for details.
#' @param ePars list. Additional parameters passed on to
#'   [simulateObservations()]
#' @param Rep integer. Optional. If specified, select specified replicate trajectory.
#' @param returnSamples logical. Optional. If true, return full results from [bayesianTrajectoryWorkflow()].
#' @inheritParams bayesianTrajectoryWorkflow
#' @inheritParams compareTrajectories
#' @param printProgress logical. Should the scenario number and parameters be
#'   printed at each step?
#'
#' @return A list similar to [compareTrajectories()] where tables for each scenario
#'   have been appended together. Plus an error log for any scenarios that
#'   failed to run.
#'   
#' @family demography
#' @export
#'
#' @examples
#' scns <- expand.grid(
#'   obsYears =c(10, 20), collarCount = c(30, 300), cowMult = 2, collarInterval = 2,
#'   iAnthro = 0,
#'   obsAnthroSlope = 0, projAnthroSlope = 0, sQuantile = 0.9,
#'   rQuantile = 0.7, N0 = 1000
#' )
#' 
#' eParsIn <- list(collarOnTime = 4, collarOffTime = 4, collarNumYears = 3)
#' simsIn <- trajectoriesFromNational()
#' scResults <- bayesianScenariosWorkflow(scns, simsIn, eParsIn,
#'                        niters = 10)# only set to speed up example. Normally keep defaults.


bayesianScenariosWorkflow <- function(scns, simInitial,ePars=list(collarOnTime=4,collarOffTime=3,collarNumYears=4),Rep=NULL,
                      printProgress = FALSE,priors="default",
                      niters=formals(bboutools::bb_fit_survival)$niters,nthin=formals(bboutools::bb_fit_survival)$nthin,
                      returnSamples=F,...) {
  
  # ePars=eParsIn;simInitial=simBig;printProgress=F;niters = formals(bboutools::bb_fit_survival)$niters)
  scns <- getScenarioDefaults(scns)
  errorLog <- vector(mode = "list", length = nrow(scns))
  rr.summary.all <- vector(mode = "list", length = nrow(scns))
  sim.all <- vector(mode = "list", length = nrow(scns))
  obs.all <- vector(mode = "list", length = nrow(scns))
  inPriors <- priors
  for (p in 1:nrow(scns)) {
    # p=2
    if (printProgress) {
      start_time <- Sys.time()
      print(paste0("start time: ", format(start_time)))
    }
    cs <- scns[p, ]
    if (printProgress) {
      print(paste0(c(p, scns[p, ]), collapse = " "))
    }
    
    if(is.element("samples",names(simInitial))){
      if(any(c("rQuantile", "sQuantile") %in% names(cs))){
        warning("scns includes parameters for National model trajectories and ",
                "simInitial includes samples. simInitial$samples will be ",
                "supplied as trajectories in simulateObservations and rQuantile ",
                "and sQuantile will be ignored")
      }
      if(cs$sSlopeMod !=1 | cs$rSlopeMod != 1){
        warning("scns includes parameters for National model trajectories and ",
                "simInitial includes samples. simInitial$samples will be ",
                "supplied as trajectories in simulateObservations and rSlopeMod ",
                "and sSlopeMod will be ignored")
      }
      
      if(is.element("lQuantile",names(cs))&&!is.na(cs$lQuantile)){
        trajectories <- subset(simInitial$samples,LambdaPercentile == round(cs$lQuantile*100))
      }else{
        trajectories <- simInitial$samples
      }
      if(!is.null(Rep)){
        trajectories <- subset(trajectories,Replicate==Rep)
      }else{
        trajectories <- subset(trajectories,Replicate==sample(unique(trajectories$Replicate),1))
      }
      
      cs$Replicate <- unique(trajectories$Replicate)
      oo <- simulateObservations(cs, trajectories, 
                                 cowCounts = ePars$cowCounts,
                                 freqStartsByYear = ePars$freqStartsByYear,
                                 collarNumYears = ePars$collarNumYears,
                                 collarOffTime = ePars$collarOffTime,
                                 collarOnTime = ePars$collarOnTime,
                                 surv_data = simInitial$surv_data, recruit_data=simInitial$recruit_data)
      
    }else{
      oo <- simulateObservations(paramTable=cs, trajectories=NULL, 
                                 cowCounts = ePars$cowCounts,
                                 freqStartsByYear = ePars$freqStartsByYear,
                                 collarNumYears = ePars$collarNumYears,
                                 collarOffTime = ePars$collarOffTime,
                                 collarOnTime = ePars$collarOnTime,
                                 surv_data = simInitial$surv_data, recruit_data=simInitial$recruit_data)
    }
    #plot(plotSurvivalSeries(oo$simSurvObs))

    if (inPriors[[1]] == "default") {
      priors <- betaNationalPriors(cs)
    }
    out <- (bayesianTrajectoryWorkflow(
      surv_data = oo$simSurvObs, recruit_data = oo$simRecruitObs,
      disturbance = oo$simDisturbance,
      priors = priors, startYear = oo$minYr, endYear = oo$maxYr,
      N0 = cs$N0,niters=niters,nthin=nthin,returnSamples=returnSamples,...))

    if (inherits(out, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, errorLog = errorLog),
              "results/temp.Rds")
      next
    }
    
    if (inherits(out$result, "try-error")) {
      errorLog[[p]] <- list(cs = cs, error = out$result)
      saveRDS(list(rr.summary.all = rr.summary.all, sim.all = sim.all,
                   obs.all = obs.all, errorLog = errorLog),
              "results/temp.Rds")
      next
    }

    outTabs <- compareTrajectories(caribouBayesDemogMod = out, startYear = minYr,
                               endYear = maxYr, simInitial = simInitial,
                               exData = oo$exData, paramTable = oo$paramTable)
    
    rr.summary.all[[p]] <- outTabs$rr.summary.all
    sim.all[[p]] <- outTabs$sim.all
    obs.all[[p]] <- outTabs$obs.all
    
    # for good measure, can slow things down but could help with memory
    gc()
    if (printProgress) {
      end_time <- Sys.time()
      print(paste0("end time: ", format(end_time)))
      print(paste0("rep time: ", format(end_time - start_time)))
    }
  }
  if (length(purrr::compact(errorLog)) > 0) {
    print(errorLog)
  }
  ret = list(rr.summary.all = bind_rows(rr.summary.all), sim.all = bind_rows(sim.all),
             obs.all = bind_rows(obs.all), errorLog = errorLog)
  if(returnSamples){
    ret$out <- out
  }
  
  return(ret)
}
