#' Get inputs from syncrosim 
#'
#' Get data for calculations from a syncrosim library
#'
#'
setGeneric("getSyncSimData", function(ssimLib, ssimProject, scnNum, rfuLUPath) standardGeneric("getSyncSimData"))

# Fastest when library and project are supplied as Ssim objects
setMethod(
  "getSyncSimData", signature(ssimLib = "SsimLibrary", ssimProject = "Project"), 
  function(ssimLib, ssimProject, scnNum, rfuLUPath){
    # TODO: expand function to pull any other data like roads
    # syncrosim connection
    
    if(!requireNamespace("rsyncrosim", quietly = T)){
      stop("the rsyncrosim package is need to use the getSyncSimData function")
    }
    
    myResult = scenario(ssimProject, scenario = scnNum,
                        forceElements = T)
    myIC = scenario(ssimProject, scenario = scnNum, 
                    forceElements = T)
    
    # mySC IDs to FUs
    scLookup <- datasheet(myResult[[1]], "STSim_StateClass", optional = TRUE) %>% 
      dplyr::mutate(FU = as.character(StateLabelXID))
    
    rfus <-  read.csv(rfuLUPath, 
                      stringsAsFactors = FALSE)[,1:3]
    names(rfus) <- c("Stratum","FU","RFU")
    
    friLU <- left_join(scLookup, rfus %>% select(-Stratum) %>% 
                         distinct(.keep_all = TRUE), by = c("FU")) %>% 
      select(ID, RFU)
    
    mySC <-  datasheetRaster(myIC, datasheet="STSim_InitialConditionsSpatial", 
                             column="StateClassFileName")
    
    # load age and convert from classes to approx age
    myAge <-  10 * datasheetRaster(myIC, datasheet="STSim_InitialConditionsSpatial", 
                                   column="AgeFileName") - 5
    
    return(lst(mySC, myAge, friLU))
})

setMethod(
  "getSyncSimData", signature(ssimLib = "character", ssimProject = "character"), 
  function(ssimLib, ssimProject, scnNum, rfuLUPath){
  # TODO: expand function to pull any other data like roads
  # syncrosim connection
    if(!requireNamespace("rsyncrosim", quietly = T)){
      stop("the rsyncrosim package is need to use the getSyncSimData function")
    }
    
  ssimLib <-  ssimLibrary(name = ssimLib) #open library
  ssimProject <-  rsyncrosim::project(ssimLib, project = ssimProject)
  
  myResult = scenario(ssimProject, scenario = scnNum,
                      forceElements = T)
  myIC = scenario(ssimProject, scenario = scnNum, 
                  forceElements = T)
  
  # mySC IDs to FUs
  scLookup <- datasheet(myResult[[1]], "STSim_StateClass", optional = TRUE) %>% 
    dplyr::mutate(FU = as.character(StateLabelXID) %>% toupper())
  
  rfus <-  read.csv(rfuLUPath, 
                    stringsAsFactors = FALSE)[,1:3]
  names(rfus) <- c("Stratum","FU","RFU")
  
  friLU <- left_join(scLookup, rfus %>% select(-Stratum) %>% 
                          distinct(.keep_all = TRUE), by = c("FU")) %>% 
    select(ID, RFU)
  
  mySC <-  datasheetRaster(myIC, datasheet="STSim_InitialConditionsSpatial", 
                         column="StateClassFileName")
  
  # load age and convert from classes to approx age
  myAge <-  10 * datasheetRaster(myIC, datasheet="STSim_InitialConditionsSpatial", 
                               column="AgeFileName") - 5
  
  return(lst(mySC, myAge, friLU))
})

  
  
  
  
  