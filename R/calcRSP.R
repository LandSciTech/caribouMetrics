#'  Calculate resource selection probability
#'
#'  Calculate the resource selection probability for northern Ontario caribou in
#'  each season and range based on Hornseth and Rempel 2016.
#'
#'  This function takes a dataframe that identifies the range and the
#'  proportions of resource types and density of linear features for each
#'  polygon. A dataframe with coefficents must also be provided. The names of
#'  columns in the \code{resourceProp} and \code{coefs} dataframes must match
#'  those expected. The following table explains the meaning of the resource
#'  type codes.
#'
#'  \tabular{ll}{ 
#'  Variable Name               \tab Code  \cr 
#'  Dense deciduous             \tab DEC   \cr
#'  Mixed deciduous and conifer \tab MIX   \cr
#'  Open peatland               \tab LGOP  \cr
#'  Dense conifer               \tab CON   \cr
#'  Conifer peatland            \tab LGTP  \cr
#'  Mixed + deciduous           \tab LGMD  \cr
#'  Sparse conifer              \tab ST    \cr
#'  Open water                  \tab LGW   \cr
#'  Natural burn                \tab DTN   \cr
#'  Gravel esker                \tab ESK   \cr
#'  Linear feature              \tab TDENLF
#'  }

#' @param resourceProp A dataframe containing the proportions/densities of each
#'    resource in the polygon and identifying the range.
#' @param coefs A dataframe containing coefficients for each season resource.
#' @param seasons A character vector containing "all" (default) or any of
#'    "Spring", "Summer", "Fall", "Winter".
#' @param polygonIdCol A character vector giving the name of the column with a
#'    unique identifier for each polygon or an integer indicating its position.
#'    By default the first column is assumed to be the polygon id.
#'
#' @return Returns a dataframe with polygon ID, and probability of use for each
#'    season.
#'
#'   

setGeneric("calcRSP", function(resourceProp, ...) standardGeneric("calcRSP"))
setMethod(
  "calcRSP", signature(resourceProp = "data.frame"),
  function(resourceProp, coefs, seasons = "all", polygonIdCol = 1){
  polygonIdColNm <- names(resourceProp)[polygonIdCol]
  resourceProp <- resourceProp %>% rename(PID = all_of(polygonIdCol))
 
  if("sf" %in% class(resourceProp)){
    resourceProp <- st_drop_geometry(resourceProp)
  }
  
  # Check col names
  expectedColNamesResProp <- c("DEC", "MIX", "LGOP", "CON", "LGTP",
                               "LGMD", "ST", "LGW", "DTN", "ESK", "TDENLF")
  
  missingColsResProp <- setdiff(expectedColNamesResProp, names(resourceProp))
  
  if(length(missingColsResProp > 0)){
    stop("The column names in resourceProp do not match the expected names: ",
         paste0(expectedColNamesResProp, sep = ", "))
  }
  
  expectedColNamesCoefs <- c("Range", "Season", "Coefficient", "Variable")
  
  missingColsCoefs <- setdiff(expectedColNamesCoefs, names(coefs))
  
  if(length(missingColsCoefs > 0)){
    stop("The column names in coefs do not match the expected names: ",
         paste0(expectedColNamesCoefs, sep = ", "))
  }
  
  if(!is.data.frame(coefs)){
    stop("coefs must be a dataframe")
  }
  
  if(!is.data.frame(resourceProp)){
    stop("resourceProp must be a dataframe")
  }
  
  if(!seasons %in% c("Spring", "Summer", "Fall", "Winter", "all")){
    stop("seasons must be one of: Spring, Summer, Fall, Winter, all")
  }
  
  # Select relevant seasons
  if (seasons == "all"){
    seasons <- c("Spring", "Summer", "Fall", "Winter")
  }
  
  # add seasons for each polygon
  resourceProp <- map_dfr(seasons, ~mutate(resourceProp, Season = .x))
  
  
  # make predictions for each polygon and season
  resourceProp <- resourceProp %>% 
    mutate(CONST = 1) %>% 
    select(PID, all_of(expectedColNamesResProp), CONST, Season) %>% 
    gather("Variable", "Value", -Season, -PID) %>% 
    left_join(coefs %>% select(Season, Range, Variable, Coefficient),
              by = c("Season", "Variable")) %>% 
    
    # some resource types are not in the model for all ranges
    mutate(Coefficient = replace_na(Coefficient, 0), 
           Pred1 = Coefficient * Value) %>% 
    group_by(PID, Season) %>% 
    
    # log odds
    summarise(Probability = sum(Pred1)) %>% 
    ungroup() %>% 
    
    # probability
    mutate(Probability = 1/(1 + exp(-Probability))) %>% 
    spread(Season, Probability) %>% 
    rename(!!polygonIdColNm := PID)
  
  return(resourceProp)
})

setMethod(
  "calcRSP", signature(resourceProp = "Raster"),
  function(resourceProp, coefs, seasons = "all"){
    
    expectedColNamesResProp <- c("DEC", "MIX", "LGOP", "CON", "LGTP",
                                 "ST", "LGW", "DTN", "ESK", "TDENLF", "LGMD",
                                 "CONST")
    
    resourceProp <- resourceProp[[expectedColNamesResProp]]
    
    if(!is.data.frame(coefs)){
      stop("coefs must be a dataframe")
    }
    
    expectedColNamesCoefs <- c("Range", "Season", "Coefficient", "Variable")
    
    missingColsCoefs <- setdiff(expectedColNamesCoefs, names(coefs))
    
    if(length(missingColsCoefs > 0)){
      stop("The column names in coefs do not match the expected names: ",
           paste0(expectedColNamesCoefs, sep = ", "))
    }
    
    if(!all(seasons %in% c("Spring", "Summer", "Fall", "Winter", "all"))){
      stop("seasons must be one of: Spring, Summer, Fall, Winter, all")
    }
    
    # Select relevant seasons
    if (any(seasons == "all")){
      seasons <- c("Spring", "Summer", "Fall", "Winter")
    }
    
    # get coefficients and explanatory vars in same order
    coefs2 <- data.frame(resType = resourceProp %>% names(),
                         stringsAsFactors = FALSE) %>% 
      mutate(order = 1:n()) %>% 
      full_join(coefs, by = c(resType = "Variable")) %>% 
      spread(Season, Coefficient, fill = 0) %>% 
      select(-`<NA>`) %>% 
      gather(Season, Coefficient, -c(resType, WinArea, Range, order)) %>% 
      arrange(order) %>% 
      split(.$Season)
    
    # calculate probability 
    out <-lapply(coefs2, function(y){
      1/(1+exp(-((resourceProp*y$Coefficient) %>% sum())))
    }) %>%
      raster::stack()
    
    return(out)
  })
