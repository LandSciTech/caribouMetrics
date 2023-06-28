#'  Calculate resource selection probability
#'
#'  Calculate the resource selection probability for northern Ontario caribou in
#'  each season and range based on Hornseth and Rempel 2016.
#'
#'  This function takes a raster or dataframe that identifies the range and the
#'  proportions of resource types and density of linear features for each
#'  cell or polygon. A dataframe with coefficents must also be provided. The names of
#'  of the `resourceProp` and `coefs` objects must match
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
#'   resource in the polygon and identifying the range.
#' @param coefs A dataframe containing coefficients for each season resource.
#'   With in the package `coefTableHR` which is supplied with the package
#'   is used and any other table must have the same column names.
#' @param seasons A character vector containing "all" (default) or any of
#'   "Spring", "Summer", "Fall", "Winter".
#' @param doScale logical. Are the coefs standardized and therefore the
#'   predictors should be scaled?
#'
#' @return Returns a RasterStack or dataframe with probability of use for each
#'   season.
#'
#' @noRd

calcRSP <- function(resourceProp, coefs, seasons = "all", doScale = FALSE){
  
    expectedColNamesResProp <- c("DEC", "MIX", "LGOP", "CON", "LGTP",
                                 "ST", "LGW", "DTN", "ESK", "TDENLF", "LGMD",
                                 "CONST")
    
    resourceProp <- resourceProp[[expectedColNamesResProp]]
    
    if(!is.data.frame(coefs)){
      stop("coefs must be a dataframe")
    }
    
    coefs <- data.table::data.table(coefs)
    
    expectedColNamesCoefs <- c("Range", "Season", "Coefficient", "Variable", "WinArea")
    
    missingColsCoefs <- setdiff(expectedColNamesCoefs, names(coefs))
    
    if(length(missingColsCoefs) > 0){
      stop("The column names in coefs do not match the expected names: ",
           paste0(expectedColNamesCoefs, sep = ", "))
    }
    
    if(!all(seasons %in% c("Spring", "Summer", "Fall", "Winter", "all"))){
      stop("seasons must be one of: Spring, Summer, Fall, Winter, all")
    }
    
    if(doScale){
      resourceProp <- terra::scale(resourceProp)
      
      # memory safe assigning values
      resourceProp[["CONST"]] <- terra::init(
        resourceProp[["CONST"]], 
        fun = function(x){rep(1, x)}, 
        filename = tempfile(pattern = "spat", 
                            tmpdir = terra::terraOptions(print = FALSE)$tempdir, 
                            fileext = ".grd"))
    }
    
    # Select relevant seasons
    if (any(seasons == "all")){
      seasons <- c("Spring", "Summer", "Fall", "Winter")
    }
    
    # get coefficients and explanatory vars in same order
    coefs2 <- data.table::data.table(resType = resourceProp %>% names(),
                                     stringsAsFactors = FALSE) %>% 
      mutate(order = 1:n()) %>% 
      full_join(coefs, by = c(resType = "Variable"), multiple = "all") %>% 
      spread(.data$Season, .data$Coefficient, fill = 0) %>% 
      select(-"<NA>") %>% 
      gather("Season", "Coefficient", -c("resType", 'WinArea', "Range", "order")) %>% 
      arrange(.data$order)
    coefs2 <- split(coefs2, coefs2$Season)
    
    # calculate probability for each season
    out <-lapply(coefs2, function(y){
      1/(1+exp(-((resourceProp*y$Coefficient) %>% sum())))
    }) %>%
      terra::rast()
    
    return(out)
  }

