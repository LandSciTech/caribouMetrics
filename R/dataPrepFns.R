#' Demographic data from google sheet
#' 
#' @param survey_url character. Google Sheet url.
#' @param shiny_progress logical. Is this inside a shiny app and called with `shiny::withProgress`?
#' @param i18n shiny.i18n translator, or NULL
#'
#' @examples
#'  
#' dataFromSheets("https://docs.google.com/spreadsheets/d/1i53nQrJXgrq3B6jO0ATHhSIbibtLq5TmmFL-PxGQNm8/edit?usp=sharing")
#' 
#' @export
dataFromSheets <- function(survey_url, shiny_progress = FALSE, i18n = NULL){
  if(!googlesheets4::gs4_has_token()){
    # this lets Google sheets work for a public sheet without authentication
    googlesheets4::gs4_deauth()
    
    sh_name <- tryCatch(
      googlesheets4::gs4_get(input$survey_url)$name,
      error = function(e) e
    )
    
    if(inherits(sh_name, "error") && sh_name$message == "Client error: (403) PERMISSION_DENIED"){
      message("Google sheet is private, attempting to authenticate")
      #Authenticate Google Sheets
      googlesheets4::gs4_auth(email = TRUE)
    }  
  }
  
  if(is.null(i18n)){
    i18n <- list(t = function(x, session = NULL)paste0(x))
  }

  if(shiny_progress && !rlang::is_installed("shiny")){
    warning("Package shiny is not installed. Setting shiny_progress to FALSE")
    shiny_progress <- FALSE
  }
  
  if(shiny_progress) shiny::setProgress(0.1, message = paste0(i18n$t("Downloading data from "), sh_name))
  survey_sh_names <- googlesheets4::sheet_names(survey_url)
  
  recruit_sh <- stringr::str_subset(survey_sh_names, "[R,r]ecruit_data")
  if(length(recruit_sh)<1){
    stop("The spreadsheet does not include a sheet named 'recruit_data'")
  }
  nms <- c("PopulationName", "Year", "Month", "Day", "Cows",
           "Bulls", "UnknownAdults", "Yearlings", "Calves")
  survey_recruit <- googlesheets4::read_sheet(survey_url, recruit_sh,
                                              na = "NA") %>%
    select(any_of(nms)) %>%
    filter(if_all(everything(), \(x)!is.na(x))) %>%
    bboudata::bbd_chk_data_recruitment(multi_pops = TRUE)
  
  # Error in make bbouSummary table if only 1 year
  #survey_recruit <- survey_recruit %>% group_by(PopulationName) %>%
  #  filter(n_distinct(Year) > 1)
  
  surv_sh <- stringr::str_subset(survey_sh_names, "[S,s]urv_data")
  if(length(surv_sh)<1){
    stop("The spreadsheet does not include a sheet named 'surv_data'")
  }
  survey_surv <- googlesheets4::read_sheet(survey_url, surv_sh,
                                           na = "NA") %>%
    bboudata::bbd_chk_data_survival(multi_pops = TRUE, allow_missing = TRUE)
  
  #survey_surv <- survey_surv %>% group_by(PopulationName) %>%
  #  filter(n_distinct(Year) > 1)
  
  pop_sh <- stringr::str_subset(survey_sh_names, "[P,p]opulation")
  if(length(pop_sh)<1){
    stop("The spreadsheet does not include a sheet named 'population'")
  }
  
  nms <- c("PopulationName", "Year", "FemalePopulationLower", "FemalePopulationUpper")
  
  survey_pop <- googlesheets4::read_sheet(survey_url, pop_sh,
                                          na = "NA") %>%
    select(any_of(nms))
  
  pop_nms <- purrr::map_lgl(nms,
                            \(x)stringr::str_detect(colnames(survey_pop), x) %>% any())
  
  if(!all(pop_nms)){
    stop("The population estimates sheet is missing the expected column names:",
         paste0(colnames(survey_pop)[!pop_nms], collapse = ", "))
  }
  
  N0 <- survey_pop %>% group_by(PopulationName) %>% filter(Year == max(Year))
  
  pops_run <- intersect(survey_recruit$PopulationName,
                        survey_surv$PopulationName) %>%
    intersect(N0$PopulationName)
  
  return(list(survey_surv = survey_surv, pops_run = pops_run,
              survey_recruit = survey_recruit, N0 = N0))
}

#' Functions used for data prep for disturbanceMetrics and caribouHabitat
#' 
#' @noRd
prepProjPoly <- function(projectPoly, landCover, bufferWidth, padProjPoly){
  if(st_crs(projectPoly) != st_crs(landCover)){
    projectPoly <- st_transform(projectPoly, crs = st_crs(landCover))
    message("projectPoly being transformed to have crs matching landCover")
  }
  
  projectPolyOrig <- projectPoly
  
  # union together multiple range polygons for raster processing
  projectPoly <- projectPoly %>% summarise()
  
  # pad projPoly to 3 times the window radius, using the larger if multiple
  if(padProjPoly){
    
    # window radius is radius of circle with winArea rounded to even number of
    # raster cells based on resolution
    # winRad <- (sqrt(max(winArea)*10000/pi)/res(landCover)[1]) %>% 
    #   round(digits = 0)*res(landCover)[1]
    
    winRad <- (bufferWidth/terra::res(landCover)[1]) %>% 
      round(digits = 0)*terra::res(landCover)[1] %>% 
      round()
    
    projectPoly <- projectPoly %>% st_buffer(winRad*3)
  }
  return(lst(projectPoly, projectPolyOrig))
}

prepRasts <- function(rastLst, landCover, projectPoly, tmplt = NULL, 
                      useTmplt = NULL){
  
  landCover <- checkAlign(landCover, projectPoly, "landCover", "projectPoly") 
  
  if(is.null(tmplt)){
    tmplt <- terra::rast(landCover) %>% terra::`res<-`(c(400, 400))
  }
    
  # remove NULLs from rastLst
  rastLst <- rastLst[which(!vapply(rastLst, function(x) is.null(x), 
                                   FUN.VALUE = TRUE))]
  rastCRSMatch <- lapply(rastLst, terra::compareGeom,
                         y = landCover, crs = TRUE, res = FALSE, ext = FALSE, 
                         rowcol = FALSE, stopOnError = FALSE) %>%
    unlist() %>% all()
  
  if(!rastCRSMatch){
    stop("all raster data sets must have matching CRS", call. = FALSE)
  }
  
  # check alignment of each raster against projectPoly and compareGeom with
  # landCover
  rastLst <- purrr::map2(rastLst, names(rastLst),
              ~checkAlign(.x, projectPoly, .y, "projectPoly")) 
  
  # # need to crop tmplt too so that it will match extent
  # tmplt <- cropIf(tmplt, projectPoly, "tmplt", "projectPoly")
  
  # tmplt is usually res 400 400 raster that linFeat and esker are rasterized to
  tmpltUse <- rep_len(list(), length(rastLst)) %>% as.list()
  if(!is.null(useTmplt)){
    tmpltUseInd <- which(names(rastLst) %in% useTmplt)
    for (i in tmpltUseInd) {
      tmpltUse[[i]] <- tmplt
    }
    
  }

  purrr::pwalk(list(rastLst, tmpltUse, names(rastLst)),
               ~checkCompRast(x = ..1, y = landCover, nmx = ..3,
                              nmy = "landCover", y2 = ..2))
  
  rastLst$refRast <- landCover
  
  return(rastLst)
}

loadFromFile <- function(indata){
  # remove NULLs from indata
  indata <- indata[which(!vapply(indata, function(x) is.null(x), 
                                 FUN.VALUE = TRUE))]
  
  if(length(indata) == 0){
    return(list())
  }
  
  charIn <-  indata %>% unlist(recursive = FALSE) %>% is.character()
  
  if(!charIn){
    stop("All data must be supplied as sf or raster objects or character
                 paths not a mixture of each", call. = FALSE)
  }  
  
  filesExist <- sapply(indata, file.exists)
  if(!all(filesExist)){
    stop("Path(s) for ",
         paste0(names(filesExist)[which(!filesExist)], collapse = ", "), 
         " do(es) not exist")
  }
  
  vect <- names(indata)[which(grepl(".shp$", indata))]
  rast <- names(indata)[which(!grepl(".shp$", indata))]
  
  neverVect <- c("refRast", "natDist", "anthroDist")
  neverRast <- c("projectPoly")
  
  if(any(vect %in% neverVect)){
    stop("Extension is .shp but a raster file must be provided for: ",
         paste0(vect[which(vect %in% neverVect)], collapse = ", "))
  }
  
  if(any(rast %in% neverRast)){
    stop("Extension is not .shp but a shapefile must be provided for: ",
         paste0(rast[which(rast %in% neverRast)], collapse = ", "))
  }
  
  
  indata[vect] <- lapply(indata[vect], sf::st_read, quiet = TRUE, agr = "constant")
  indata[rast]<- lapply(indata[rast], terra::rast)
  
  return(indata)
}