

#' Make a categorical raster map 
#'
#' @param rast a raster 
#' @param lbls a vector of labels for each class in the raster
#' @param ttle a map title
#' @param pal a brewer palette name
#'
#' @return a tmap object
#' @export
#'
tmapCatRast <- function(rast, lbls = NULL, ttle = NULL, pal = "Accent"){
  if(!requireNamespace("tmap", quietly = TRUE)){
    stop("please install.packages(tmap) to use this function")
  }
  if(!requireNamespace("tmaptools", quietly = TRUE)){
    stop("please install.packages(tmaptools) to use this function")
  }
  
  if(raster::is.factor(rast)){
    rast <- rast %>% 
      raster::deratify(att = 1, drop = TRUE, complete = TRUE)
  }
  
  tmap::tm_shape(rast)+
    tmap::tm_raster(style = "cat", 
              palette = getPalette(rast, pal = pal, ncolour = length(lbls)), 
              labels = lbls, title = ttle)+
    tmap::tm_legend(bg.color = "white", legend.outside = TRUE)
}

getPalette <- function(rast, pal = "Accent", ncolour = 9){
  vals <- rast %>% raster::values() %>% unique() %>% sort()
  if(min(vals) == 0){
    vals <- vals+1
  }
  pal <- tmaptools::get_brewer_pal(pal, n = ncolour, plot = FALSE)
  tibble(allVals = 1:ncolour, pal) %>% filter(allVals %in% vals) %>% 
    pull(pal)
}