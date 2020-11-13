getPalette <- function(rast, pal = "Accent", ncolour = 9){
  vals <- rast %>% raster::values() %>% unique() %>% sort()
  if(min(vals) == 0){
    vals <- vals+1
  }
  pal <- tmaptools::get_brewer_pal(pal, n = ncolour, plot = FALSE)
  tibble(allVals = 1:ncolour, pal) %>% filter(allVals %in% vals) %>% 
    pull(pal)
}

tmapCatRast <- function(rast, lbls = NULL, ttle = NULL, pal = "Accent"){
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