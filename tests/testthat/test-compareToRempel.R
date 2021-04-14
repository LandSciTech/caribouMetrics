if(interactive()){
  context("Compare results to Rempel's results")
  
  pthBase <- "../ChurchillAnalysis/"
  #pthBase <- ""
  
  plcD <- raster(paste0(pthBase, "inputNV/HornsethRempelInfo/Provincial-Landcover-2000/studyAreaMerged250.tif"))
  plcD <- reclassPLC(plcD)
  
  eskerD <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp"))
  distTypeD <- raster(paste0(pthBase, "inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Churchill_Dist_Type.tif"))
  distTypeD[is.na(distTypeD)] <- 0
  
  natDistD <- distTypeD == 1
  anthroDistD <- distTypeD == 2
  
  natDistD <- raster::resample(natDistD, plcD, method = "ngb")
  anthroDistD <- raster::resample(anthroDistD, plcD, method = "ngb")
  
  roadD <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Roads_Churchill.shp"))
  railD <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Railyway_Churchill.shp"))
  utilD <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/UtilLines_Churchill.shp"))
  
  projectPolyD <- st_read(paste0(pthBase, "inputNV/caribouRanges/Caribou_Range_Boundary.shp"), 
                          quiet = TRUE) %>% 
    filter(RANGE_NAME == "Churchill") %>% 
    st_transform(crs = st_crs(plcD))

  caribouResults <- caribouHabitat(
    landCover = plcD, 
    esker = eskerD, 
    natDist = natDistD, 
    anthroDist = anthroDistD,
    linFeat = list(roads = roadD, rail = railD,
                   utilities = utilD),
    projectPoly = projectPolyD, 
    caribouRange = "Churchill", 
    padProjPoly = TRUE 
  ) 
  # Rempel's results
  rempelResults592 <- st_read(paste0(pthBase, "inputNV/HornsethRempelInfo",
                                     "/Churchill_4_592_1_Oct3/Churchill_4_592_1.shp"), quiet = TRUE) %>%
    st_transform( crs = st_crs(caribouResults@landCover))
  
  rempelResults592 <- rempelResults592 %>% mutate(LGMD_S5 = DEC_S5+MIX_S5)
  
  # Summarise to hexs using the rempel Hexs
  caribouResult592 <- raster::extract(raster::stack(caribouResults@habitatUse,
                                                    caribouResults@processedData),
                                      rempelResults592 %>% select(HEXID),
                                      fun = mean, df = TRUE, sp = TRUE) %>%
    st_as_sf()
  
  
  compar <- full_join(rempelResults592 %>%
                        st_drop_geometry() %>%
                        select(HEXID, contains("_S5"), matches("P.._USE")),
                      caribouResult592 %>% st_drop_geometry(), by = "HEXID",
                      suffix = c("_remp", "endi")) %>%
    filter(!is.na(CON)) %>%
    select(-CONST, -other)
  
  comparLong <- compar %>%
    select(HEXID, matches("P.._USE"), Spring, Summer, Fall, Winter) %>%
    gather(season, Pred, -HEXID) %>%
    mutate(Source = ifelse(stringr::str_detect(season, "P.._USE"),
                           "Rempel", "Endicott"),
           season = case_when(season == "PFA_USE" ~ "Fall",
                              season == "PSP_USE" ~ "Spring",
                              season == "PSU_USE" ~ "Summer",
                              season == "PWI_USE" ~ "Winter",
                              TRUE ~ season)) %>%
    spread(Source, Pred)
  
  comparLM <- comparLong %>% group_by(season) %>% nest() %>%
    mutate(linMod  = map(data, ~lm(Rempel ~ Endicott, data = .x)),
           r.sq = map_dbl(linMod, ~broom::glance(.x)$r.squared),
           slope = map_dbl(linMod, ~broom::tidy(.x)$estimate[2]))
  
  test_that("slope and rsq within 0.25 of 1 and r ",{
    expect_true(max(abs(comparLM$slope-1)) < 0.06)
    expect_true(max(abs(comparLM$r.sq-1)) < 0.06)
  })
  
}


