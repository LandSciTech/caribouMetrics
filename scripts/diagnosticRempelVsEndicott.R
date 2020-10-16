# Diagnostic comparing our results to Rempel's
library(tmap)
library(ggplot2)
library(purrr)
devtools::load_all(".")

tmap_mode("plot")
# Report parameters #===========================================================
# To compare to version with adjusted rfu to resType run below

# # adjust rfuToResType so that OcLow and SbLow are labeled ST and PjMx1 and SbMx1
# # are labeled MIX
# rfuToResType <- mutate(rfuToResType,
#                        ResourceType = case_when(
#                          RegionalForestUnit == "SbLow" ~ "ST",
#                          RegionalForestUnit == "OcLow" ~ "ST",
#                          RegionalForestUnit == "PjMx1" ~ "MIX",
#                          RegionalForestUnit == "SbMx1" ~ "MIX",
#                          TRUE ~ ResourceType)
# )


rmdTitle <- "Compare with Rempel fixed data"
rmdDescrip <- "New version with Rempel disturbance, roads, rail, utilities. Try rastize to 16 ha. Slightly narrower offset"
rmdConclu <- "Changing to rasters instead of hexagons did not change the agreement with Rempel's results. Adding offseting improved the agreement.
\n\n We concluded after consulting Rob Rempel that the difference in the calculation given the same inputs was not significant due to the random distribtion around zero of the differences in each season. However, we still don't have a staisfactory answer for why there are differences or why these differences seem to be larger for the Fall model."

# slightly narrower offset improved slope to closer to 1 for all except TDENLF which got worse

# based on folder of .rmd file as working directory
rmdOutFile <- "../../outputs/diagnosticRVsE_wRempelDataFixed6.html"

# Load Data #===================================================================
 
# Our results
sSimData <- getSyncSimData(ssimLib = "./inputNV/Libraries/ChurchillBC",
                           ssimProject = "2008/09 Plans",
                           rfuLUPath = "inputTables/OLTBorealCaribouFURegionalForestUnitMapping.csv",
                           scnNum = 6)

# data as rasters and hexgrid provided
plcD <- raster("./inputNV/intermediaryData/plc_aligned.tif")
#plcD <- raster("inputNV/Provincial-Landcover-2000/Provincial Landcover 2000/z15-27class.tif")
eskerD <- st_read("inputNV/HornsethRempelInfo/Eskers_Ontario/Eskers_Ontario.shp")
distTypeD <- raster("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Churchill_Dist_Type.tif")
distTypeD[is.na(distTypeD)] <- 0


natDistD <- distTypeD == 1
anthroDistD <- distTypeD == 2

roadD <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Roads_Churchill.shp")
railD <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/Railyway_Churchill.shp")
utilD <- st_read("inputNV/HornsethRempelInfo/Dist_UtilLinesEtc_Churchill/UtilLines_Churchill.shp")

# linFeatDen <- raster("inputNV/HornsethRempelInfo/linFeatDen.tif")
# eskerDen <- raster("inputNV/HornsethRempelInfo/eskerDen.tif")

#linFeatD <- raster("./inputNV/intermediaryData/linFeatRastDen.tif")

projectPolyD <- st_read("./inputNV/caribouRanges/Caribou_Range_Boundary.shp", 
                        quiet = TRUE) %>% 
  filter(RANGE_NAME == "Churchill") %>% 
  st_transform(crs = st_crs(plcD))

# plcD <- plcD %>% raster::crop((raster::extent(projectPolyD) + 100000)) %>% 
#   raster::aggregate(fact = 10, fun = raster::modal)

caribouResults <- caribouHabitat(
    plc = plcD, 
    esker = eskerD, 
    fri = plcD %>% raster::setValues(NA), 
    age = plcD %>% raster::setValues(NA), 
    natDist = natDistD, 
    anthroDist = anthroDistD,
    harv = plcD %>% raster::setValues(NA),
    # linFeat = linFeatDen,
    linFeat = list(roads = roadD, rail = railD,
                   utilities = utilD),
    projectPoly = projectPolyD, 
    friLU = sSimData$friLU, 
    winArea = 5000,
    caribouRange = "Churchill", 
    padProjPoly = TRUE 
    # eskerSave = "inputNV/HornsethRempelInfo/eskerDen.tif",
    # linFeatSave = "inputNV/HornsethRempelInfo/linFeatDen.tif"
  ) 
  beepr::beep()

# Rempel's results
rempelResults592 <- st_read(paste0("./inputNV/HornsethRempelInfo", 
                                   "/Churchill_4_592_1_Oct3/",
                                   "Churchill_4_592_1.shp")) %>% 
  st_transform( crs = st_crs(caribouResults@plc))

rempelResults592 <- rempelResults592 %>% mutate(LGMD_S5 = DEC_S5+MIX_S5)

# Summarise to hexs using the rempel Hexs

caribouResult592 <- raster::extract(raster::stack(caribouResults@habitatUse, 
                                                  caribouResults@processedData),
                                    rempelResults592 %>% select(HEXID), 
                                    fun = mean, df = TRUE, sp = TRUE) %>% 
  st_as_sf() %>% select(-CONST, -other)

# test summary of anthroDist
# anthroDist592 <- raster::extract(raster::stack(caribouResults@anthroDist),
#                                     rempelResults592 %>% select(HEXID), 
#                                     fun = mean, df = TRUE, sp = TRUE) %>% 
#   st_as_sf()
# 
# natDist592 <- raster::extract(raster::stack(caribouResults@natDist),
#                                  rempelResults592 %>% select(HEXID), 
#                                  fun = mean, df = TRUE, sp = TRUE) %>% 
#   st_as_sf()

# Add category 2 
caribouResult592 <- caribouResult592 %>% 
  mutate(Cat2 = calcBinaryUse(HEXID, Spring, Summer, Fall, Winter, 
                                   "Churchill")) %>% 
  left_join(caribouResult592 %>% st_drop_geometry() %>% 
              select(PID = HEXID, Spring, Summer, Fall, Winter) %>% 
              calcBinaryUse(caribouRange = "Churchill", bySeason = TRUE),
            by = c(HEXID = "PID"))

# Visual comparison of results #================================================
mapCompare <- function(x, y, xColVar, yColVar){
  if(max(x[[xColVar]]) > 1){
    brks <- seq(from = 0, to = max(max(x[[xColVar]], na.rm = TRUE), 
                                   max(y[[yColVar]], na.rm = TRUE)),
                length.out = 10) %>% ceiling()
  } else {
    brks <- 0:10/10
  }
  tmap_arrange(tm_shape(x)+
                 tm_polygons(col = xColVar, border.col = NULL,
                             breaks = brks),
               tm_shape(y)+
                 tm_polygons(col = yColVar, border.col = NULL,
                             breaks = brks))
}

resCompMaps <- map2(c("PSP_USE", "PSU_USE", "PFA_USE", "PWI_USE"),
                    c("Spring", "Summer", "Fall", "Winter"), mapCompare, 
                    x = rempelResults592, y = caribouResult592)

# Compare each explanatory variable to locate issue
rempNames <- stringr::str_subset(names(rempelResults592), "_S5") %>% sort()

cariNames <- stringr::str_subset(names(caribouResult592),
                                 "[[:upper:]]{3}\\D+|[[:upper:]]{3}$|[[:upper:]]{2}$") %>% 
  sort() %>% stringr::str_subset("HEXID", negate = TRUE)

varCompMaps <- map2(rempNames, cariNames, mapCompare, 
                    x = rempelResults592, y = caribouResult592)

# Compare Cat 2 agreement
combinedCat2 <- caribouResult592 %>% select(HEXID, contains("Cat2")) %>% 
  full_join(rempelResults592 %>% select(HEXID, contains("CF_"), CAT2) %>% 
              st_drop_geometry(), 
            by = "HEXID") %>% 
  mutate(Cat2_agree = case_when(CAT2 == Cat2 & Cat2 == 1 ~ "Both 1",
                                CAT2 == Cat2 & Cat2 == 0 ~ "Both 0",
                                CAT2 == 1 ~ "Rempel",
                                Cat2 == 1 ~ "Endicott"))

mapCombinedCat2 <- qtm(combinedCat2, fill = "Cat2_agree", borders = NULL, fill.title = "Category 2",
    fill.palette =  tmaptools::get_brewer_pal("Paired", n = 12, plot = FALSE)[c(9:12)])

tblCombinedCat2 <- combinedCat2 %>% st_drop_geometry() %>% 
  summarise(Rempel = sum(CAT2)/n(), 
            Endicott = sum(Cat2, na.rm = TRUE)/n(),
            agreement = sum(CAT2 == Cat2)/n())

mapDiff <- function(x, xColVar, yColVar){
  xColVar <- xColVar
  yColVar <- yColVar
  xColVarE <- rlang::ensym(xColVar)
  yColVarE <- rlang::ensym(yColVar)
  
  x <- mutate(x, dif = !!xColVarE - !!yColVarE)
  
  # if(max(x$dif) > 1){
    brks <- seq(from = (min(x$dif, na.rm = TRUE) %>% round(2)) - 0.01,
                to = (max(x$dif, na.rm = TRUE) %>% round(2)) + 0.01,
                length.out = 10) 
  # } else {
  #   brks <- -5:5/10
  # }

  tm_shape(x)+
    tm_polygons(col = "dif", border.col = NULL, breaks = brks, 
                title = paste0(xColVar," - ", yColVar))+
    tm_legend(outside = TRUE)
}

comprResults <- rempelResults592 %>%
  select(HEXID, contains("_S5"), contains("USE")) %>%  
  left_join(caribouResult592 %>% st_drop_geometry(), by = "HEXID")

difMapVars <- map2(rempNames, cariNames,
                   ~mapDiff(comprResults, .x, .y))

difMapRes <- map2(c("PSP_USE", "PSU_USE", "PFA_USE", "PWI_USE"),
                   c("Spring", "Summer", "Fall", "Winter"),
                   ~mapDiff(comprResults, .x, .y))



# Quantitative comparison # ====================================================

rempelResultData <- st_drop_geometry(rempelResults592)

caribouResultData <- st_drop_geometry(caribouResult592)

compar <- full_join(rempelResultData %>% 
                      select(HEXID, contains("_S5"), matches("P.._USE")), 
                    caribouResultData, by = "HEXID", 
                    suffix = c("_remp", "endi"))

cmprRes <- compar %>% 
  mutate(pdifSP = abs((PSP_USE-Spring)/PSP_USE),
         pdifSU = abs((PSU_USE-Summer)/PSU_USE),
         pdifFA = abs((PFA_USE-Fall)/PFA_USE),
         pdifWI = abs((PWI_USE-Winter)/PWI_USE),
         difSP = (PSP_USE-Spring),
         difSU = (PSU_USE-Summer),
         difFA = (PFA_USE-Fall),
         difWI = (PWI_USE-Winter)
  )
pdf()
# Check which variables have biggest differences 
cmprRes %>% select(-contains("_"), -contains("dif"),
                   -matches("[[:upper:]][[:lower:]]", ignore.case = F), 
                   pdifWI, -HEXID) %>% 
  gather(variables, value, -pdifWI) %>% 
  ggplot(aes(value, pdifWI))+
  geom_point()+
  facet_wrap(~variables, scales = "free_x")

# Check which variables have biggest differences 
cmprRes %>% select(contains("_S5"), pdifWI) %>% 
  gather(variables, value, -pdifWI) %>% 
  ggplot(aes(value, pdifWI))+
  geom_point()+
  facet_wrap(~variables, scales = "free_x")
dev.off()

# brks <- cmprRes %>% 
#   summarise_at(vars(matches("^dif")), list(max, min)) %>% 
#   gather() %>% 
#   {c(min(.$value), max(.$value))} %>% `*`(100) %>% 
#   {c(floor(.[1]), ceiling(.[2]))}%>% `/`(100)
# 
# brks <- seq(from = brks[1], to = brks[2], length.out = 8) %>% round(3)
# 
# difResMap <- left_join(caribouResult592 %>% select(HEXID),
#                        cmprRes %>% select(HEXID, contains("dif")),
#                        by = "HEXID") %>% 
#   tm_shape()+
#   tm_fill(col = c("difSP", "difSU", "difFA", "difWI"), breaks = brks, 
#           title = "Rempel - Endicott")+
#   tm_layout(legend.outside = TRUE, panel.show = TRUE,
#             panel.labels = c("difSP", "difSU", "difFA", "difWI"))


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
         r.sq = map_dbl(linMod, ~broom::glance(.x)$r.squared) %>% round(3), 
         x = 1,
         y = 0.55)

lmRvsEPlot <- comparLong %>% 
  ggplot(aes(Rempel, Endicott))+
  geom_abline(slope = 1, intercept = 0, col = "grey40")+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_text(data = comparLM, aes(x, y, label = paste0("R sq =\n", r.sq)), 
            hjust = "right", vjust = "top")+
  facet_wrap(~season)+
  theme_bw()

graphCompare <- function(cmpr, xVar, yVar){
  linMod <- lm(as.formula(paste0(xVar, "~", yVar)), cmpr)
  mx <- max(max(cmpr[xVar]), max(cmpr[yVar]))
  plt <- cmpr %>% 
    ggplot(aes_string(xVar, yVar))+
    geom_abline(slope = 1, intercept = 0, col = "grey40")+
    geom_point()+
    geom_smooth(method = "lm")+
    annotate("text", mx, mx, label = paste0("R sq =\n", 
                                            round(summary(linMod)$r.squared,
                                            3)))+
    lims(x = c(NA, mx), y =  c(NA, mx))+
    theme_bw()
  return(plt)  
}

varCompGraphs <- map2(rempNames, cariNames, ~graphCompare(compar, .x, .y))

# Compare results of calcRSP function based on Rempel expalantory vars #========
coefTableChill <- coefTableHR %>% filter(stringr::str_detect(Range, "Churchill"))

comparCalc <- rempelResults592 %>% select(HEXID, contains("_S5")) %>% 
  rename_all(~stringr::str_remove(.x, "_S5")) %>% 
  mutate(LGMD = DEC + MIX, 
         Range = "Churchill") %>% 
  st_drop_geometry() %>% 
  calcRSP(coefTableChill) %>% 
  right_join(rempelResults592, by = "HEXID") %>% 
  select(HEXID,contains("_S5"), matches("P.._USE"),
         Spring, Summer, Fall, Winter, geometry) %>% 
  mutate(difSP = (PSP_USE-Spring),
         difSU = (PSU_USE-Summer),
         difFA = (PFA_USE-Fall),
         difWI = (PWI_USE-Winter)) %>% 
  st_as_sf()

# Map difference
brks <- comparCalc %>% st_drop_geometry() %>% 
  summarise_at(vars(matches("^dif")), list(max, min)) %>% 
  gather() %>% 
  {c(min(.$value), max(.$value))} %>% `*`(100) %>% 
  {c(floor(.[1]), ceiling(.[2]))}%>% `/`(100)

brks <- seq(from = brks[1], to = brks[2], length.out = 8) %>% round(3)

difMapCalc <- tm_shape(comparCalc)+
  tm_fill(col = c("difSP", "difSU", "difFA", "difWI"), breaks = brks)

# Graph difference by season
hexidSeasonDiffGraph <- comparCalc %>% st_drop_geometry() %>% 
  select(HEXID, matches("^dif")) %>% 
  gather(season, diff, -HEXID) %>% 
  ggplot(aes(HEXID, diff))+
  geom_point()+
  facet_wrap(~season)+
  theme_bw()

# Graph distribution of residuals in each season
difDistSeasonGraph <- comparCalc %>% st_drop_geometry() %>% 
  select(HEXID, matches("^dif")) %>% 
  gather(season, diff, -HEXID) %>% 
  ggplot(aes(diff))+
  geom_density()+
  facet_wrap(~season)+
  theme_bw()

# Graph difference across variables
difFAVarsGraph <- comparCalc %>% select(HEXID, difFA, contains("_S5")) %>% 
  st_drop_geometry() %>% 
  gather("variable", "value", -HEXID, -difFA) %>% 
  ggplot(aes(value, difFA))+
  geom_point(alpha = 0.1, shape = 16)+
  facet_wrap(~ variable, scales = "free_x")+
  theme_bw()

# linear model to check if interaction of vars explains residuals 

# Note: would need to do more in-depth modelling given non-normal distribution
# of explanatory vars

# difLm <- lm(abs(difFA) ~ (.)^2 , data = comparCalc %>%
# st_drop_geometry() %>% select(difFA, contains("_S5")))
#
# difLmCoef <- difLm %>% broom::tidy()
#
# difLmCoef %>% filter(p.value < 0.05) %>% View()
#
# eff <- ggeffects::ggeffect(difLm, terms = c("LGOP_S5", "CON_S5")) %>% plot()
#
# ggeffects::ggeffect(difLm, terms = c("CON_S5", "LGTP_S5")) %>% plot()

rmarkdown::render("scripts/RempelHabitatDefs/diagnosticRempelVsEndicott.Rmd",
                  output_file = rmdOutFile)
beepr::beep()


