## code to prepare `coefTableStd` dataset goes here

# 10000 ha scale used in all models

# Note: Models were based on standardized data, thus the magnitude of
# standardized coefficients approximately represents relative effect size. Data
# were z-deviate transformed and Bayesian logistic regression was used to
# develop models. *In the Missisa and James Bay ranges, mixed and deciduous
# classes were combined to reduce collinearity. †In spring, fall, and winter,
# the Kinloch models were applied to the Spirit range and the Nipigon summer
# model was applied to the Pagwachuan range (AUC > 0.7 for evaluated models,
# where AUC is area under the curve).

# Table 3. Top model coefficients from resource selection probability functions
# as selected through Akaike’s information criterion (see Supplementary Table
# S1)1 for six ranges and four seasons.
nms <- c( "LGW", "LGOP", "LGTP", "ST", "CON",  "DEC", "MIX", "DTN", 
          "TDENLF", "ESK")
Spring <- data.frame(
  Nipigon = c(0.590, 0.050, 0.340, 0.050, 0.308, -1.374, 0.518, -0.204, -0.175, 0.329), 
  Pagwachuan = c(0.188, 1.183, 0.205, 0.114, 0.308, -1.250, 0.021, -0.859, -0.554, 0.465),
  Kinloch = c(0.499, -0.060, 0.720, 0.377, 0.275, -0.267, -0.539, -0.250, -0.329, 0.001),
  Missisa = c(0.325, 0.407, 0.110, 0.153, 0.342, NA, -0.028, -1.102, -1.367, NA),
  James_Bay = c(-0.046, -0.780, 0.387, 0.122, -0.451, NA, -0.222, 0.057, -0.013, NA))

Summer <- data.frame(
  Nipigon = c(1.381, 0.315, 0.388, 0.652, 0.854, -0.601, 0.612, 0.152, -0.176, 0.439),
  Kinloch = c(0.685, 0.008, 0.866, 0.221, 0.569, 0.198, 0.130, -0.195, -0.254, 0.066),
  Spirit =  c(0.894, 0.800, 0.848, 0.474, 0.530, -0.472, 0.131, 0.812, -0.287, -0.175),
  Missisa = c(1.130, 1.195, 1.934, 1.080, 1.178, NA, 0.272, -1.869, 0.228, NA),
  James_Bay = c(-0.061, -0.728, 0.315, 0.052, -0.411, NA, -0.417, -0.065, -0.031, NA)
)

Fall <- data.frame(
  Nipigon = c(0.850, 0.321, 0.537, 0.348, 0.531, -1.311, 0.510, -0.059, -0.001, 0.211),
  Pagwachuan = c(0.160, 1.167, 0.838, 0.133, 1.099, 0.939, -0.890, -0.625, -1.142, 0.593),
  Kinloch = c(0.328, -0.374, 0.409, -0.014, -0.102, -0.504, -1.099, -0.625, -0.268, 0.117),
  Missisa = c(0.464, 0.656, 0.191, 0.591, 0.336, NA, -0.570, -0.301, -0.548, NA),
  James_Bay = c(-0.158, 0.103, 0.394, 0.563, -0.526, NA, -0.148, -1.253, -0.410, NA)
)

Winter <- data.frame(
  Nipigon  = c(-0.405, 0.264, 0.018, 0.188, 0.109, -0.365, -0.042, 0.559, -0.175, -0.023),
  Pagwachuan = c(0.143, 0.326, 0.201, 0.564, 0.015, -0.235, -0.718, 0.290, -1.342, 0.136),
  Kinloch = c(-0.405, 0.264, 0.018, 0.188, 0.109, -0.365, -0.042, -0.559, -0.351, -0.023),
  Missisa = c(-0.003, 1.563, 0.003, 1.492, -0.109, NA, 0.246, -1.305, 0.546, NA),
  James_Bay =  c(-0.242, 0.429, 0.040, 0.434, -0.465, NA, -0.032, 0.041, -0.300, NA)
)

coefTableStd <- bind_rows(lst(Spring, Summer, Fall, Winter ), .id = "Season") %>% 
  mutate(Variable = rep(nms, 4), 
         Spirit = ifelse(is.na(Spirit), Kinloch, Spirit),
         Pagwachuan = ifelse(is.na(Pagwachuan), Nipigon, Pagwachuan)) %>% 
  pivot_longer(c(-Season, -Variable), names_to = "Range", 
               values_to = "Coefficient") %>% 
  mutate(Variable = case_when(Range %in% c("Missisa", "James_Bay") & 
                                Variable == "MIX" ~ "LGMD", 
                              TRUE ~ Variable),
         Range = ifelse(Range == "James_Bay", "James Bay", Range), 
         WinArea = 10000) %>% 
  arrange(Range, Season) %>% 
  select(Season, Variable, Range, Coefficient, WinArea) %>% 
  filter(!is.na(Coefficient))

usethis::use_data(coefTableStd, overwrite = TRUE)
