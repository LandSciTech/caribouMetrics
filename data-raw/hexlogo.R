## code to prepare `hexlogo`
library(ggplot2)
library(png)
library(grid)
library(ggforce)
library(dplyr)

temp_fn <- tempfile()

img <- readPNG("~/My Pictures/caribouMetricsLogo.png", native = TRUE)

meter <- readPNG("~/../Downloads/noun-meter-2414673.png") 

meter[which(meter < 1)] <- 0

meter <- meter[1:600, , ] 

meter[200:412,290:500,4] <- 0

meter <- meter %>% 
  rasterGrob(interpolate=TRUE, gp = gpar(fill = "white"))

boo <- readPNG("~/../Downloads/noun-caribou-930989.png")
boo <- boo[1:600, , ] %>% 
  rasterGrob(interpolate=TRUE)

g <- rasterGrob(img, interpolate=TRUE)

  # scale_y_continuous( limits = c(0,200))+
  # coord_polar(theta ='y', start= -pi/2)+
  # theme(panel.grid.major.x =  element_line(colour = "black"), 
  #       axis.ticks = element_line(colour = "black"))

p <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(meter, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  annotation_custom(boo, xmin=2.5, xmax=8.5, ymin=3.4, ymax=7.5)+
  theme_void()+
  theme(rect = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA))
p

s <- hexSticker::sticker(p, package = "caribouMetrics", 
                         p_size = 12, # This seems to behave very differently on a Mac vs PC
                         p_y = 0.5, 
                         p_color = "black", 
                         s_x = 1, s_y = 1,
                         s_width = 2.5, s_height = 1.75,
                         h_fill = "#449448",
                         h_color = "black",
                         white_around_sticker = FALSE,
                         filename = paste0(temp_fn, ".png"))

plot(s)

use_logo(paste0(temp_fn, ".png"))
