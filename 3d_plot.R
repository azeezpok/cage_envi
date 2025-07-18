#_____________________________________________________________________
######### Packages #######################
#_____________________________________________________________________
library("plot3D")
library(latticeExtra)

#_____________________________________________________________________
######### Data #######################
#_____________________________________________________________________

setwd("D:/")
data<- read.csv("merged_phyto.csv")
head(data)

data$variable <- as.factor(data$variable)
data$species <- as.factor(data$species)

#_____________________________________________________________________
######### 3D barplot #######################
#_____________________________________________________________________

# Define a custom color palette
custom_colors <- rainbow(length(levels(data$species)))
names(custom_colors) <- levels(data$species)

# Function to map factor levels to colors
color_map <- function(factor_levels) {
  return(custom_colors[as.character(factor_levels)])
}


#jpeg("plantkons.jpg", res = 600,height = 7,width = 12,units = "in")
#par(mar=c(5, 15, 3, 7) + 0.5)
cloud(Density ~ variable + species, data = data, panel.3d.cloud = panel.3dbars, 
      col.facet = color_map(data$species), 
      xbase = 0.4, ybase = 0.4, xlab = "",ylab = "",zlab = "",
      scales = list(arrows = FALSE, col = 1), 
      par.settings = list(axis.line = list(col = "transparent")),
      screen = list(z = 50, x = -70),
      aspect = c(3, 1, 0.5), zoom= 0.9)
dev.off()

#jpeg("planktons1_new.jpg", res = 600,height = 21,width = 32,units = "in")
# Plot with increased frame size and adjusted aspect ratio
cloud(
  Density ~ variable + species, 
  data = data, 
  panel.3d.cloud = panel.3dbars, 
  col.facet = color_map(data$species), 
  xbase = 0.4, 
  ybase = 0.4, 
  xlab = "", 
  ylab = "", 
  zlab = "",
  scales = list(
    arrows = FALSE, 
    col = 1, 
    distance = c(1, 1, 1),  # Increase scale distance
    cex = 0.8               # Adjust text size for better fit
  ), 
  par.settings = list(
    axis.line = list(col = "transparent"),
    layout.heights = list(
      top.padding = 2, 
      main.key.padding = 1, 
      key.axis.padding = 1, 
      axis.xlab.padding = 1, 
      xlab.key.padding = 1, 
      key.sub.padding = 1, 
      bottom.padding = 2
    ),
    layout.widths = list(
      left.padding = 2, 
      key.ylab.padding = 1, 
      ylab.axis.padding = 1, 
      axis.key.padding = 1, 
      right.padding = 2
    )
  ),
  screen = list(z = 50, x = -70),#(50,-70; -50,-70; 50,-20; -50,-20 combinations for different plots )
  aspect = c(3, 1, 0.5),  # Adjust aspect ratio
  zoom = 0.9
)
dev.off()
#_____________________________________________________________________
######### THE END #######################
#_____________________________________________________________________
