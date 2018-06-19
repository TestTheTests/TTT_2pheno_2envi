######################################
#
# interpolating_adaptree_env.R
# Brett Ford
# Created 20180507
#
# This script takes the original environmental
# variables used in the Adaptree project
# and interpolates them on a 360by360 matrix
#
#
# Usage: Rscript interpolating_adaptree_env.R
#
###################################### 

#Load libraries
library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(ncf)
library(fields)
library(tidyr)

#Set paths to directories
root_path <- ("/Users/brettford/Desktop/Northeastern/slim/TTT_2pheno_2envi")
input_path <- paste0(root_path, "/inputfiles/")
### Load data


### Adaptree environments
env_AT <- read.csv(paste0(input_path, "MAT06-BLUEs&Climate-pine.csv"))
head(env_AT)

### Subset to only include environmental varibles; not phenotype variables too
env_AT <- env_AT[complete.cases(env_AT),20:ncol(env_AT)]

#Plot individuals
plot(env_AT$Longitude, env_AT$Latitude)

### Transform X and Y to match (0 to 360) while preserving
X1 <- env_AT$Longitude - min(env_AT$Longitude, na.rm=TRUE)
Y1 <- env_AT$Latitude - min(env_AT$Latitude, na.rm=TRUE)
plot(X1, Y1)

# multiply both sides by a factor to get the landscape to fit on 360 x 360
X2 <- X1*360/max(X1, na.rm=TRUE) 
Y2 <- Y1*360/max(Y1, na.rm=TRUE) 
plot(X2, Y2)

### Interpolate

#For each variable, also compare (i) distributions and (ii) autocorrelation between real environments and resampled environment from the 1R model.
env_AT$x <- X2
env_AT$y <- Y2

#create a spatial coordinates object
coordinates(env_AT) = ~x + y

grd <- expand.grid(x = 1:360, y = 1:360)  # expand points to grid
#expand.grid just creates a data base with all of the combinations of the numbers supplied to it
coordinates(grd) <- ~x + y #creates a spatialdataframe 
gridded(grd) <- TRUE #converts to a SpatialGrid-class or SpatialGridDataFrame-class
dev.off()

#Create grid
plot(grd, cex = 1.5, col = "grey", axes=T)
#plot points of lodgepole pine
points(env_AT, pch = 1, col = "red", cex = 1)

for (j in 1:22){
  # interpolate!
  thisvarname <- names(env_AT)[j]
  print(c(j, thisvarname))
  thisdat <- env_AT[,j]
  print(class(thisdat))
  print(names(env_AT[,j]))
  print(colnames(thisdat))
  print(names(thisdat))
  names(thisdat) = "e"
  head(thisdat)
  idw <- idw(formula = e ~ 1, locations = thisdat, newdata = grd)  
  
  idw.output = as.data.frame(idw)  # output is defined as a data table
  
  names(idw.output)[1:3] <- c("long", "lat", "e.pred")  # give names to the modelled variables
  
  head(idw.output)
  ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = e.pred)) + geom_point(data = env_AT2, aes(x = X2, y = Y2), shape = 21, 
                                                                                              colour = "red") #+ geom_point(data = sample_locs, aes(x = X_Pops, y = Y_Pops), shape = 25,
                                                                                                                           #colour = "yellow")
  idw.output$std_e.pred <- (idw.output$e.pred*2)/(max(idw.output$e.pred, na.rm=TRUE))-1 #put map on a scale of -1 to 1
  idw.output2 <- idw.output[,c(1:2, 5)] #make dataframe of only lat, long, and standardized e
  e_grid <- spread(idw.output2, key = long, value = std_e.pred)
  e_grid <- e_grid[order(-e_grid$lat),]
  e_grid <- e_grid[,2:ncol(e_grid)]
  write.table(e_grid, file = paste0(input_path, names(env_AT)[j], "_adaptree_map_360x360.csv"), sep=",", row.names = F, col.names = F)
}