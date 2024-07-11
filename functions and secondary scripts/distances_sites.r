###########################################
###########################################
#' Check for packages and if necessary install into library 
#+ message = FALSE
rm(list=ls())
pkgs <- c("data.table", "dplyr","terra") 

inst <- pkgs %in% installed.packages()
if (any(inst)) install.packages(pkgs[!inst])
pkg.out <- lapply(pkgs, require, character.only = TRUE)

setwd(dir="C:/Users/Duchenne/Documents/evolution_pheno_morpho/")
coords=fread("Coordinates_Beefun.csv")

coords=vect(coords,geom=c("Long", "Lat"),crs="EPSG:4326")

di=distance(coords,unit="km")
max(di)
mean(di)