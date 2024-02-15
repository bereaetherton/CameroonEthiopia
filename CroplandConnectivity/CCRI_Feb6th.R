
library(raster)
library(igraph)
library(splines)
library(rworldmap)
library(sp)
library(maps)
library(geosphere)
library(RColorBrewer)
library("colorspace") 
data("countriesLow")
library(viridis)
library(tidyverse)
library(dplyr)
library(terra)
library(geohabnet)

#use geohabnet function to download potato raster

crop<-"potato"
#OR
crop<-"sweet potato"

potato<-cropharvest_rast(crop,"mapspam2017Africa") #may take a few minutes
eth_lims<-c(33,48,3,15) #crop boundaries for analysis

Ethiopia<- raster::getData("GADM",country = "Ethiopia",level = 1)

palette1 <- c( "#F4E156FF", "#F6D746FF", "#F8CD37FF", "#FAC329FF", "#FBB91EFF", "#FCAF13FF", 
               "#FCA50BFF", "#FB9C06FF", "#FA9207FF", "#F8890CFF", "#F68013FF", "#F37819FF",
               "#F06F20FF", "#EC6727FF", "#E85F2EFF", "#E25834FF", "#DD5139FF", "#D74B3FFF",
               "#D04545FF", "#CA404AFF", "#C33B4FFF", "#BC3754FF", "#B43359FF", "#AC305EFF",
               "#A42C60FF", "#9B2964FF", "#932667FF", "#922568FF", "#902568FF", "#8F2469FF",
               "#8D2369FF", "#8C2369FF", "#8A226AFF", "#88226AFF", "#87216BFF", "#85216BFF",
               "#84206BFF", "#82206CFF", "#801F6CFF", "#7F1E6CFF", "#7D1E6DFF", "#7C1D6DFF",
               "#7A1D6DFF", "#781C6DFF", "#771C6DFF", "#751B6EFF", "#741A6EFF", "#721A6EFF",
               "#71196EFF", "#6E196EFF", "#6D186EFF", "#6B186EFF", "#6A176EFF", "#68166EFF",
               "#66166EFF", "#65156EFF", "#63156EFF", "#61136EFF", "#60136EFF", "#5E126EFF",
               "#5C126EFF", "#5B126EFF", "#59106EFF", "#58106EFF", "#560F6DFF", "#540F6DFF",
               "#530E6DFF", "#510E6CFF", "#500D6CFF", "#4D0D6CFF", "#4C0C6BFF", "#4A0C6BFF",
               "#490B6AFF", "#470B6AFF", "#450A69FF", "#440A68FF", "#420A68FF", "#400A67FF",
               "#3E0966FF", "#3D0965FF", "#3B0964FF", "#390963FF", "#380962FF", "#360961FF",
               "#340A5FFF", "#320A5EFF", "#310A5CFF", "#2F0A5BFF", "#2D0B59FF", "#2B0B57FF",
               "#290B55FF", "#280B53FF", "#250C51FF", "#240C4EFF", "#230C4BFF", "#200C49FF",
               "#1F0C47FF", "#1D0C44FF", "#1C0C42FF", "#1A0C40FF", "#190C3DFF", "#170C3BFF",
               "#150B38FF", "#150B36FF", "#130A33FF", "#110A31FF", "#11092EFF", "#0F092CFF",
               "#0D082AFF", "#0C0827FF", "#0B0725FF", "#0A0723FF", "#090620FF", "#08051EFF",
               "#07051CFF", "#060419FF", "#050418FF", "#040315FF", "#040312FF", "#030210FF",
               "#02020EFF", "#02020CFF", "#02010AFF", "#010108FF", "#010106FF", "#010005FF",
               "#000004FF", "#000004FF", "#000004FF")

#----------- East hemisphere-------------------
east_ext <- extent(eth_lims)#extent(8, 18, 1, 16)
## 1.2 Customize crop and values of parameters
beta0<-0.5                                       ###
beta<-1                                          ###
beta1<-1.5                                       ###
gamma00<-0.05                                    ###
gamma0<-0.1                                      ###
gamma<-0.2                                       ###
gamma1<-0.3                                      ###
gamma2<-1                                        ###                                
cutoff1<- 0.001 #9.943357e-07                    ###
cutoff2 <- 0.01 # cutoff of adjancecy matrix     ###

cropharvestA <- potato
cropharvestA<-crop(cropharvestA,east_ext)
Resolution <- 2 # Set aggregated resolution, for example, assign 12 for 1 degree.

#----------- total mean aggregration -----------------------------
cropharvestAGGTM_crop <- aggregate(cropharvestA, fact = Resolution, fun=mean, na.rm = TRUE)
plot(cropharvestAGGTM, col = palette1)#,zlim= c(0, 2))
plot(countriesLow, add=TRUE, border = "black")
#----------- Extract cropland density data -----------------------
CropValues <- values(cropharvestAGGTM)
CropValuesAzero <- which(CropValues > cutoff) # find the cells with value > 0.0001
cropValue <- CropValues[CropValuesAzero]
#----------- Extract xy corrdination for "povalue" cells ---------
lon <- NULL # xmin
lat <- NULL # ymax

for(i in 1:length(CropValuesAzero)){
  temp <- extentFromCells(cropharvestAGGTM_crop, CropValuesAzero[i])
  AVxminO <- temp[1]
  lon <- c(lon, AVxminO)
  AVymaxO <- temp[4]
  lat <- c(lat, AVymaxO)
}

cropdata1 <- data.frame(lon, lat, cropValue)
latilongimatr <- cropdata1[ ,c(1:2)]# save the latitude and longitude as new matrix  
#---- use Geosphere package, function distVincentyEllipsoid() is used to calculate the distance, defult distance is meter
dvse <- distVincentyEllipsoid(c(0,0), cbind(1, 0)) # reference of standard distance in meter for one degree
latilongimatr <- as.matrix(latilongimatr)
TemMat <- matrix(-999, nrow( latilongimatr),nrow(latilongimatr))

for (i in 1:nrow(latilongimatr)) {
  TemMat[i, ] <- distVincentyEllipsoid(latilongimatr[i,], latilongimatr)/dvse
}
distance_matrix <- TemMat
map_grey_background <- raster("map_grey_background.tif")
#----------------------------------------------------------

map_grey_background_east <- crop(map_grey_background, east_ext)
plot(map_grey_background_east, col = "grey75",  xaxt='n',  yaxt='n', axes=F, box=F, legend = F, 
     main=paste('crop density: banana-plantain'), cex.main=0.7)
plot(cropharvestAGGTM_crop, col = palette1, xaxt='n',
     yaxt='n', axes=F, box=F, add = TRUE)
plot(countriesLow, add=TRUE, border = "darkblue")
plot(Ethiopia, col = NA, border = "darkblue", add=TRUE)

CCRI_powerlaw_function <- function(beta, cutoffadja, distance_matrix, lon, lat, cropValue, cropRaster, CellNumber)   {
  ##############################################
  distancematr <- distance_matrix # pairwise distance matrix
  distmat<-graph_from_adjacency_matrix(distancematr)
  #---- end of code
  distancematrexp <- distancematr^(-beta) #use function C=AX^(-beta), here A=1, X=distancematr
  cropmatr <- cropValue # complete gravity model with crop data
  cropmatr1 <- matrix(cropmatr, , 1 )
  cropmatr2 <- matrix(cropmatr, 1, )
  
  cropmatrix <- cropmatr1 %*% cropmatr2
  cropmatrix <- as.matrix(cropmatrix)
  cropdistancematr <- distancematrexp * cropmatrix # adjacecy matrix
  logicalmatr <- cropdistancematr > cutoffadja # adjacency matrix after threshold
  stan <- cropdistancematr * logicalmatr
  stan <- round(stan, 6) # use round() because betweenness() may have problem when do the calculation
  stan[which(stan=="Inf")]<-0 #####ADDED IN NOW
  cropdistancematrix <- graph.adjacency(stan,mode=c("undirected"),diag=F,weighted=T)#create adjacency matrix
  plot(cropdistancematrix,layout=layout_nicely(cropdistancematrix))
  edge_density(cropdistancematrix)
  ##############################################
  ## sum of nearest neighbors degree
  knnpref0<-knn(cropdistancematrix,weights=NA)$knn
  knnpref0[is.na(knnpref0)]<-0
  degreematr<-degree(cropdistancematrix)
  knnpref<-knnpref0*degreematr
  if(max(knnpref)==0){knnprefp=0}else
    if(max(knnpref)>0){knnprefp=knnpref/max(knnpref)/6}
  
  ##############################################
  #### node degree, node strengh 

  nodestrength<-graph.strength(cropdistancematrix) 
  nodestrength[is.na(nodestrength)]<-0
  if(max(nodestrength)==0){nodestr=0}else
    if(max(nodestrength)>0){nodestr=nodestrength/max(nodestrength)/6}
 
  ##############################################
  #### betweenness centrality
  
  vec<-(max(E(cropdistancematrix)$weight)*1.00001-E(cropdistancematrix)$weight)
  dt<-graph_from_adjacency_matrix(cropdistancematrix)
  between2<-betweenness(dt,weights = vec)
  between[is.na(between)]<-0
  if(max(between)==0){betweenp=0}else
    if(max(between)>0){betweenp=between/max(between)/2}
  ##############################################
  #### eigenvector and eigenvalues
  
  eigenvectorvalues<-evcent(cropdistancematrix)
  ev<-eigenvectorvalues$vector
  ev[is.na(ev)]<-0
  if(max(ev)==0){evp=0}else
    if(max(ev)!=0){evp=ev/max(ev)/6}
  
  ##############################################
  #### CCRI is a weighted mean of 4 network metric
  index<-knnprefp+evp+betweenp+nodestr
  indexpre<-cropRaster
  indexpre[]<-0
  indexpre[CellNumber]<- index
  indexv<-indexpre
  return(indexv)
}

# CCRI calculated by negative exponential function 

CCRI_negExponential_function <-function(gamma,cutoffadja, distance_matrix, lon, lat, cropValue, cropRaster, CellNumber)   {
  ##############################################
  #### create adjacency matrix
  ####
  distancematr <- distance_matrix
  eulernumber<-exp(1)
  distancematrexponential <- eulernumber ^ (-gamma * distancematr)# exponential model
  cropmatr <- cropValue # complete gravity model with crop data
  cropmatr1 <- matrix(cropmatr,,1) # complete gravity model with crop data
  cropmatr2 <- matrix(cropmatr,1,)
  cropmatrix <- cropmatr1 %*% cropmatr2
  cropmatrix <- as.matrix(cropmatrix)
  cropdistancematr <- distancematrexponential #* cropmatrix
  logicalmatr <- cropdistancematr > cutoffadja
  stan <- cropdistancematr * logicalmatr
  stan <- round(stan, 6) # use round() because betweenness() may have problem when do the calculation
  cropdistancematrix<-graph.adjacency(stan,mode=c("undirected"),diag=F,weighted=T)#create adjacency matrix
  ##############################################
  #### create network for all the selected nodes
  ####
  #V(cropdistancematrix)$color=colororder
  V(cropdistancematrix)$label.cex=0.7
  edgeweight<-E(cropdistancematrix)$weight*4000
  E(cropdistancematrix)$color="red"
  
  knnpref0<-graph.knn(cropdistancematrix,weights=NA)$knn
  knnpref0[is.na(knnpref0)]<-0
  degreematr<-degree(cropdistancematrix)
  knnpref<-knnpref0*degreematr
  if(max(knnpref)==0){knnprefp=0}else
    if(max(knnpref)>0){knnprefp=knnpref/max(knnpref)/6}
  
  ##############################################
  #### node degree, node strengh 
  ####
  nodestrength<-graph.strength(cropdistancematrix) 
  nodestrength[is.na(nodestrength)]<-0
  if(max(nodestrength)==0){nodestr=0}else
    if(max(nodestrength)>0){nodestr=nodestrength/max(nodestrength)/6}
  
  ##############################################
  #### betweenness centrality
  #### 
  
  vec<-(max(E(cropdistancematrix)$weight)*1.00001-E(cropdistancematrix)$weight)
  dt<-graph_from_adjacency_matrix(distance_matrix)
  between<-betweenness(cropdistancematrix, weights = vec)
  between[is.na(between)]<-0
  if(max(between)==0){betweenp=0}else
    if(max(between)>0){betweenp=between/max(between)/2}
  ##############################################
  #### eigenvector and eigenvalues
  #### 
  eigenvectorvalues<-evcent(cropdistancematrix)
  ev<-eigenvectorvalues$vector
  ev[is.na(ev)]<-0
  if(max(ev)==0){evp=0}else
    if(max(ev)!=0){evp=ev/max(ev)/6}
  ##############################################
  #### plot index layer
  ####    
  index<-knnprefp+evp+betweenp+nodestr
  
  indexpre<-cropRaster
  indexpre[]<-0
  indexpre[CellNumber] <- index
  indexv<-indexpre
  return(indexv)
  
}


index1 <- CCRI_powerlaw_function(beta0, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index2 <- CCRI_powerlaw_function(beta, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index3 <- CCRI_powerlaw_function(beta1, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index4 <- CCRI_negExponential_function(gamma00, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index5 <- CCRI_negExponential_function(gamma0, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index6 <- CCRI_negExponential_function(gamma, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index7 <- CCRI_negExponential_function(gamma1, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

index8 <- CCRI_negExponential_function(gamma2, cutoff1, distance_matrix, lon, lat, cropValue, cropharvestAGGTM_crop, CropValuesAzero)

mean_index_raster <- mean (index1, index2, index3, index4, index5, index6, index7, index8)
plot(mean_index_raster,col=c("white",palette1))
#mean_index_raster<- as.raster(mean_index_raster)

density<-cropharvestAGGTM_crop$potato_harv_area_all[which(cropharvestAGGTM_crop$potato_harv_area_all[]>0)]
dim(density)
ccri<-mean_index_raster$potato_harv_area_all[which(mean_index_raster$potato_harv_area_all[]>0)]
dim(ccri)
val<-mean_index_raster$potato_harv_area_all[]
pts<-which(val>0.3)
val[pts]
plot(density$potato_harv_area_all,ccri$potato_harv_area_all,
     xlab="Potato Crop Desnity",ylab="CCRI Score")
#plot(CCRI)

#---------------------------------------------------------
map_grey_background <- raster("map_grey_background.tif")
map_grey_background_east <- crop(map_grey_background, east_ext)
#indices <- which(val > 0.4)#, arr.ind = TRUE)
coords<-crds(mean_index_raster)[pts,]

plot(map_grey_background_east, col = "grey75",  xaxt='n',  yaxt='n', axes=F, box=F, legend = F, 
     main=paste('Mean in cropland connectivity risk index from sensitivity analysis:', crop), cex.main=0.7)
plot(mean_index_raster , col = c("NA",palette1),
                                      xaxt='n',
     yaxt='n', axes=F, box=F, add = TRUE,zlim=c(0,1))
points(coords,col="white")
plot(countriesLow, add=TRUE, border = "darkblue")
plot(Ethiopia, col = NA, border = "darkblue", add=TRUE)
