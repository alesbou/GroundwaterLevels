library(raster)
library(rgdal)
library(sp)
library(gstat)
library(lubridate)
library(dismo)

#Central Valley and California shapefiles
cv = readOGR(dsn = "boundaries", layer = "cv_temp")
proj <- crs(cv)
ca = readOGR(dsn = "boundaries", layer = "california")
ca <- spTransform(ca, proj)
boundbuf = readOGR(dsn = "boundaries", layer = "study_area_boundary")
boundbuf <- spTransform(boundbuf, proj)
r <- raster(boundbuf,res=5000)

#Data groundwater
stations <- read.csv("./datadownload/stations.csv", header = TRUE, sep = ",")
measurements <- read.csv("./datadownload/measurements.csv", header = TRUE, sep = ",")

#Stations in the Central Valley
stationspt <- SpatialPoints(stations[,6:5], proj4string=CRS("+init=epsg:4269"))
stationspt <- SpatialPointsDataFrame(stationspt,stations)
stationsptutm <- spTransform(stationspt,proj)
stationscv <- stationsptutm[boundbuf,]

#Merge data and measurements
data <- merge(measurements,stationscv,by="STN_ID")

#Obtain measurement 1998 - 2017
data$date <- as.Date(data$MSMT_DATE,format=("%Y-%m-%d"))
data <- subset(data, date>as.Date("1997-12-31"))
data <- subset(data, date<as.Date("2018-01-01"))
data$fall <- 0
data$fall[month(data$date)>9]<-1
data$spring <- 0
data$spring[month(data$date)>2 & month(data$date)<7] <- 1

#Data fall
datafall <- subset(data, fall==1)

#Selecting only one measurement per year
datafall$year <- format(as.Date(datafall$date, format="%d/%m/%Y"),"%Y")
datafall <- transform(datafall, cumsumyear = ave(fall, STN_ID, year, FUN=cumsum))
datafall <- subset(datafall, cumsumyear==1)

#Selection only stations with more than 5 observations
datafall <- transform(datafall, cumsumstations = ave(fall, STN_ID, FUN=cumsum))
datafallstations <- subset(datafall, cumsumstations==5)
datafallstations <- datafallstations["STN_ID"]

#Merge to obtain data with 1 observation in fall for stations with more than 5 annual observations
datastations5obs <- merge(datafall,datafallstations, by = "STN_ID", all.y=TRUE)

# RMSE function
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

idwrmse <- krigrmse <- rep(NA, 5)

listrasters = list()
for (yearn in 1998:2017){
  print(yearn)
  #Data for fall
  datay <- subset(datastations5obs, fall==1)
  datay <- subset(datay, year(datay$date)==yearn)
  datay <- subset(datay, GSE_WSE>0)
  datay <- subset(datay, GSE_WSE<750)
  #Obtain logs
  datay$log_gw <- log(datay$GSE_WSE)

  dataypt <- SpatialPoints(datay[,34:35], proj4string=proj)
  dataypt <- SpatialPointsDataFrame(dataypt,datay)
  dataypt <- dataypt[-zerodist(dataypt)[,1],] #This removes duplicate entries which make krigging not work

  nfolds <- 5
  k <- kfold(dataypt, nfolds)

  for (i in 1:nfolds) {
    train <- dataypt[k!=i,]
    test <- dataypt[k==i,]

    #Inverse distance weighted
    gs <- gstat(formula=log_gw~1,locations=train)
    p1 <- predict(gs, newdata=test, debug.level=0)$var1.pred
    idwrmse[i] <- RMSE(test$log_gw, p1)

    #Ordinary krigging
    gs <- gstat(formula=log_gw~1,locations=train)
    v<-variogram(gs,width=5000)
    fve <- fit.variogram(v,vgm(c("Exp","Mat","Sph"),range=100000))
    kr<-gstat(formula=log_gw~1,locations=train,model=fve)
    p2 <- predict(kr, newdata=test, debug.level=0)$var1.pred
    krigrmse[i] <- RMSE(test$log_gw, p2)
  }

  #Average RMSE
  avrmseidw <- mean(idwrmse)
  avrmsekrig <- mean(krigrmse)
  avgrmse <- c(avrmseidw, avrmsekrig) #create a vector of rmse
  weights <- ((1/avgrmse)/(sum(1/avgrmse)))  #create normalized weights

  # Obtain interpolated masked raster for IDW
  gs <- gstat(formula=log_gw~1,locations=dataypt)
  idw <- interpolate(r,gs)
  idw <- mask(idw,boundbuf)

  # Obtain interpolated masked raster for OK
  gs <- gstat(formula=log_gw~1,locations=dataypt)
  v<-variogram(gs,width=5000)
  fve <- fit.variogram(v,vgm(c("Exp","Mat","Sph"), range=100000))
  kr<-gstat(formula=log_gw~1,locations=dataypt,model=fve)
  g<-as(r,'SpatialGrid')
  kp<-predict(kr,g)
  ok<-brick(kp)
  krig <- mask(ok$var1.pred,boundbuf)

  #Obtain weighted raster
  ln_weighted_raster <- sum(weights[1]*idw, weights[2]*krig)
  listrasters[yearn-1997] <- exp(ln_weighted_raster)
}

#Writes raster results as rasters
for (i in 1998:2017){
  namyear <- paste("r_", i, sep = "")
  writeRaster(listrasters[[i-1997]] , namyear, format="GTiff")

}


#Calculates lineal fitting coefficients
listresultsperyear = list()
counter = 0
for (yearinit in seq(from=1997, to=2008, by=5)){
  print(yearinit)
  listresults <- matrix(ncol=2,nrow=13674) #
  for (ncell in 1:13674){
    output <- matrix(ncol=2, nrow=(2017-yearinit))
    for (yearn2 in ((yearinit+1):2017)){
      output[yearn2-yearinit,1]<-yearn2
      output[yearn2-yearinit,2]<-listrasters[[yearn2-1997]][ncell]
      output <- as.data.frame(output)
    }
    print(ncell)
    if (is.na(output[1,2])){
      listresults[ncell,1]<-NA
      listresults[ncell,2]<-NA
    }else{
      linearMod <- lm(V2 ~ V1, data=output)
      listresults[ncell,1]<-linearMod$coefficients[1]
      listresults[ncell,2]<-linearMod$coefficients[2]
    }
  }
  counter <- counter + 1
  listresultsperyear[[counter]] <- listresults
}


#Extrapolates lineal fitting into the future for each cell and writes the results as rasters
counter = 0
for (i in seq(from=20, to=10, by=-5)) {
  print(i)
  counter <- counter +1

  r2040 <- listrasters[[1]]
  r2040gp <- listrasters[[1]]
  r2020 <- listrasters[[1]]
  rslope <- listrasters[[1]]


  for (ncell in 1:13674){
    r2020[ncell]<-listresultsperyear[[counter]][ncell,1]+2020*listresultsperyear[[counter]][ncell,2]
    r2040[ncell]<-listresultsperyear[[counter]][ncell,1]+2040*listresultsperyear[[counter]][ncell,2]
    r2040gp[ncell]<-r2020[ncell] + listresultsperyear[[counter]][ncell,2]*190/20
    rslope[ncell] <- listresultsperyear[[counter]][ncell,2]
  }

  namvar <- paste("r2020_", i, sep = "")
  namvar2 <- paste("r2040_", i , sep = "")
  namvar3 <- paste("r2040gp_", i , sep = "")
  namvar4 <- paste("rslope_", i , sep = "")

  writeRaster(r2020 ,namvar, format="GTiff")
  writeRaster(r2040 ,namvar2, format="GTiff")
  writeRaster(r2040gp ,namvar3, format="GTiff")
  writeRaster(rslope ,namvar4, format="GTiff")
}