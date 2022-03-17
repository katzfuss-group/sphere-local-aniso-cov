################################################################################################
## Mollweide projection plot
################################################################################################
library(oce)

# full data
# z.all
zlim=c(min(z.all),max(z.all))
z.all.mat <- matrix(z.all,nrow = length(lon))

cm <- colormap(zlim = zlim,name="gmt_gebco")
drawPalette(colormap=cm)
mapPlot(longitude=c(-180,180,0,0), latitude=c(0,0,-90,90), projection="+proj=moll", 
        grid=TRUE, col="lightgray",drawBox=FALSE,
        longitudelim=c(-180,180),type="n") # defaults to moll projection
mapImage(longitude=lon*180/pi, latitude=lat*180/pi, z=z.all.mat,zlim=zlim, colormap=cm,missingColor=NA)


        
# training data
z.train = z.all
z.train[mask.test]=NA
z.train.mat <- matrix(z.train,nrow = length(lon))

cm <- colormap(zlim = zlim,name="gmt_gebco")
drawPalette(colormap=cm)
mapPlot(longitude=c(-180,180,0,0), latitude=c(0,0,-90,90), projection="+proj=moll", 
        grid=TRUE, col="lightgray",drawBox=FALSE,
        longitudelim=c(-180,180),type="n") # defaults to moll projection
mapImage(longitude=lon*180/pi, latitude=lat*180/pi, z=z.train.mat,zlim=zlim, colormap=cm,missingColor=NA)


# test data
z.test = z.all
z.test[mask.train]=NA
z.test.mat <- matrix(z.test,nrow = length(lon))

cm <- colormap(zlim = zlim,name="gmt_gebco")
drawPalette(colormap=cm)
mapPlot(longitude=c(-180,180,0,0), latitude=c(0,0,-90,90), projection="+proj=moll", 
        grid=TRUE, col="lightgray",drawBox=FALSE,
        longitudelim=c(-180,180),type="n") # defaults to moll projection
mapImage(longitude=lon*180/pi, latitude=lat*180/pi, z=z.test.mat,zlim=zlim, colormap=cm,missingColor=NA)


# prediction
# preds$mu.pred
z.test.pred = z.test
z.test.pred[mask.test]=preds$mu.pred
z.test.pred.mat = matrix(z.test.pred,nrow = length(lon))

cm <- colormap(zlim = zlim,name="gmt_gebco")
drawPalette(colormap=cm)
mapPlot(longitude=c(-180,180,0,0), latitude=c(0,0,-90,90), projection="+proj=moll", 
        grid=TRUE, col="lightgray",drawBox=FALSE,
        longitudelim=c(-180,180),type="n") # defaults to moll projection
mapImage(longitude=lon*180/pi, latitude=lat*180/pi, z=z.test.pred.mat,zlim=zlim, colormap=cm,missingColor=NA)


################################################################################################


