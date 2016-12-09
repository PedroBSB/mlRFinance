#Page 160 - Small
library(mlRFinance)
library(ggplot2)

#Read the data frame
data("Crime")

#Read the shapefile
shp <- maptools::readShapePoly("data/columbus.shp", IDvar="POLYID")

#Create the weight matrix
nb <- spdep::poly2nb(shp)
wMat <- spdep::nb2mat(nb, style="W", zero.policy =TRUE)

#Data matrix
A<-as.matrix(scale(Crime[,c(2:4)]))

#Run the SpatialWOCSCM
svc<-SpatialWOCSCM(A, wMat, 0.1, 2, 0.01, 1, 100, "Gaussian", c(0.5))

#Plot the log-likelihood
plot(svc$LogLikelihood, type="l")

#Get the probability matrix
iClasse <- apply(svc$Zmat,1,function(x) which(x==max(x)[1]))

#Create the thematic map
Crime$Cluster<-as.factor(iClasse)

#Spatial Polygon Data Frame
shp2 <- sp::merge(shp, Crime, by.x="POLYID", by.y="ID")

#Colors
col <- rainbow(length(levels(Crime$Cluster)))

#Make tha map
sp::spplot(shp2, "Cluster", col.regions=col, main="Clusters",
       colorkey = FALSE, lwd=.4, col="white")




