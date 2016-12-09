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
svc<-SpatialWOCSCM(A, wMat, 0.5, 2, 1, 1, 100, "Gaussian", c(0.5))

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









#Pequeno exemplo
A<-matrix(c(1,2,5,6,
            5,5,2,1,
            8,1,1,7),nrow=4,ncol=3)

wMat<-matrix(c(0,1,0,5,
               1,0,0,0,
               0,0,0,2,
               5,0,2,0),byrow = TRUE, nrow=4, ncol=4)


svc<-SpatialWOCSCM(A,wMat, 1, 2, 1, 100, "Gaussian", c(0.5))

svc

#Pequeno grande
data(iris)
df<-as.matrix(iris[,1:4])


#WOCSCM(X, C, k, sigma, inter, kernel, parms)
svc<-WOCSCM(df, 0.05, 2, 1.2, 100, "Gaussian", c(0.5))
svc

head(iris)
teste<-cbind(iris$Species,svc$Zmat)
teste[,2:3]<-round(teste[,2:3])
table(teste[,1],teste[,2])

