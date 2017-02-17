#Page 160 - Small
library(mlRFinance)

#Pequeno exemplo
A<-matrix(c(1,2,5,6,
            5,5,2,1,
            8,1,1,7),nrow=4,ncol=3)
svc<-CSVC(A, 1.0, "Gaussian", c(0.5))

#Load the data
data("circle")

#Plot the data
plot(circle)

Xdata<- as.matrix(circle)

#SVC
svc<-CSVC(Xdata, 1.0, "Gaussian", c(0.5))















#Pequeno exemplo
A<-matrix(c(1,2,5,6,
            5,5,2,1,
            8,1,1,7),nrow=4,ncol=3)

svc<-WOCSCM(A, 1, 2, 1, 100, "Gaussian", c(0.5))
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
