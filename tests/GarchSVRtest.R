library(mlRFinance)
kernel <- "Gaussian"
kernelGarch <- "Polynomial"
Cm<-seq(1,10)
Cg<-seq(3,7)
epsilonM<-seq(0.1,0.5,length.out=7)
epsilonG<-seq(0.4,0.9,length.out=4)

parmMat<-matrix(c(1.0,2,
                  1.5,5,
                  2.0,7),ncol=2,nrow=3,byrow=T)

parmMatGarch<-matrix(c(0.5,0.3,
                       0.1,0.4),ncol=2,nrow=2,byrow=T)


teste<-GarchSVR(train,valid,Cm,epsilonM,kernel,parmMat,Cg,epsilonG,kernelGarch,parmMatGarch)
