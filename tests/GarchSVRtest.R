#MlRFinance
rm(list=ls())
library(mlRFinance)
library(quantmod)
#Cria um novo ambiente para armazenar os dados
stockData <- new.env()

#Especifica as datas de interesse
startDate = as.Date("2011-01-01")
endDate = as.Date("2011-12-31")

#Obtêm os dados do ativo PETR4 e PETR3
getSymbols("^BVSP", src="yahoo",from=startDate,to=endDate)

#Calcula o log-retorno
retorno<-na.omit(diff(log(Cl(BVSP))))

#Training set
train <- as.numeric(retorno[1:180])
#Validation set
valid <- as.numeric(retorno[181:216])
#Cost parameter - Mean Equation
Cmean<-seq(0.01,0.05,length.out = 5)
#Epsilon parameter - Mean Equation
epsilonMean <-seq(0.04,0.1,length.out=7)
#Kernel mean equation
kernelMean <- "Gaussian"
#parameters kernel - Mean Equation
parmsMean <-matrix(c(1.0,
                     1.5,
                     2.0),ncol=1,nrow=3,byrow=T)

#Cost parameter - Volatility Equation
Cvola <-seq(0.5,0.7,length.out = 4)
#Epsilon parameter - Volatility Equation
epsilonVola <-seq(0.01,0.1,length.out=4)
#Kernel - Volatility Equation
kernelVolat <- "Polynomial"
#parameters kernel - Volatility Equation
parmsVola <- matrix(c(2,1,
                      3,5,
                      4,7),ncol=2,nrow=3,byrow=T)

#Do the cross-validation

train<-train
valid<-valid
Cm<-Cmean
epsilonM<-epsilonMean
kernel<-kernelMean
parmMat<-parmsMean
Cg<-Cvola
epsilonG<-epsilonVola
kernelGarch<-kernelVolat
parmMatGarch<-parmsVola

ptm <- proc.time()
  teste<-GarchSVR(train,valid,Cmean,epsilonMean,kernelMean,parmsMean,
                  Cvola,epsilonVola,kernelVolat,parmsVola)
proc.time() - ptm

