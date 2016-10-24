#MlRFinance
library(mlRFinance)

#Habilita o pacote quantmod
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

train <- as.numeric(retorno[1:180])
valid <- as.numeric(retorno[181:216])
Cmean <- 0.01
epsilonMean <- 0.04
kernelMean <- "Gaussian"
parmsMean <- 1

Cvola <- 0.5
epsilonVola <- 0.01
kernelVolat <- "Polynomial"
parmsVola <- c(2,1)

svm<-GARCHCSVRL1(train, valid, Cmean, epsilonMean, Cvola, epsilonVola, kernelMean, parmsMean, kernelVolat, parmsVola)
svm$PredictedMeanTraining
