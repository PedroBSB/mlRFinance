#Page 160 - Small
library(mlRFinance)
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

rank<-order(d)
svm2<- CSVWQR(d,rank, A, 50,0.5 ,0.5, "Polynomial", c(2,1))
svm2

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

#Cria os objetos
y<-train[2:length(train)]
X<-matrix(train[1:(length(train)-1)],ncol=1)

#SVR
rank<-order(y)
svm2<- CSVWQR(y, rank, X, 1, 0.5, 1.0, "Gaussian", c(0.5))
