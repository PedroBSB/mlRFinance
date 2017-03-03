#Habilita a biblioteca
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

train <- as.numeric(retorno[1:216])

#Cria os objetos
y<-train[2:length(train)]
X<-matrix(train[1:(length(train)-1)],ncol=1)

#Define quais são as variáveis dependentes
train.y<- y

#Define quais são as variáveis independentes
train.X <- X

#Lista de custos
C <- seq(2^-15,2^15, length.out = 10)
#Lista de epsilon's
tau <- 0.025
gamma <- seq(2^-15,2^15, length.out = 10)

#Tipo do kernel
kernel<-"Gaussian"
#Matriz com os parâmetros do kernel
parmMat<-matrix(seq(2^-15,2^15, length.out = 10),nrow=10,ncol=1)
#Ensina a máquina
res<-LearningSVWQR1(train.y, train.X,
                   C, tau, gamma, kernel, parmMat)
