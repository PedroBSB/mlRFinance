#Habilita a biblioteca
library(mlRFinance)
#Carrega os dados
data("sinc")
#Fixa a semente
set.seed(3838)
#Plota os dados
plot(sinc[,1],sinc[,2],type="l",
     xlab="Variável independente", ylab="Variável dependente")

#Separa os dados em dois grupos
ids <-sample(1:nrow(sinc),30)
train <- as.matrix(sinc[ids,])
#Ordena os dados
train <- train[order(train[,1]),]
#Conjunto de validação
valid <- as.matrix(sinc[-ids,])
#Ordena os dados
valid <- valid[order(valid[,1]),]

#Define quais são as variáveis dependentes
train.y<- train[,2]
valid.y<- valid[,2]
#Define quais são as variáveis independentes
train.X <- as.matrix(train[,1])
valid.X <- as.matrix(valid[,1])
#Lista de custos
C <- seq(2^-15,2^15, length.out = 10)
#Lista de epsilon's
epsilon <- c(0.1,0.2,0.5)
#Tipo do kernel
kernel<-"Gaussian"
#Matriz com os parâmetros do kernel
parmMat<-matrix(c(0.1,0.2,0.3),nrow=3,ncol=1)
#Ensina a máquina
res<-LearningSVRL1(train.y, train.X,
                   valid.y, valid.X,
                   C, epsilon, kernel, parmMat)
