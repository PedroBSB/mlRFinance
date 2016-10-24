#Fixa a semente
set.seed(2929)

#Gera algum conjunto de dados
#Numero de variaveis
p<-10
#Numero de observacoes
n<-1000
dados<-matrix(rnorm(p*n),ncol=p,nrow=n)
Kmat<-KernelMatrixComputation(dados,"Gaussian",c(2.5))
