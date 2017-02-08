#Page 160 - Small
library(mlRFinance)


#Page 160 - Small
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

n <- length(d)
C <- 0.05
epsilon <- 0.01
#### Solve by hand
K <- mlRFinance::KernelMatrixComputation(A,"Gaussian",0.5)
K
# Q matrix
Q <- rbind(cbind(K,-K),cbind(-K,+K))
# Near Positive Definite
Q <- Matrix::nearPD(Q)
Q <- Q$mat
#Dvec
dvec <- c(+d-epsilon,-d-epsilon)
#Constraints Equality
Emat <- c(rep(1,n),rep(-1,n))
#Constraints Inequality
Amat <- cbind(Emat,diag(1,2*n),diag(-1,2*n))
bvec <- c(0, rep(0,2*n),rep(-C,2*n))
#Solve the Quadratic Problem
quadprog::solve.QP(D= Q,dvec =dvec,Amat = Amat, bvec = bvec, meq= 1)



#####
#Dvec
Q <- as.matrix(Q)
dvec <- c(+d-epsilon,-d-epsilon)
#Constraints Equality
Emat <- matrix(c(rep(1,n),rep(-1,n)),ncol=1,nrow=2*n)
ce <- 0
#Constraints Inequality
Amat <- as.matrix(cbind(diag(1,2*n),diag(-1,2*n)))
bvec <- c(rep(0,2*n),rep(-C,2*n))
solveTest(Q, dvec, Amat, bvec, Emat, ce)

#####
A<-matrix(c(1,2,5,6),nrow=4,ncol=1)
d<-c(-1,-1,+1,-1)

n <- length(d)
C <- 0.05
epsilon <- 0.01

svm2<- CSVRL1(d, A, C,epsilon, "Gaussian", c(0.5))
svm2



svm2<- CSVRL1(d, A, 50,0.5, "Polynomial", c(2,1))
svm2



#Create a Linear regression
X<-matrix(rnorm(5000,0,1),ncol=5,nrow=1000)
betas<-matrix(c(-1,2,3,-2,1),ncol=1)
y<- X%*%betas+rnorm(1000,0,1)

#Create the SVR
svm<- CSVRL1(y, X, 0.5,0.1, "Linear", c(0))






Apred<-matrix(c(1.8,3,6,5),nrow=4,ncol=1)
PredictedCSVRL1(svm2, A, A)
R2PredictedCSVRL1(svm2, A)

minimumCSVRL1(d, A, 0.5,"Polynomial", c(2,1))


###################################################################

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
svm2<- CSVRL1(y, X, 0.001,0.01, "Polynomial", c(2,1))

svm3<- CSVRL1(y, X, 1,0.5, "Gaussian", c(1))


#Faz analise de sensibilidade
Q<-minimumCSVRL1(y, X, 0.5, "Gaussian", 1)
Q2<-Matrix::nearPD(Q)
teste<-eigen(Q)

sum(Q2$mat-Q)


solve(Q2$mat)









svm2<- CSVRL1(y, X, 1, 1, "Gaussian", 1)
svm2


svmFn<-function(parms){
  C<-parms[1]
  eps<-parms[2]
  gamma<-parms[3]
  svm2<- CSVRL1(y, X, C, eps, "Gaussian", gamma)
  res<-sum(svm2$SupportVectors)
  if(is.nan(res)){
    return(1e+10);
  }
  else{
    return(exp(res))
  }
}

library(RcppDE)
de<-RcppDE::DEoptim(fn=svmFn, lower=c(0.0001,0.0001,0.0001),
                upper=c(100,100,100),
                control = DEoptim.control(NP = 50, trace = FALSE,
                                          itermax = 1000, F = 1.2, CR = 0.7))


C<-seq(0.001,2,length.out = 1000)
for(i in 1:length(C)){
  svm2<- CSVRL1(y, X, 1, 1, "Gaussian", C[i])
  if(!any(is.na(svm2$SupportVectors)==T)){
    print(svm2$Parameters)
  }
}













