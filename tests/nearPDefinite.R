library(mlRFinance)
A<-matrix(0,ncol=3,nrow = 3)
A[1,1]<-A[2,2]<-A[3,3]<-1
A[1,2]<-A[2,1]<-3
A[1,3]<-A[3,1]<-2
A[2,3]<-A[3,2]<-4

nearPDefinite(A)


#Example Kernel (SVR)
K<-matrix(c(    4,     9,    36,    49,    -4,    -9,   -36,    -49,
                9,    25,   121,   169,    -9,   -25,  -121,   -169,
                36,   121,   676,   961,   -36,  -121,  -676,   -961,
                49,   169,   961,  1369,   -49,  -169,  -961,   -1369,
                -4,    -9,   -36,   -49,     4,     9,    36,    49,
                -9,   -25,  -121,  -169,     9,    25,   121,   169,
                -36,  -121,  -676,  -961,    36,   121,   676,   961,
                -49,  -169,  -961, -1369,    49,   169,   961,  1369),nrow=8,ncol=8,byrow=T)


nearPDefinite(K)

Matrix::nearPD(K)

#########################################

#Load the data
data("sinc")
#Separate the data
ids <-sample(1:nrow(sinc),30)
train <- as.matrix(sinc[ids,])
#Order the dataset
train <- train[order(train[,1]),]
#Validation
valid <- as.matrix(sinc[-ids,])
#Order the dataset
valid <- valid[order(valid[,1]),]
#Parameters
n <- nrow(train)
C <- 1.0
epsilon <- 0.5
d <- train[,2]
A <- as.matrix(train[,1])
#### Solve by hand
K <- mlRFinance::KernelMatrixComputation(A,"Gaussian",0.5)
# Q matrix
Q <- rbind(cbind(K,-K),cbind(-K,+K))

##
Qnear <- Matrix::nearPD(Q, corr = FALSE, keepDiag = TRUE, do2eigen = TRUE,
                        doSym = FALSE, doDykstra = FALSE, only.values = FALSE,
                        ensureSymmetry = FALSE,
                        eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                        maxit = 100, conv.norm.type = "I", trace = FALSE)
Qnear <- Qnear$mat


#############################################
n <- 3
vecD <- c(-1,-2,-3)
X <-matrix(rnorm(9,100),ncol=3,nrow=3)

vecD*X*rep(vecD, each = n)



for(r in 1:n){
  for(c in 1:n){
    X[r,c] = vecD[r]*X[r,c]*vecD[c];
  }
}

