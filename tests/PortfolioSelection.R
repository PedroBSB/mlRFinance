#Read data
data("GermanBank")

#Sample the data
set.seed(3831)
ids<-sample(1:nrow(GermanBank),700,F)

#Split between training, validation
train <- as.matrix(GermanBank[ids,])
valid <- as.matrix(GermanBank[-ids,])

train.y<- ifelse(train[,32]==1,+1,-1)
valid.y<- ifelse(valid[,32]==1,+1,-1)

train.X <- as.matrix(train[,3:11])
valid.X <- as.matrix(valid[,3:11])

C <- seq(2^-15,2^15, length.out = 10)
kernel<-"Gaussian"
parmMat<-matrix(c(0.1,0.2,0.3),nrow=3,ncol=1)
typePredict <- 0

res<-LearningSVML1(train.y,train.X, valid.y, valid.X, C, kernel, parmMat, typePredict)

##############################################################

#Read data
data("MackeyGlass")

#Lagging the data frame
library(DataCombine)
DataSlid1 <- slide(MackeyGlass, Var = "V1", slideBy = c(-1,-2,-3,-4), reminder=F)
DataSlid1<-DataSlid1[-1:-4,]

#Split between training, validation
train <- as.matrix(DataSlid1[1:1000,])
valid <- as.matrix(DataSlid1[1001:1496,])

train.y<- train[,1]
valid.y<- valid[,1]

train.X <- as.matrix(train[,2:5])
valid.X <- as.matrix(valid[,2:5])

C <- seq(2^-15,2^15, length.out = 10)
epsilon <- c(0.1,0.2,0.5)
kernel<-"Gaussian"
parmMat<-matrix(c(0.1,0.2,0.3),nrow=3,ncol=1)


res<-LearningSVRL1(train.y,train.X, valid.y, valid.X, C, epsilon, kernel, parmMat, typePredict)
