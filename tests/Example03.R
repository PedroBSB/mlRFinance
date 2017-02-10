#Load the data
data("sinc")
set.seed(3838)
#Plot the data
plot(sinc[,1],sinc[,2],type="l")

#Separate the data
ids <-sample(1:nrow(sinc),30)
train <- as.matrix(sinc[ids,])

#Order the dataset
train <- train[order(train[,1]),]

#Validation
valid <- as.matrix(sinc[-ids,])

#Order the dataset
valid <- valid[order(valid[,1]),]


parms<-matrix(rep(0.5,3),nrow=3)

#Do the SVR
svm2<- CSVRL1(train[,2], as.matrix(train[,1]), 0.5, 0.05, "Gaussian", parms[1,])
#svm2<- CSVRL1(train[,2], as.matrix(train[,1]), 0.5,0.05, "Polynomial", c(2,1))
svm2

#Do the forecast
pred<-PredictedCSVRL1(svm2,train[,2], as.matrix(train[,1]), as.matrix(train[,1]))

#Sort dataset
plot(train[,1],train[,2],type="l",col="blue")
lines(train[,1],pred,col="red")

#Exercise: Predict with the validation data.set

#Exercise: Change the kernel

#Exercise: Improve the results


