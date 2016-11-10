#Load the data
data("SpiralDataSet")

#Plot the data
colors<-ifelse(SpiralDataSet[,3]==-1,"red","blue")
plot(SpiralDataSet[,1:2],col=colors,pch=16)

#Separate the data
train <- as.matrix(SpiralDataSet[1:100,])
valid <- as.matrix(SpiralDataSet[101:194,])

#Train
svm1 <- CSVML1(train[,3], train[,1:2], 50, "Polynomial", c(2,1))

#Predict values
pred <- PredictedCSVML1(svm1, train[,3], train[,1:2], train[,1:2], "Polynomial", c(2,1), 0)

#Plot the prediction
colors<-ifelse(pred==-1,"red","blue")
plot(train[,1:2],col=colors,pch=16)

#Confusion matrix
table(train[,3],pred)

#Predict based on the valid dataset
predValid <- PredictedCSVML1(svm1, train[,3], train[,1:2], valid[,1:2], "Polynomial", c(2,1), 0)

#Confusion matrix
table(valid[,3],predValid)

#Plot the results
colors<-ifelse(predValid==-1,"red","blue")
plot(valid[,1:2],col=colors,pch=16)

#Exercise: Change the kernel


