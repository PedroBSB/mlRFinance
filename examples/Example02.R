#Load the data
data("SpiralDataSet")

#Plot the data
colors<-ifelse(SpiralDataSet[,3]==-1,"red","blue")
plot(SpiralDataSet[,1:2],col=colors,pch=16)

#Sample the data
ids<-sample(1:nrow(SpiralDataSet),100,F)

#Separate the data
train <- as.matrix(SpiralDataSet[ids,])
valid <- as.matrix(SpiralDataSet[-ids,])

#Train
svm1 <- CSVML1(train[,3], train[,1:2], 50, "Polynomial", c(2,1), FALSE)

#Predict values
pred <- PredictedCSVML1(svm1, train[,3], train[,1:2], train[,1:2], 0, FALSE)

#Plot the prediction
colors<-ifelse(pred==-1,"red","blue")
plot(train[,1:2],col=colors,pch=16)

#Confusion matrix
table(train[,3],pred)

#Predict based on the valid dataset
predValid <- PredictedCSVML1(svm1, train[,3], train[,1:2], valid[,1:2], 0, FALSE)

#Confusion matrix
table(valid[,3],predValid)

#Plot the results
colors<-ifelse(predValid==-1,"red","blue")
plot(valid[,1:2],col=colors,pch=16)

#Exercise: Change the kernel


