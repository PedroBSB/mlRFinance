#Load the data
data("SpiralDataSet")

#Plot the data
colors<-ifelse(SpiralDataSet[,3]==-1,"red","blue")
plot(SpiralDataSet[,1:2],col=colors,pch=16)

#Separate the data
train <- as.matrix(SpiralDataSet[1:100,])
X_train <- train[,1:2]
y_train <- train[,3]

valid <- as.matrix(SpiralDataSet[101:194,])
X_valid <- valid[,1:2]
y_valid <- valid[,3]

#Automatic cross-validation
cross<- PortfolioSelectionCSVML1(y_train, X_train,
                                 y_valid, X_valid,
                                 1,
                                 "Polynomial", c(2,1), 0)

#Exercise: Change the kernel

#Exercise: Improve the results

