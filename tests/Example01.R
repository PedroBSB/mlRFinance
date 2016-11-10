#Load the data
data("circle")

#Plot the data
plot(circle)

#SVC
svc<-WOCSCM(circle, 0.05, 2, 1.2, 100, "Gaussian", c(0.5))
svc

#Plot the log-likelihood
plot(svc$LogLikelihood,type="l")

#Define the class
class<-apply(svc$Zmat,1,function(x) which(x==max(x)))

#Join with the data
res<-cbind(circle,class)

#Plot the results
plot(res[,1:2],col=res[,3])

#Exercise: Change the kernel

#Exercise: Improve the results
