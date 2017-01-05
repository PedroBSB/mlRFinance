#Loda data
data("circle")
#Plot data
plot(circle)
#Compute the Kernel PCA
kpca<-KPCA(circle, "Gaussian", 0.5)
#Plot variability
plot(kpca, type = "l")
# summary method
summary(kpca)
#Get the first two components
Loadings <- as.data.frame(kpca$rotation[,1:2])
plot(Loadings)
#Separate the data
Loadings$Type<-ifelse(Loadings$PC1< -0.02,1,2)
plot(Loadings$PC1,Loadings$PC2,col=Loadings$Type)
circle<-as.data.frame(circle)
circle$Type<-ifelse(Loadings$PC1< -0.02,1,2)
plot(circle[,1],circle[,2],col=circle$Type)

#Faz a projecao
X2<-matrix(rnorm(5000),nrow=1000,ncol=5)
kpca2<-KPCApredict(kpca,X2,"Gaussian", 0.5)






Xcov<-cov(X)
test<-PCA(Xcov, "Gaussian", 0.5)

e.cov<-eigen(cov(X))

eigenvalues<-e.cor$values
eigenvectors<-e.cor$vectors



X2=scale(X) #New scaled dataset with mean of 0 and sd of 1
scores<-X2%*%eigenvectors #New Variables
total.var<-sum(diag(cov(X2))) #Calculate total variance in scaled data
prop.var<-rep(NA,ncol(X));cum.var<-rep(NA,ncol(X)) #Create empty vectors
#Calculate proportion of variance explained and cumulative variance explained
for(i in 1:ncol(X)){prop.var[i]<-var(scores[,i])/total.var}
for(i in 1:ncol(X)){cum.var[i]<-sum(prop.var[1,i])
sdev=sqrt(eigenvectors) #Std of the new components




ir.pca <- prcomp(Xcov,
                 center = TRUE,
                 scale. = TRUE)





print(ir.pca)
test$Componentes

kpca<-KPCA(X, "Gaussian", 0.5)
tt<-kpca$Componentes
plot(kpca$VariabilityExplained)
kpca$TotalVariability
