X<-matrix(rnorm(5000),nrow=1000,ncol=5)

kpca<-KPCA(X, "Gaussian", 0.5)
tt<-kpca$Componentes
plot(kpca$VariabilityExplained)
kpca$TotalVariability
