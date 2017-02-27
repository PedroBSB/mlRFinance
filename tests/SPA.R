#Habilita a biblioteca
library(mlRFinance)
#Gera algumas performances relativas aleatórias (LxM)
D<-matrix(runif(10000,-1,+1),ncol=10,nrow = 10)
#Gera a performance relativa do benchmark
b<-runif(1000,-1,+1)
#Call the function
white.spa(Dmat=D,bVec=b,typeFunc=1,B=500,geomMean=20)


#Gera alguns valores de perda para os modelos
model <-runif(1000,-1,1)
bench <-runif(1000,-1,1)
#Diebold Mariano tradicional
DM.epa(e.model=model, e.bench=bench, M=NA)
#Diebold Mariano corrigido
DM.epa.corrected(e.model=model, e.bench=bench, M=NA)

