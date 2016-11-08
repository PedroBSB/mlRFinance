#Generate the sample
Dmat<-matrix(runif(10000,ncol=10,nrow=1000),nrow=1000,ncol=10)
bvec<-runif(1000)
typeFunc<-1
B<-500
geomMean<-20
bandwidth<-0.5
alpha<-0.05
k<-1
gamma<-0.1
tt<-hansen.spa(Dmat,bVec,typeFunc=1,B=1000,geomMean=20,bandwidth=0.5, alpha=0.05, k=1, gamma=0.1)
fdp<-tt$FDP
fdp<-tt$FWERk



#####################################################################################################
