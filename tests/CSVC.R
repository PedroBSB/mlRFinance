#Page 160 - Small
library(mlRFinance)

#Pequeno exemplo
A<-matrix(c(1,2,5,6,
            5,5,2,1,
            8,1,1,7),nrow=4,ncol=3)
svc<-CSVC(A, 1.0, "Gaussian", c(0.5))

#Load the data
data("circle")

#Plot the data
plot(circle)

Xdata<- as.matrix(circle)

#Generate a sample
ids<-sample(1:nrow(Xdata),100,F)
Xdata <- Xdata[ids,]
plot(Xdata)

#SVC
svc<-CSVC(Xdata, 1.0, "Gaussian", c(0.5))

A<-svc$AdjacencyMatrix
image(A)


#use igraph for example
library(igraph)
#plot it
graph <- graph.adjacency(A, mode = "undirected")
plot(graph,  layout=layout.fruchterman.reingold)
#Community
wc <- fastgreedy.community(graph)
#wc <- edge.betweenness.community(graph)
#wc <- spinglass.community(graph)
#wc <- leading.eigenvector.community(graph)
#wc <- label.propagation.community(graph)
#wc <- walktrap.community(graph)
#Plot community
plot(graph, vertex.label=NA, vertex.size=5,edge.width=0.5,
     vertex.color=membership(wc), layout=layout.fruchterman.reingold)
#Colors
col <- membership(wc)
plot(Xdata, col=col, pch = 19)

#Pequeno exemplo
A<-matrix(c(1,2,5,6,
            5,5,2,1,
            8,1,1,7),nrow=4,ncol=3)

svc<-WOCSCM(A, 1, 2, 1, 100, "Gaussian", c(0.5))
svc

#Pequeno grande
data(iris)
df<-as.matrix(iris[,1:4])


#WOCSCM(X, C, k, sigma, inter, kernel, parms)
svc<-WOCSCM(df, 0.05, 2, 1.2, 100, "Gaussian", c(0.5))
svc

head(iris)
teste<-cbind(iris$Species,svc$Zmat)
teste[,2:3]<-round(teste[,2:3])
table(teste[,1],teste[,2])
