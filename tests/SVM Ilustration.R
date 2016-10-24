library(ggplot2);library(kernlab)
n=30; x1=matrix(runif(2*n), n, 2); y1=rep(1,n); x2=1.3 + matrix(runif(2*n),n,2); y2=rep(-1,n);
data=data.frame(x=c(x1[,1], x2[,1]), y=c(x1[,2], x2[,2]), class=factor(c(y1,y2)))
head(data); tail(data)


ggplot(data, aes(x,y, color=class)) +
  geom_point(size=3) +
  geom_abline(intercept=2, slope=-0.5) +
  geom_abline(intercept=2.3, slope=-0.9) + geom_abline(intercept=1.4, slope=-0.3) +                   geom_abline(intercept=2.3, slope=-0.8) + geom_abline(intercept=2.1, slope=-0.7)


svm_model= ksvm(class ~ x+ y, data=data, kernel='vanilladot')

##  Setting default kernel parameters
alpha(svm_model)

coef(svm_model)
b(svm_model)


w = colSums(coef(svm_model)[[1]] * data[alphaindex(svm_model)[[1]],c('x','y')])
b = b(svm_model)
data[, 1:2] = sapply(data[,1:2], scale)

ggplot(data,aes(x, y, color=class)) +
  geom_point(size=3) + geom_abline(intercept=b/w[1], slope=-w[2]/w[1]) +
  geom_abline(intercept=(b+1)/w[1], slope=-w[2]/w[1], linetype=2)


plot(svm_model, data=data)


##Redo with smaller cost factor C=0.1
svm_model= ksvm(class ~ x+ y, data=data, kernel='vanilladot')


##  Setting default kernel parameters
alpha(svm_model); coef(svm_model); b(svm_model)


w = colSums(coef(svm_model)[[1]] * data[alphaindex(svm_model)[[1]],c('x','y')]); b = b(svm_model)
data[, 1:2] = sapply(data[,1:2], scale)
plot(svm_model, data=data)
