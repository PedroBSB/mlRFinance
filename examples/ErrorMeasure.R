set.seed(9292)
y<-rnorm(1000)^2
yPred<-rnorm(1000)^2

num<-sum(((y[2:1000]-yPred[2:1000])^2)/y[2:1000])
den<-sum(((y[1:999]-y[2:1000])^2)/y[1:999])
sqrt(num/den)

TheilUfunction(y, yPred)
