library(rJava)
library(xlsxjars)
library(xlsx)
library(readr)
library(dplyr)
library(car)
library(MASS)
library(forecast)
library(glmnet)
library(Rcpp)
library(optimization)

# 
gene<-read.csv(file = "gene.csv",header = TRUE)

# GMC: input a matrix consists of (X,Y), output a vector
# consists of (GMC(X|Y),GMC(Y|X)).					 
GMC<-function(data,nlocal=25){
  x=data[,1];y=data[,2]
  n=length(x);ln=nlocal
  xdata=data[order(data[,1]),];ydata=data[order(data[,2]),]
  E_xy=rep(0,n);E_yx=rep(0,n)
  X=t(matrix(rep(xdata[,1],n),ncol=n)); X=1/(1+abs(X-t(X)))
  Y=t(matrix(rep(ydata[,2],n),ncol=n)); Y=1/(1+abs(Y-t(Y)))
  for(i in 1:n){
    li=max(1,i-ln)
    ui=min(n,i+ln)
    E_yx[i]=sum(X[i,li:ui]*xdata[li:ui,2])/sum(X[i,li:ui]); 
    E_xy[i]=sum(Y[i,li:ui]*ydata[li:ui,1])/sum(Y[i,li:ui])
  }
  GMC=c(var(E_xy)/var(x),var(E_yx)/var(y))
  return(GMC)
}

#
GMCcalculator<-function(vec){
  data<-cbind(response.variable,scale(vec))
  result<-GMC(data)[1]
}

#
gene.new <- gene[-c(541,542,552,553,398,399,394,395,396,397,
                    408,409,138,139,100,101,163,164,420,421),]

row.names(gene.new) <- 1:nrow(gene.new)

#
mod1<-lm(response~.,data = gene.new)
outlierTest(mod1)
dat1<-gene.new[-480,]
mod2<-lm(response~.,data = dat1)
outlierTest(mod2)
dat2<-dat1[-308,]
mod3<-lm(response~.,data = dat2)
plot(mod3,which = 4)
dat3<-dat2[-316,]
row.names(dat3) <- 1:nrow(dat3)
dat4<-dat3[-449,]

k<-ceiling((544/log(544)))

GMCscreening<-c(0)
response.variable<-scale(dat4$response)

GMCscreening<-apply(scale(dat4[,-1]),MARGIN = 2,FUN = GMCcalculator)

# GMCkeep result for response
GMCkeep<-sort(GMCscreening,decreasing = TRUE)[1:87]

# 87 variables 
sub<-cbind(dat4[,1],dat4[,colnames(dat4)%in%names(GMCkeep)])
colnames(sub)[1]<-"response"

# transformation
sub$response <- BoxCox(sub$response,BoxCox.lambda(sub$response))

sub.scale <- as.data.frame(scale(sub))
model_lasso <- glmnet(as.matrix(sub.scale[,-1]),sub.scale$response,alpha = 1)
par(mfrow = c(1,1))
plot(model_lasso, xvar = "lambda")
crossval <-  cv.glmnet(as.matrix(sub.scale[,-1]),sub.scale$response)
plot(crossval)
penalty <- crossval$lambda.min #optimal lambda
penalty #minimal shrinkage
lasso.fit <- glmnet(as.matrix(sub.scale[,-1]),sub.scale$response, alpha = 1, lambda = penalty ) #estimate the model with that
(result = coef(lasso.fit))

dat5 <- sub.scale[,as.vector(which(as.vector(result) != 0) - 1)[-1]]
dat5$response <- sub.scale$response
mod4 <- lm(response ~ .,data = dat5)
kk <- summary(mod4)
var(kk$residuals)
par(mfrow = c(2,2))
plot(mod4)

######### glm1 x form
xp<-as.data.frame(sub.scale[,2:88])
response.scale<-scale(sub.scale[,1])

MAX<-function(beta){
  fit<-beta[1]+as.matrix(xp)%*%beta[2:88]
  data<-cbind(response.scale,fit)
  -GMC(data)[1]+lambda*sum(abs(beta[2:88]))
}

linearmod<-lm(response.scale~.,data = xp)
par<-linearmod$coefficients

#lambdachoice<-c(0.0001,0.001,0.01,0.1,1,10)
lambdachoice<-seq(0.001,0.1,by=0.001)
for(j in 1:length(lambdachoice)){
  lambda<-lambdachoice[j]
  list<-optim(par,fn=MAX,method ="Nelder-Mead")
  print(c(list$value,lambda))
}
lambda<-0.01
list<-optim(par,fn=MAX,method ="Nelder-Mead")
for(i in 1:88){
  index<-which(abs(list$par)>=sort(abs(list$par),decreasing = TRUE)[i])
  names(index)<-NULL
  beta<-vector("numeric",length = 88)
  beta[index]<-list$par[index]
  beta[1]<-list$par[1]
  fit<-beta[1]+as.matrix(xp)%*%beta[2:88]
  data<-cbind(response.scale,fit)
  print(c(GMC(data),i))
}
for(i in 26){
  index<-which(abs(list$par)>=sort(abs(list$par),decreasing = TRUE)[i])
  print(index)
  print(names(index))
  names(index)<-NULL
  beta<-vector("numeric",length = 88)
  beta[index]<-list$par[index]
  beta[1]<-list$par[1]
  fit<-beta[1]+as.matrix(xp)%*%beta[2:88]
  data<-cbind(response.scale,fit)
  print(c(GMC(data),i))
}
