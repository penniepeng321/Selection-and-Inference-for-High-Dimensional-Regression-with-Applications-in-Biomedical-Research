rm(list=ls())
setwd("/mnt/home/hpeng2")
library(lattice) 
library(scalreg)
library(glmnet)
library(glasso)
library(camel)
library(clime)
library(hdi)
library(mvtnorm)
library(Matrix)


pvalue_de=function(X,y){
  n=nrow(X)
  p=ncol(X)
  model=glmnet(X,y,family="gaussian")
  linda=cv.glmnet(X,y,nfolds = 10)
  cvpenalty=linda$lambda.min
  betahat=coef(model,s=cvpenalty)[-1]
  
  temp=X%*%betahat
  
  samplecov=matrix(0,((p+1)/2),((p+1)/2))
  for(cx in 1:n){
    samplecov=samplecov+(X[cx,((p+1)/2):p])%*%t(X[cx,((p+1)/2):p])*(y[cx]-temp[cx])^2
  }
  samplecov=samplecov/n
  
  ## Precision matrix
  temp=camel.tiger(X[,((p+1)/2+1):p], lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL,
                   method = "slasso", sym = "or", shrink=NULL, prec = 1e-4, mu = 0.01,
                   max.ite = 1e4, standardize = T, correlation = FALSE,
                   perturb = TRUE, verbose = F)
  temp1=matrix(as.numeric(unlist(temp$icov[1])),((p+1)/2-1),((p+1)/2-1))
  precisionmatrix=as.matrix(bdiag(4,temp1))
  
  aspvarmatrix=precisionmatrix%*%samplecov%*%t(precisionmatrix)
  sigmadiag=diag(aspvarmatrix)
  biasdelta=rep(0,((p+1)/2))
  bhat=betahat[((p+1)/2):p]+(precisionmatrix%*%(t(X[,((p+1)/2):p]))%*%(y-X%*%betahat))/n
  bhat[bhat==-Inf]=0
  bhat[bhat==Inf]=0
  bhat[is.nan(bhat)]=0
  pvalue_de=rep(NA,((p+1)/2))
  for(k in 1:((p+1)/2)){
    if(sigmadiag[k]<=0 & bhat[k]==0){
      pvalue_de[k]=1
    }
    if(sigmadiag[k]<=0 & bhat[k]!=0){
      pvalue_de[k]=0
    }
    if(sigmadiag[k]>0){
      pvalue_de[k]=(1-pnorm(abs(bhat[k]/(sqrt(sigmadiag[k])/sqrt(n))),0,1))*2
    }
  }
  order_pvalue=sort(pvalue_de,decreasing = F,index.return=T)$x
  index_pvalue=sort(pvalue_de,decreasing = F,index.return=T)$ix
  err_var=1/(n-1)*sum((y-X%*%betahat)^2)
  
  upperci=rep(NA,((p+1)/2))
  lowerci=rep(NA,((p+1)/2))
  lengthci=rep(NA,((p+1)/2))
  for(k in 1:((p+1)/2)){
    if(sigmadiag[k]<=0){
      lowerci[k]=bhat[k]
      upperci[k]=bhat[k]
      lengthci[k]=0
    }
    if(sigmadiag[k]>0){
      lowerci[k]=bhat[k]-1.96*sqrt(sigmadiag[k])/sqrt(n)
      upperci[k]=bhat[k]+1.96*sqrt(sigmadiag[k])/sqrt(n)
      lengthci[k]=upperci[k]-lowerci[k]
    }
  }
  
  pvalue_delist=list(upperci=upperci,lowerci=lowerci,lengthci=lengthci,
                     pvalue=order_pvalue,sequence=index_pvalue,
                     error_variance=err_var,betahat=betahat,bhat=bhat)
  return(pvalue_delist)
}




n=200
p=200

############################################################
sim=1000
Vm=rep(0,sim)
set.seed(9871)
for(hj in 1:sim){
  temp=runif(p,0,1)
  Um=sort(temp,decreasing = F)
  Vm[hj]=max((((1:p)/(p)-Um)/sqrt(Um*(1-Um)))[1:round(p/2)])
}
Vm=sort(Vm)


pro=0.5
rhocorr=0.5
covmat=rhocorr^{abs(outer(1:p,1:p,'-'))}
gamma_par=c(rep(1,4),rep(0,p-4))
beta_par=c(1,rep(2,4),rep(0,p-4))
signal_index=c(2:5)
s=4

#simulation starts here!
sim_size=30
BIC=rep(0,sim_size)
sBIC=rep(0,sim_size)
LCV=rep(0,sim_size)
sLCV=rep(0,sim_size)
BON=rep(0,sim_size)
sBON=rep(0,sim_size)
FDR=rep(0,sim_size)
sFDR=rep(0,sim_size)
NPE=rep(0,sim_size)
sNPE=rep(0,sim_size)
set.seed(3241)
for(k in 1:sim_size){
  
  
  X = rmvnorm(n, mean = rep(0, nrow(covmat)), sigma = covmat,method=c("chol"))
  X_intercept=cbind(rep(1,n),X)
  
  trt=rbinom(n,size=1,prob=pro)
  trt_X_intercept=trt*X_intercept
  
  my=1+X%*%gamma_par+trt_X_intercept%*%beta_par
  y=rep(NA,n)
  for (hfg in 1:n){
    y[hfg]=rnorm(1,my[hfg],1)
  }
  
  newX=cbind(X,(trt-pro)*X_intercept)
  modelfit=pvalue_de(newX,y)
  
  order_pvalue=modelfit$pvalue
  
  #######################################BIC cut
  modelr=glmnet(newX,y)
  betamatrix=modelr$beta
  zeroone=matrix(as.numeric(betamatrix!=0),nrow(betamatrix),ncol(betamatrix))
  action=0
  for(oppo in 1:(ncol(betamatrix)-1)){
    temp=(1:(2*p+1))[(zeroone[,oppo+1]-zeroone[,oppo])!=0]
    if(length(temp)>0){
      action=append(action,temp)
    }
  }
  action=action[-1]
  knotr=modelr$lambda
  lenr=length(knotr)
  BICc=rep(0,lenr+1)
  sigmag=modelfit$error_variance
  BICc[1]=sum((y-mean(y))^2)/(sigmag)
  intX=cbind(rep(1,n),newX)
  for(lp in 2:(lenr+1)){
    temp=matrix(as.numeric(unlist(coef(modelr,s=knotr[lp-1]))),2*p+2,1)
    BICc[lp]=(1/(sigmag))*sum((y-intX%*%temp)^2)+lp*log(2*p+2)
  }
  BIC[k]=which.min(BICc)
  if(BIC[k]>=lenr){
    BIC[k]=length(action[action>p])
    sBIC[k]=s-length(setdiff(signal_index,action-p))
  }
  if(BIC[k]>0 & BIC[k]<lenr){
    temp=action[1:BIC[k]]
    select_set_BIC=temp[temp>p]
    BIC[k]=length(select_set_BIC)
    sBIC[k]=s-length(setdiff(signal_index,select_set_BIC-p))
  }
  if(BIC[k]==0){
    sBIC[k]=0
  }
  
  
  ###############################lasso cross validation
  selected_set=lasso.cv(newX,y)
  LCV[k]=length(selected_set)
  if(LCV[k]>=lenr){
    LCV[k]=length(action[action>p])
    sLCV[k]=s-length(setdiff(signal_index,action-p))
  }
  if(LCV[k]>0 & LCV[k]<lenr){
    temp=action[1:LCV[k]]
    select_set_LCV=temp[temp>p]
    LCV[k]=length(select_set_LCV)
    sLCV[k]=s-length(setdiff(signal_index,select_set_LCV-p))
  }
  if(LCV[k]==0){
    sLCV[k]=0
  }
 
  
  ###################################### Bon cut
  BON[k]=sum(((order_pvalue-0.05/(p+1))<0)==T)
  
  if(BON[k]>0){
    sBON[k]=s-length(setdiff(signal_index,modelfit$sequence[1:BON[k]]))
  }
  if(BON[k]==0){
    sBON[k]=0
  }
  
  ###################################### FDR cut
  FDR[k]=sum(((order_pvalue-0.05*(1:(p+1))/(p+1))<0)==T)
  
  if(FDR[k]>0){
    sFDR[k]=s-length(setdiff(signal_index,modelfit$sequence[1:FDR[k]]))
  }
  if(FDR[k]==0){
    sFDR[k]=0
  }
  
  ###################################### New procedure
  #bounding sequence
  ppp=1-1/sqrt(log(p+1))
  Vm=as.numeric(unlist(Vm))
  cbd=as.numeric(quantile(Vm,max(ppp,0)))
  lpm=p+1
  qts1=order_pvalue
  NPE[k]=max(ceiling(max(((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1)[1:round(lpm/2)])*lpm),
             0)
  
  if(NPE[k]>0){
    sNPE[k]=s-length(setdiff(signal_index,modelfit$sequence[1:NPE[k]]))
  }
  if(NPE[k]==0){
    sNPE[k]=0
  }
  
}

c(mean(BIC,na.rm=TRUE),sd(BIC,na.rm=TRUE),mean(sBIC,na.rm=TRUE),sd(sBIC,na.rm=TRUE),mean(sBIC/s, na.rm=TRUE),sd(sBIC/s, na.rm=TRUE),mean((BIC-sBIC)/BIC, na.rm=TRUE),sd((BIC-sBIC)/BIC, na.rm=TRUE))

c(mean(LCV,na.rm=TRUE),sd(LCV,na.rm=TRUE),mean(sLCV,na.rm=TRUE),sd(sLCV,na.rm=TRUE),mean(sLCV/s, na.rm=TRUE),sd(sLCV/s, na.rm=TRUE),mean((LCV-sLCV)/LCV, na.rm=TRUE),sd((LCV-sLCV)/LCV, na.rm=TRUE))

c(mean(BON,na.rm=TRUE),sd(BON,na.rm=TRUE),mean(sBON,na.rm=TRUE),sd(sBON,na.rm=TRUE),mean(sBON/s, na.rm=TRUE),sd(sBON/s, na.rm=TRUE),mean((BON-sBON)/BON, na.rm=TRUE),sd((BON-sBON)/BON, na.rm=TRUE))

c(mean(FDR,na.rm=TRUE),sd(FDR,na.rm=TRUE),mean(sFDR,na.rm=TRUE),sd(sFDR,na.rm=TRUE),mean(sFDR/s, na.rm=TRUE),sd(sFDR/s, na.rm=TRUE),mean((FDR-sFDR)/FDR, na.rm=TRUE),sd((FDR-sFDR)/FDR, na.rm=TRUE))

c(mean(NPE,na.rm=TRUE),sd(NPE,na.rm=TRUE),mean(sNPE,na.rm=TRUE),sd(sNPE,na.rm=TRUE),mean(sNPE/s, na.rm=TRUE),sd(sNPE/s, na.rm=TRUE),mean((NPE-sNPE)/NPE, na.rm=TRUE),sd((NPE-sNPE)/NPE, na.rm=TRUE))

############################################################
temp=matrix(c(BIC,sBIC,LCV,sLCV,BON,sBON,FDR,sFDR,NPE,sNPE),nrow=sim_size,ncol=10,byrow=F)
write.table(temp,file="sim13.csv",row.names=F)

mean((BIC-sBIC)/(rep(p+1-s,sim_size)),na.rm=TRUE);sd((BIC-sBIC)/(rep(p+1-s,sim_size)),na.rm=TRUE)
mean((LCV-sLCV)/(rep(p+1-s,sim_size)),na.rm=TRUE);sd((LCV-sLCV)/(rep(p+1-s,sim_size)),na.rm=TRUE)
mean((BON-sBON)/(rep(p+1-s,sim_size)),na.rm=TRUE);sd((BON-sBON)/(rep(p+1-s,sim_size)),na.rm=TRUE)
mean((FDR-sFDR)/(rep(p+1-s,sim_size)),na.rm=TRUE);sd((FDR-sFDR)/(rep(p+1-s,sim_size)),na.rm=TRUE)
mean((NPE-sNPE)/(rep(p+1-s,sim_size)),na.rm=TRUE);sd((NPE-sNPE)/(rep(p+1-s,sim_size)),na.rm=TRUE)


upperci=modelfit$upperci
lowerci=modelfit$lowerci
lengthci=modelfit$lengthci
betahat=modelfit$betahat
bhat=modelfit$bhat

mean((bhat-beta_par)^2,na.rm=T)
mean((betahat-beta_par)^2,na.rm=T)

mean(upperci[1:s]-lowerci[1:s],na.rm=T)
mean((upperci[1:s]>=beta_par[1:s])&(lowerci[1:s]<=beta_par[1:s]),na.rm=T)
mean(upperci[(s+1):(p+1)]-lowerci[(s+1):(p+1)],na.rm=T)
mean((upperci[(s+1):(p+1)]>=beta_par[(s+1):(p+1)])&(lowerci[(s+1):(p+1)]<=beta_par[(s+1):(p+1)]),
     na.rm=T)

jpeg(paste('s',s,rhocorr,'sim13.jpg'))
plot(upperci,type='l',col='grey',ylab='coef',xlab='predictor index',
     ylim=c(-4,4),
     main=paste('n=',n,';p=',p,';s=',s))
lines(lowerci,col='green')
lines(beta_par,col='black',lwd=2)
legend('topright',c('upper','lower','est'),lwd=c(1,1,2),col=c('grey','green',
                                                              'black'))
dev.off()


modelfit
BIC;sBIC
LCV;sLCV
BON;sBON
FDR;sFDR
NPE;sNPE
action