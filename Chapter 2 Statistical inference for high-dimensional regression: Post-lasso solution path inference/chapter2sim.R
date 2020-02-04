#rm(list=ls())
#setwd("/mnt/home/hpeng2")
library(lattice) 
library(scalreg)
library(glmnet)
library(glasso)
library(camel)
library(clime)
library(hdi)
library(mvtnorm)


pvalue_de=function(X,y){
  n=nrow(X)
  p=ncol(X)
  model=scalreg(X,y)
  betahat=model$coefficients
  sighat=model$hsigma
  barX=colMeans(X)
  samplecov=matrix(0,p,p)
  for(cx in 1:n){
    samplecov=samplecov+t(X[cx,]-matrix(barX,1,p))%*%(X[cx,]-matrix(barX,1,p))
  }
  samplecov=samplecov/(n-1)
  
  #adaptive thresholding to estimate sample covariance matrix
  #   theta_adj=matrix(0,p,p)
  #   for(op1 in 1:p){
  #     for(op2 in 1:p){
  #       for(op3 in 1:n){
  #         theta_adj[op1,op2]=theta_adj[op1,op2]+((X[op3,op1]-barX[op1])*
  #                                                  (X[op3,op2]-barX[op2])-
  #           samplecov[op1,op2])^2
  #       }
  #       theta_adj[op1,op2]=theta_adj[op1,op2]/n
  #     }
  #   }
  
  theta_adj=matrix(0,p,p)
  for(op1 in 1:p){
    for(op2 in 1:p){
      theta_adj[op1,op2]=sum(((X[,op1]-barX[op1])*(X[,op2]-barX[op2])-
                                samplecov[op1,op2])^2)/n
    }
  }
  
  lamda_adj=2*sqrt(theta_adj*log(p)/n)
  
  samplecov_adj=samplecov*(samplecov>lamda_adj)
  
  
  ## Precision matrix
  temp=camel.tiger(X, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL,
                   method = "slasso", sym = "or", shrink=NULL, prec = 1e-4, mu = 0.01,
                   max.ite = 1e4, standardize = FALSE, correlation = FALSE,
                   perturb = TRUE, verbose = F)
  temp1=matrix(as.numeric(unlist(temp$icov[1])),p,p)
  
  precisionmatrix=temp1
  aspvarmatrix=(sighat^2)*precisionmatrix%*%samplecov%*%t(precisionmatrix)
  sigmadiag=diag(aspvarmatrix)
  #biasdelta=sqrt(n)*((precisionmatrix%*%samplecov_adj-diag(p))%*%betahat)
  biasdelta=rep(0,p)
  bhat=betahat+precisionmatrix%*%(t(X)%*%(y-X%*%betahat))/n
  bhat[bhat==-Inf]=0
  bhat[bhat==Inf]=0
  bhat[is.nan(bhat)]=0
  pvalue_de=rep(NA,p)
  for(k in 1:p){
    if(sigmadiag[k]<=0){
      pvalue_de[k]=0
    }
    if(sigmadiag[k]>0){
      pvalue_de[k]=(1-pnorm(abs(bhat[k]/(sqrt(sigmadiag[k])/sqrt(n))),0,1))*2
    }
  }
  order_pvalue=sort(pvalue_de,decreasing = F,index.return=T)$x
  index_pvalue=sort(pvalue_de,decreasing = F,index.return=T)$ix
  main_index=0
  interaction_index=0
  for(h in 1:p){
    if(index_pvalue[h]>(p-1)/2){
      interaction_index=append(interaction_index,index_pvalue[h]-(p-1)/2)
    }
    if(index_pvalue[h]<=(p-1)/2){
      main_index=append(interaction_index,index_pvalue[h])
    }
  }
  interaction_index=interaction_index[-1]
  main_index=main_index[-1]
  err_var=sighat^2
  pvalue_delist=list(pvalue=order_pvalue,sequence=index_pvalue,
                     main_sequence=main_index,
                     interaction_sequence=interaction_index,
                     error_variance=err_var)
  return(pvalue_delist)
}




n=200
p=500

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
rhocorr=0.9
covmat=rhocorr^{abs(outer(1:p,1:p,'-'))}
gamma_par=c(rep(1,20),rep(0,p-20))
beta_par=c(1,rep(1,30),rep(0,p-30))
signal_index=c(1:20,500+(1:31))
s=51

#simulation starts here!
sim_size=50
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
set.seed(9876)
for (hfg in 1:n){
  y[hfg]=rnorm(1,my[hfg],1)
}

newX=cbind(X,trt_X_intercept)
modelfit=pvalue_de(newX,y)

order_pvalue=modelfit$pvalue

#######################################BIC cut
modelr=glmnet(newX,y)
betamatrix=modelr$beta
zeroone=matrix(as.numeric(betamatrix!=0),nrow(betamatrix),ncol(betamatrix))
action=0
for(k in 1:(ncol(betamatrix)-1)){
  temp=(1:p)[(zeroone[,k+1]-zeroone[,k])!=0]
  if(length(temp)>0){
    action=append(action,temp)
  }
}
action=action[-1]
step=0
for(k in 1:(ncol(betamatrix)-1)){
  temp=(rep(k,p))[(zeroone[,k+1]-zeroone[,k])!=0]
  if(length(temp)>0){
    step=append(step,temp)
  }
}
step=step[-1]
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
  sBIC[k]=s-length(setdiff(signal_index,action))
}
if(BIC[k]>0 & BIC[k]<lenr){
  select_set_BIC=action[1:BIC[k]]
  sBIC[k]=s-length(setdiff(signal_index,select_set_BIC))
}
if(BIC[k]==0){
  sBIC[k]=0
}


###############################lasso cross validation
selected_set=lasso.cv(newX,y)
LCV[k]=length(selected_set)
if(LCV[k]>=lenr){
  sLCV[k]=s-length(setdiff(signal_index,action))
}
if(LCV[k]>0 & LCV[k]<lenr){
  select_set_LCV=action[1:LCV[k]]
  sLCV[k]=s-length(setdiff(signal_index,select_set_LCV))
}
if(LCV[k]==0){
  sLCV[k]=0
}

###################################### Bon cut
BON[k]=sum(((order_pvalue-0.05/(2*p+1))<0)==T)

if(BON[k]>0){
  sBON[k]=s-length(setdiff(signal_index,modelfit$sequence[1:BON[k]]))
}
if(BON[k]==0){
  sBON[k]=0
}

###################################### FDR cut
FDR[k]=sum(((order_pvalue-0.05*(1:(2*p+1))/(2*p+1))<0)==T)

if(FDR[k]>0){
  sFDR[k]=s-length(setdiff(signal_index,modelfit$sequence[1:FDR[k]]))
}
if(FDR[k]==0){
  sFDR[k]=0
}

###################################### New procedure
#bounding sequence
ppp=1-1/sqrt(log(2*p+1))
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
lpm=2*p+1
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
write.table(temp,file="Al11.csv",row.names=F)

mean((BIC-sBIC)/(rep(2*p+1-s,sim_size)),na.rm=TRUE);sd((BIC-sBIC)/(rep(2*p+1-s,sim_size)),na.rm=TRUE)
mean((LCV-sLCV)/(rep(2*p+1-s,sim_size)),na.rm=TRUE);sd((LCV-sLCV)/(rep(2*p+1-s,sim_size)),na.rm=TRUE)
mean((BON-sBON)/(rep(2*p+1-s,sim_size)),na.rm=TRUE);sd((BON-sBON)/(rep(2*p+1-s,sim_size)),na.rm=TRUE)
mean((FDR-sFDR)/(rep(2*p+1-s,sim_size)),na.rm=TRUE);sd((FDR-sFDR)/(rep(2*p+1-s,sim_size)),na.rm=TRUE)
mean((NPE-sNPE)/(rep(2*p+1-s,sim_size)),na.rm=TRUE);sd((NPE-sNPE)/(rep(2*p+1-s,sim_size)),na.rm=TRUE)
