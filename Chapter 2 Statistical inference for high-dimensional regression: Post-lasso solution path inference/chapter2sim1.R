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

adaplasso_inter=function(newX,y){
  n=nrow(newX)
  p=(ncol(newX)-1)/2
  initial_value=lm(y~newX)$coefficients
  gamma_init_value=initial_value[1:p]
  beta_init_value=initial_value[(p+1):(2*p+1)]
  weights=abs(beta_init_value)^(-1)
  compute_y=y-newX[,1:p]%*%gamma_init_value
  compute_X=newX[,(p+1):(2*p+1)]
  pls_model=glmnet(compute_X,compute_y,
                   family='gaussian',
                   penalty.factor=weights)
  beta_lambda_est=pls_model$beta
  #ncol(beta_lambda_est)
  lambda_seq=pls_model$lambda
  length_lambda_seq=length(lambda_seq)
  BIC_crit=rep(0,length_lambda_seq)
  for(lambda_index in 1:length_lambda_seq){
    temp1=c(beta_lambda_est[,lambda_index])
    temp2=sum(temp1!=0)
    BIC_crit[lambda_index]=
      sum((compute_y-compute_X%*%temp1)^2)/
      sum((compute_y-compute_X%*%beta_init_value)^2)+
      temp2*log(n)/n
  }
  pick_BIC_index=which.min(BIC_crit)
  BIC_minimum=BIC_crit[pick_BIC_index]
  pick_lambda=lambda_seq[pick_BIC_index]
  pick_beta=c(beta_lambda_est[,pick_BIC_index])
  select_interaction_terms=(1:(p+1))[pick_beta!=0]
  interaction_selection_results=list(
    lam=pick_lambda,
    selected_set=select_interaction_terms,
    coef=pick_beta
  )
  return(interaction_selection_results)
}



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
  temp=camel.tiger(X[,1:((p-1)/2)],method = "slasso", standardize = T)
  temp1=matrix(as.numeric(unlist(temp$icov[1])),((p-1)/2),((p-1)/2))
  precisionmatrix=4*as.matrix(bdiag(1,temp1))
  
  aspvarmatrix=precisionmatrix%*%samplecov%*%t(precisionmatrix)
  sigmadiag=diag(aspvarmatrix)
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
                     error_variance=err_var,betahat=betahat,bhat=bhat,
                     sigmadiag=sigmadiag)
  return(pvalue_delist)
}




n=200
p=100

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

s=5
inten=4
gamma_par=c(rep(1,4),rep(0,p-4))
beta_par=c(rep(inten,s),rep(0,p-s+1))
signal_index=c(1:s)

#simulation starts here!
sim_size=1000
ALS=rep(0,sim_size)
sALS=rep(0,sim_size)
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
cov_rate_noise=rep(0,sim_size)
cov_rate_true=rep(0,sim_size)
lenci_noise=rep(0,sim_size)
lenci_true=rep(0,sim_size)
power=rep(0,sim_size)
typeIerr=rep(0,sim_size)
typeIIerr=rep(0,sim_size)
bhat1=rep(0,sim_size)
bhat2=rep(0,sim_size)
betahat1=rep(0,sim_size)
betahat2=rep(0,sim_size)
sebhat1=rep(0,sim_size)
sebhat2=rep(0,sim_size)

set.seed(3241)
for(k in 1:sim_size){
  
  
  X = rmvnorm(n, mean = rep(0, nrow(covmat)), sigma = covmat,method=c("chol"))
  X_intercept=cbind(rep(1,n),X)
  
  trt=rbinom(n,size=1,prob=pro)
  trt_X_intercept=trt*X_intercept
  
  my=1+0.05*exp(X%*%gamma_par)+trt_X_intercept%*%beta_par
  y=rep(NA,n)
  for (hfg in 1:n){
    y[hfg]=rnorm(1,my[hfg],1)
  }
  
  newX=cbind(X,(trt-pro)*X_intercept)
  modelfit=pvalue_de(newX,y)
  
  order_pvalue=modelfit$pvalue
  
  cov_rate_true[k]=(1/(s+1))*sum((modelfit$upperci[1:(s+1)]>=inten)&
                                   (modelfit$lowerci[1:(s+1)]<=inten),na.rm=T)
  cov_rate_noise[k]=(1/(p-s))*sum((modelfit$upperci[(s+2):(p+1)]>=0)&
                                    (modelfit$lowerci[(s+2):(p+1)]<=0),na.rm=T)
  lenci_true[k]=(1/(s+1))*sum(modelfit$lengthci[1:(s+1)],na.rm=T)
  lenci_noise[k]=(1/(p-s))*sum(modelfit$lengthci[(s+2):(p+1)],na.rm=T)
  
  power[k]=1-(1/(p))*sum((modelfit$upperci>=0)&
                           (modelfit$lowerci<=0),na.rm=T)
  typeIerr[k]=1-(1/(p-s))*sum((modelfit$upperci[(s+2):(p+1)]>=0)&
                                (modelfit$lowerci[(s+2):(p+1)]<=0),na.rm=T)
  typeIIerr[k]=(1/(s+1))*sum((modelfit$upperci[1:(s+1)]>=0)&
                               (modelfit$lowerci[1:(s+1)]<=0),na.rm=T)
  
  bhat1[k]=modelfit$bhat[1]
  bhat2[k]=modelfit$bhat[2]
  betahat1[k]=modelfit$betahat[1]
  betahat2[k]=modelfit$betahat[2]
  sebhat1[k]=sqrt(modelfit$sigmadiag[1])/sqrt(n)
  sebhat2[k]=sqrt(modelfit$sigmadiag[2])/sqrt(n)
  
  
  ####################adaptive lasso interaction inference ALS
  adap_model=adaplasso_inter(newX,y)
  ALS_lambda=adap_model$lam
  ALS_coef=adap_model$coef
  ALS_set=adap_model$selected_set
  
  ALS[k]=length(ALS_set)
  if(ALS[k]>0){
    sALS[k]=s-length(setdiff(signal_index,ALS_set))
  }
  if(ALS[k]==0){
    sALS[k]=0
  }
  
  
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

c(mean(cov_rate_true,na.rm=T),sd(cov_rate_true,na.rm=T),
  mean(cov_rate_noise,na.rm=T),sd(cov_rate_noise,na.rm=T),
  mean(lenci_true,na.rm=T),sd(lenci_true,na.rm=T),
  mean(lenci_noise,na.rm=T),sd(lenci_noise,na.rm=T))

c(mean(power,na.rm=T),sd(power,na.rm=T),
  mean(typeIerr,na.rm=T),sd(typeIerr,na.rm=T),
  mean(typeIIerr,na.rm=T),sd(typeIIerr,na.rm=T))

c(mean(bhat1,na.rm=T),sd(bhat1,na.rm=T),mean((bhat1-inten)^2,na.rm=T),
  mean(sebhat1,na.rm=T))
c(mean(bhat2,na.rm=T),sd(bhat2,na.rm=T),mean((bhat2-inten)^2,na.rm=T),
  mean(sebhat2,na.rm=T))
c(mean(betahat1,na.rm=T),sd(betahat1,na.rm=T),mean((betahat1-inten)^2,na.rm=T))
c(mean(betahat2,na.rm=T),sd(betahat2,na.rm=T),mean((betahat2-inten)^2,na.rm=T))

c(mean(ALS,na.rm=T),sd(ALS,na.rm=T),mean(sALS,na.rm=T),sd(sALS,na.rm=T),mean(sALS/s, na.rm=T),sd(sALS/s, na.rm=T),mean((ALS-sALS)/ALS, na.rm=T),sd((ALS-sALS)/ALS, na.rm=T),mean((ALS-sALS)/(p+1-s), na.rm=T),sd((ALS-sALS)/(p+1-s), na.rm=T))

c(mean(BON,na.rm=T),sd(BON,na.rm=T),mean(sBON,na.rm=T),sd(sBON,na.rm=T),mean(sBON/s, na.rm=T),sd(sBON/s, na.rm=T),mean((BON-sBON)/BON, na.rm=T),sd((BON-sBON)/BON, na.rm=T),mean((BON-sBON)/(p+1-s), na.rm=T),sd((BON-sBON)/(p+1-s), na.rm=T))

c(mean(FDR,na.rm=T),sd(FDR,na.rm=T),mean(sFDR,na.rm=T),sd(sFDR,na.rm=T),mean(sFDR/s, na.rm=T),sd(sFDR/s, na.rm=T),mean((FDR-sFDR)/FDR, na.rm=T),sd((FDR-sFDR)/FDR, na.rm=T),mean((FDR-sFDR)/(p+1-s), na.rm=T),sd((FDR-sFDR)/(p+1-s), na.rm=T))

c(mean((ALS-sALS)/ALS,na.rm=T),mean((s-sALS)/(p+1-ALS),na.rm=T))

c(mean((BON-sBON)/BON,na.rm=T),mean((s-sBON)/(p+1-BON),na.rm=T))

c(mean((FDR-sFDR)/FDR,na.rm=T),mean((s-sFDR)/(p+1-FDR),na.rm=T))



c(mean(BIC,na.rm=T),sd(BIC,na.rm=T),mean(sBIC,na.rm=T),sd(sBIC,na.rm=T),mean(sBIC/s, na.rm=T),sd(sBIC/s, na.rm=T),mean((BIC-sBIC)/BIC, na.rm=T),sd((BIC-sBIC)/BIC, na.rm=T),mean((BIC-sBIC)/(p+1-s), na.rm=T),sd((BIC-sBIC)/(p+1-s), na.rm=T))

c(mean(LCV,na.rm=T),sd(LCV,na.rm=T),mean(sLCV,na.rm=T),sd(sLCV,na.rm=T),mean(sLCV/s, na.rm=T),sd(sLCV/s, na.rm=T),mean((LCV-sLCV)/LCV, na.rm=T),sd((LCV-sLCV)/LCV, na.rm=T),mean((LCV-sLCV)/(p+1-s), na.rm=T),sd((LCV-sLCV)/(p+1-s), na.rm=T))

c(mean(NPE,na.rm=T),sd(NPE,na.rm=T),mean(sNPE,na.rm=T),sd(sNPE,na.rm=T),mean(sNPE/s, na.rm=T),sd(sNPE/s, na.rm=T),mean((NPE-sNPE)/NPE, na.rm=T),sd((NPE-sNPE)/NPE, na.rm=T),mean((NPE-sNPE)/(p+1-s), na.rm=T),sd((NPE-sNPE)/(p+1-s), na.rm=T))


############################################################
temp=matrix(c(ALS,sALS,BIC,sBIC,LCV,sLCV,BON,sBON,FDR,sFDR,NPE,sNPE,
              cov_rate_true,cov_rate_noise,lenci_true,lenci_noise,power,typeIerr,typeIIerr),
            nrow=sim_size,ncol=19,byrow=F)
write.table(temp,file="vs6_2.csv",row.names=F)

modelfit$pvalue
modelfit$sequence