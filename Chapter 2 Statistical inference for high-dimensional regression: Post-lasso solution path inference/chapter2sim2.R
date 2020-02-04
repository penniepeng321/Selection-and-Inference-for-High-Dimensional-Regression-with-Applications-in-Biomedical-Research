#simulation 
rm(list=ls())

setwd("/mnt/home/hpeng2")

library(lars)
library(MASS)
library(hdi)
library(scalreg)

AR=function(n,p,rhoo){
  X=matrix(0,nrow=n,ncol=p)
  X[,1]=matrix(rnorm(n),n,1)
  for(xp in 2:p){
    X[,xp]=rhoo*X[,(xp-1)]+sqrt(1-rhoo^2)*matrix(rnorm(n),n,1)
  }
  return(X)
}

#covariance test p values
pvalue3=function(X,y){
  n=nrow(X)
  p=ncol(X)
  X=X-matrix(colMeans(X),n,p,byrow=T)
  y=y-matrix(mean(y),n,1)
  l1=lars(X,y,type="lar",use.Gram=F)
  #sequential index of added varibles
  LARpath=as.numeric(l1$actions)
  #knots
  knot=l1$lambda
  num_ts=length(knot)-1
  #store test statistic
  covTS=rep(0,num_ts)
  covtscomp=function(i){
    knot[i]*(knot[i]-knot[i+1]) 
  }
  covTS=sapply(1:num_ts,covtscomp)
  qstar=rep(0,num_ts)
  qstarcomp=function(i){
    exp(-sum(covTS[i:num_ts]))
  }
  qstar=sapply(1:num_ts,qstarcomp)
  pvalue=list(model=l1,LAR=LARpath,q_value=qstar,covT=covTS)
  return(pvalue)
}

############################################################
n=200
p=2000
############################################################
m=min(n-1,p)
sim=1000
Vm=rep(0,sim)
set.seed(9871)
for(hj in 1:sim){
  X=matrix(rnorm(n*p,0,1),n,p)
  y=matrix(rnorm(n,0,1),n,1)
  L1=pvalue3(X,y)
  Um=L1$q_value
  Vm[hj]=max((((1:(m-1))/(m-1)-Um)/sqrt(Um*(1-Um)))[1:round(m/2)])
}
Vm=sort(Vm)

############################################################


#simulation setting
s=40
b_s=0.7
b_v=rep(0,p)
signal_index=c(ceiling((1:s)*p/s))
b_v[signal_index]=b_s
############################################################
#simulation starts here!
sim_size=1000
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
path_length=rep(0,sim_size)
num_path_sig=rep(0,sim_size)
m0vec=rep(0,sim_size)
m1vec=rep(0,sim_size)
count=rep(0,sim_size)
countlb=rep(0,sim_size)
countub=rep(0,sim_size)
set.seed(3241)
for(k in 1:sim_size){
  #X=matrix(rnorm(n*p,0,1),n,p,byrow=T)
  X=AR(n,p,0.5)
  ymean=X%*%b_v
  y=rep(0,n)
  for(yin in 1:n){
    y[yin]=rnorm(1,ymean[yin],1)
  }
  Lk=pvalue3(X,y)
  qts1=Lk$q_value
  act1=Lk$LAR
  TS1=Lk$covT
  num_path_sig[k]=s-length(setdiff(signal_index,act1))
  temp=setdiff(signal_index,setdiff(signal_index,act1))
  hi=length(temp)
  positions=rep(0,hi)
  link=rep(0,hi)
  for(pos in 1:hi){
    positions[pos]=match(temp[pos],act1)
  }
  for(pos in 1:hi){
    link[pos]=sum(positions<=pos, na.rm=TRUE)/pos
  }
  if(sum(positions,na.rm=TRUE)>0){
    m0vec[k]=sum(link==1, na.rm=TRUE)
    m1vec[k]=max(positions, na.rm=TRUE)
  }
  else{
    m0vec[k]=sum(link==1, na.rm=TRUE)
    m1vec[k]=min(n-1,p)
  }
  
  #######################################BIC cut
#   m=length(TS1)
#   lpm=length(qts1)
  modelr=glmnet(X,y)
  knotr=modelr$lambda
  lenr=length(knotr)
  BICc=rep(0,lenr+1)
#   for(lp in 1:lpm){
#     temp=lm(y~X[,(act1[1:lp])])
#     BICc[lp]=sum((temp$residuals)^2)+lp*log(p)
#   }
  modelg=scalreg(X,y)
  sigmag=modelg$hsigma
  BICc[1]=sum((y-mean(y))^2)/(sigmag^2)
  intX=cbind(rep(1,n),X)
  for(lp in 2:(lenr+1)){
    temp=matrix(as.numeric(unlist(coef(modelr,s=knotr[lp-1]))),p+1,1)
    BICc[lp]=(1/(sigmag^2))*sum((y-intX%*%temp)^2)+lp*log(p)
  }
  BIC[k]=which.min(BICc)
  if(BIC[k]>=(lpm-1)){
    sBIC[k]=s-length(setdiff(signal_index,act1))
  }
  if(BIC[k]<(lpm-1) & BIC[k]>0){
    sBIC[k]=s-length(setdiff(signal_index,act1[1:BIC[k]]))
  }
  if(BIC[k]==0){
    sBIC[k]=0
  }
  
  
  ###############################lasso cross validation
  selected_set=lasso.cv(X,y)
  LCV[k]=length(selected_set)
  if(LCV[k]>=(lpm-1)){
    sLCV[k]=s-length(setdiff(signal_index,act1))
  }
  if(LCV[k]<(lpm-1) & LCV[k]>0){
    sLCV[k]=s-length(setdiff(signal_index,act1[1:LCV[k]]))
  }
  if(LCV[k]==0){
    sLCV[k]=0
  }
  
  ###################################### Bon cut
  BON[k]=sum(((qts1-0.05/lpm)<0)==T)
  if(BON[k]>=(lpm-1)){
    sBON[k]=s-length(setdiff(signal_index,act1))
  }
  if(BON[k]<(lpm-1) & BON[k]>0){
    sBON[k]=s-length(setdiff(signal_index,act1[1:BON[k]]))
  }
  if(BON[k]==0){
    sBON[k]=0
  }
  
  ###################################### FDR cut
  FDR[k]=sum(((qts1-0.05*(1:lpm)/lpm)<0)==T)
  if(FDR[k]>=(lpm-1)){
    sFDR[k]=s-length(setdiff(signal_index,act1))
  }
  if(FDR[k]<(lpm-1) & FDR[k]>0){
    sFDR[k]=s-length(setdiff(signal_index,act1[1:FDR[k]]))
  }
  if(FDR[k]==0){
    sFDR[k]=0
  }
  
  ###################################### New procedure
  #bounding sequence
  ppp=1-1/sqrt(log(lpm))
  Vm=as.numeric(unlist(Vm))
  cbd=as.numeric(quantile(Vm,max(ppp,0)))
  NPE[k]=max(ceiling(max((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1)/(1-qts1))[1:round(lpm/2)])*m),0)
  if(NPE[k]>=(lpm-1)){
    sNPE[k]=s-length(setdiff(signal_index,act1))
  }
  if(NPE[k]<(lpm-1) & NPE[k]>0){
    sNPE[k]=s-length(setdiff(signal_index,act1[1:NPE[k]]))
  }
  if(NPE[k]==0){
    sNPE[k]=0
  }
  count[k]=as.numeric(((NPE[k]>=m0vec[k])&(NPE[k]<=m1vec[k])))
  countlb[k]=as.numeric(((NPE[k]<=m1vec[k])))
  countub[k]=as.numeric(((NPE[k]>=m0vec[k])))
  path_length[k]=length(TS1)
}
############################################################
c(mean(BIC),sd(BIC),mean(sBIC),sd(sBIC),mean(sBIC/num_path_sig, na.rm=TRUE),sd(sBIC/num_path_sig, na.rm=TRUE),mean((BIC-sBIC)/BIC, na.rm=TRUE),sd((BIC-sBIC)/BIC, na.rm=TRUE))

c(mean(LCV),sd(LCV),mean(sLCV),sd(sLCV),mean(sLCV/num_path_sig, na.rm=TRUE),sd(sLCV/num_path_sig, na.rm=TRUE),mean((LCV-sLCV)/LCV, na.rm=TRUE),sd((LCV-sLCV)/LCV, na.rm=TRUE))

c(mean(BON),sd(BON),mean(sBON),sd(sBON),mean(sBON/num_path_sig, na.rm=TRUE),sd(sBON/num_path_sig, na.rm=TRUE),mean((BON-sBON)/BON, na.rm=TRUE),sd((BON-sBON)/BON, na.rm=TRUE))

c(mean(FDR),sd(FDR),mean(sFDR),sd(sFDR),mean(sFDR/num_path_sig, na.rm=TRUE),sd(sFDR/num_path_sig, na.rm=TRUE),mean((FDR-sFDR)/FDR, na.rm=TRUE),sd((FDR-sFDR)/FDR, na.rm=TRUE))

c(mean(NPE),sd(NPE),mean(sNPE),sd(sNPE),mean(sNPE/num_path_sig, na.rm=TRUE),sd(sNPE/num_path_sig, na.rm=TRUE),mean((NPE-sNPE)/NPE, na.rm=TRUE),sd((NPE-sNPE)/NPE, na.rm=TRUE))

c(mean(num_path_sig),sd(num_path_sig))
c(mean(m0vec, na.rm=TRUE),sd(m0vec, na.rm=TRUE),mean(m1vec, na.rm=TRUE),sd(m1vec, na.rm=TRUE))
c(mean(count, na.rm=TRUE),mean(countlb, na.rm=TRUE),mean(countub, na.rm=TRUE))


temp=matrix(c(BIC,sBIC,LCV,sLCV,BON,sBON,FDR,sFDR,NPE,sNPE,m0vec,m1vec,countlb,countub,count,num_path_sig),nrow=sim_size,ncol=16,byrow=F)
write.table(temp,file="s3.csv",row.names=F)

mean((BIC-sBIC)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE);sd((BIC-sBIC)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE)
mean((LCV-sLCV)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE);sd((LCV-sLCV)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE)
mean((BON-sBON)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE);sd((BON-sBON)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE)
mean((FDR-sFDR)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE);sd((FDR-sFDR)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE)
mean((NPE-sNPE)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE);sd((NPE-sNPE)/(rep(199,sim_size)-num_path_sig),na.rm=TRUE)

############################################################
