#libraries attached

library(lattice) 
library(scalreg)
library(glmnet)
library(glasso)
library(camel)
library(clime)
library(hdi)
library(mvtnorm)
library(Matrix)

#p value func


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
  temp=camel.tiger(X[,1:((p-1)/2)], lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL,
                   method = "slasso", sym = "or", shrink=NULL, prec = 1e-4, mu = 0.01,
                   max.ite = 1e4, standardize = T, correlation = FALSE,
                   perturb = TRUE, verbose = F)
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




#real data work
readin=read.table('/media/huimin/F67428F17428B5F1/meeting notes/stard.csv',
                  sep=',',header=T)
str(readin)
readinmatrix=matrix(c(unlist(readin)),319,308)
sum(is.na(readinmatrix))
n=319
p=305
X=readinmatrix[,4:308]

scaled.dat <- scale(X)

# check that we get mean of 0 and sd of 1
#colMeans(scaled.dat)  # faster version of apply(scaled.dat, 2, mean)
#apply(scaled.dat, 2, sd)


X=scaled.dat

X_intercept=cbind(rep(1,n),X)

trt=readinmatrix[,2]-1
trt_X_intercept=trt*X_intercept

sum(trt==0)/n
sum(trt==1)/n


#original data fit

pro=0.5
newX=cbind(X,(trt-pro)*X_intercept)
y=readinmatrix[,3]
set.seed(1234)
modelfit=pvalue_de(newX,y)
order_pvalue=modelfit$pvalue
order_pvalue

hist(order_pvalue,breaks=20,main='p value histogram',xlab='ordered p value')

bhat1=modelfit$bhat[1]
bhat1
bhat2=modelfit$bhat[2]
bhat2
betahat1=modelfit$betahat[1]
betahat1
betahat2=modelfit$betahat[2]
betahat2
sebhat1=sqrt(modelfit$sigmadiag[1])/sqrt(n)
sebhat1
sebhat2=sqrt(modelfit$sigmadiag[2])/sqrt(n)
sebhat2

int_names=c(colnames(readin)[-c(1,3)])

selected_set=lasso.cv(newX,y)
select_set_LCV=selected_set[selected_set>p]-p
LCV=length(select_set_LCV)
int_names[select_set_LCV]
select_set_LCV
LCV

order_pvalue[match(select_set_LCV,modelfit$sequence)]

BON=sum(((order_pvalue-0.05/(p+1))<0)==T)

BON
int_names[modelfit$sequence[1:BON]]
modelfit$sequence[1:BON]

FDR=sum(((order_pvalue-0.05*(1:(p+1))/(p+1))<0)==T)

FDR
int_names[modelfit$sequence[1:FDR]]
modelfit$sequence[1:FDR]

modelfit$sequence

int_names[modelfit$sequence]

modelfit$bhat[208]
sqrt(modelfit$sigmadiag[208])/sqrt(n)
modelfit$betahat[208]

#sig 0.01 var sel

match('DKPBM',int_names)
EX=readinmatrix[,match('DKPBM',int_names)+2]
EY=readinmatrix[,3]
Etrt=readinmatrix[,2]
Eindicator1=(EX==1)
Eindicator0=(EX==0)
EindicatorBUP=(Etrt==1)
EindicatorSER=(Etrt==2)
mean(EY[Eindicator1&EindicatorBUP]);sd(EY[Eindicator1&EindicatorBUP])
sum(Eindicator1&EindicatorBUP)
mean(EY[Eindicator1&EindicatorSER]);sd(EY[Eindicator1&EindicatorSER])
sum(Eindicator1&EindicatorSER)
mean(EY[Eindicator0&EindicatorBUP]);sd(EY[Eindicator0&EindicatorBUP])
sum(Eindicator0&EindicatorBUP)
mean(EY[Eindicator0&EindicatorSER]);sd(EY[Eindicator0&EindicatorSER])
sum(Eindicator0&EindicatorSER)

#sig 0.05 var sel


(int_names[modelfit$sequence])[order_pvalue<0.05]
order_pvalue[order_pvalue<0.05]

order_pvalue[match('tx',(int_names[modelfit$sequence]))]
match('tx',(int_names[modelfit$sequence]))

#refit model after var sel


refit_X=readinmatrix[,c(match('DKPBM',int_names)+2,
                   match('NVCRD',int_names)+2,
                   match('RENAL',int_names)+2,
                   match('EMSOC',int_names)+2,
                   match('student0',int_names)+2,
                   match('DSMLE',int_names)+2,
                   match('HINTR',int_names)+2,
                   match('student1',int_names)+2,
                   match('ANWOR',int_names)+2)]
refit_Y=readinmatrix[,3]
refit_trt=readinmatrix[,2]
scaled.dat <- scale(refit_X)
refit_X=scaled.dat
refit_X_intercept=cbind(rep(1,n),refit_X)
refit_trt=readinmatrix[,2]-1
refit_trt_X_intercept=refit_trt*refit_X_intercept
pro=0.5
refit_newX=cbind(refit_X,(refit_trt-pro)*refit_X_intercept)
refit_y=readinmatrix[,3]
set.seed(1234)
refit_modelfit=pvalue_de(refit_newX,refit_y)
refit_order_pvalue=refit_modelfit$pvalue
refit_order_pvalue
refit_modelfit$betahat
refit_modelfit$bhat

#lasso refit after var sel


EY=readinmatrix[,3]
Etrt=readinmatrix[,2]
Eindicator1=(refit_X_intercept%*%(refit_modelfit$betahat[ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])>0)
Eindicator0=(refit_X_intercept%*%(refit_modelfit$betahat[ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])<0)
EindicatorBUP=(Etrt==1)
EindicatorSER=(Etrt==2)
mean(EY[Eindicator1&EindicatorBUP]);sd(EY[Eindicator1&EindicatorBUP])
sum(Eindicator1&EindicatorBUP)
mean(EY[Eindicator1&EindicatorSER]);sd(EY[Eindicator1&EindicatorSER])
sum(Eindicator1&EindicatorSER)
mean(EY[Eindicator0&EindicatorBUP]);sd(EY[Eindicator0&EindicatorBUP])
sum(Eindicator0&EindicatorBUP)
mean(EY[Eindicator0&EindicatorSER]);sd(EY[Eindicator0&EindicatorSER])
sum(Eindicator0&EindicatorSER)


#refit_X_intercept%*%(refit_modelfit$betahat[ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])

sum(refit_X_intercept%*%(refit_modelfit$betahat[ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])>0)

sum(refit_X_intercept%*%(refit_modelfit$betahat[ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])<0)

#desparsified lasso refit after var sel


EY=readinmatrix[,3]
Etrt=readinmatrix[,2]
Eindicator1=(refit_X_intercept%*%(refit_modelfit$bhat)>0)
Eindicator0=(refit_X_intercept%*%(refit_modelfit$bhat)<0)
EindicatorBUP=(Etrt==1)
EindicatorSER=(Etrt==2)
mean(EY[Eindicator1&EindicatorBUP]);sd(EY[Eindicator1&EindicatorBUP])
sum(Eindicator1&EindicatorBUP)
mean(EY[Eindicator1&EindicatorSER]);sd(EY[Eindicator1&EindicatorSER])
sum(Eindicator1&EindicatorSER)
mean(EY[Eindicator0&EindicatorBUP]);sd(EY[Eindicator0&EindicatorBUP])
sum(Eindicator0&EindicatorBUP)
mean(EY[Eindicator0&EindicatorSER]);sd(EY[Eindicator0&EindicatorSER])
sum(Eindicator0&EindicatorSER)


sum(refit_X_intercept%*%(refit_modelfit$bhat)>0)

sum(refit_X_intercept%*%(refit_modelfit$bhat)<0)

nrow(X)


#LS refit after var sel
# lasso and desparsified lasso fit both the main and interaction according
# to the definitions in the paper
# LS only fits Y~interaction with intercept
# decision rule is X_intercept%*%coef_est

LS_coef=lm(refit_y~refit_newX[,ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)])$coefficients
LS_coef=as.numeric(LS_coef)
temp=refit_X_intercept%*%LS_coef[-1]

LS_coef[-1]

EY=readinmatrix[,3]
Etrt=readinmatrix[,2]
Eindicator1=(temp>0)
Eindicator0=(temp<0)
EindicatorBUP=(Etrt==1)
EindicatorSER=(Etrt==2)
mean(EY[Eindicator1&EindicatorBUP]);sd(EY[Eindicator1&EindicatorBUP])
sum(Eindicator1&EindicatorBUP)
mean(EY[Eindicator1&EindicatorSER]);sd(EY[Eindicator1&EindicatorSER])
sum(Eindicator1&EindicatorSER)
mean(EY[Eindicator0&EindicatorBUP]);sd(EY[Eindicator0&EindicatorBUP])
sum(Eindicator0&EindicatorBUP)
mean(EY[Eindicator0&EindicatorSER]);sd(EY[Eindicator0&EindicatorSER])
sum(Eindicator0&EindicatorSER)

sum(temp>0)

sum(temp<0)


#cross-validation with half half split and 100 repetitions

refit_newX_AX=refit_newX[,ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)]
refit_Etrt=readinmatrix[,2]
receive_BUP=(refit_Etrt==1)
receive_SER=(refit_Etrt==2)



cv_reps_num=200

All_SER_mean=rep(NA,cv_reps_num)
All_SER_sd=rep(NA,cv_reps_num)
All_SER_count=rep(NA,cv_reps_num)

All_BUP_mean=rep(NA,cv_reps_num)
All_BUP_sd=rep(NA,cv_reps_num)
All_BUP_count=rep(NA,cv_reps_num)

Opt_Regime_mean=rep(NA,cv_reps_num)
Opt_Regime_sd=rep(NA,cv_reps_num)
Opt_Regime_count=rep(NA,cv_reps_num)

sample_size_org=nrow(newX)

set.seed(1234)
for(cv_rep_index in 1:cv_reps_num){
  
  rand_draw_index_train=sample(1:sample_size_org,200,replace=F)
  rand_draw_index_test=(1:sample_size_org)[-rand_draw_index_train]
  train_newX=refit_newX_AX[rand_draw_index_train,]
  test_newX=refit_newX_AX[rand_draw_index_test,]
  train_y=refit_y[rand_draw_index_train]
  test_y=refit_y[rand_draw_index_test]
  LS_coef=lm(train_y~refit_newX[rand_draw_index_train,])$coefficients
  LS_coef=as.numeric(LS_coef)
  temp=refit_X_intercept%*%LS_coef[(ncol(refit_X_intercept)+1):(2*ncol(refit_X_intercept))]
  decision_SER_test=(temp<0)[rand_draw_index_test]
  decision_BUP_test=(temp>0)[rand_draw_index_test]
  receive_SER_test=receive_SER[rand_draw_index_test]
  receive_BUP_test=receive_BUP[rand_draw_index_test]
  
  All_SER_index=receive_SER_test
  All_SER_mean[cv_rep_index]=mean(test_y[All_SER_index],na.rm=T)
  All_SER_sd[cv_rep_index]=sd(test_y[All_SER_index],na.rm=T)
  All_SER_count[cv_rep_index]=sum(All_SER_index,na.rm=T)
  
  All_BUP_index=receive_BUP_test
  All_BUP_mean[cv_rep_index]=mean(test_y[All_BUP_index],na.rm=T)
  All_BUP_sd[cv_rep_index]=sd(test_y[All_BUP_index],na.rm=T)
  All_BUP_count[cv_rep_index]=sum(All_BUP_index,na.rm=T)
  
  Opt_Regime_index=((decision_SER_test*receive_SER_test)+
    (decision_BUP_test*receive_BUP_test)==1)
  Opt_Regime_mean[cv_rep_index]=mean(test_y[Opt_Regime_index],na.rm=T)
  Opt_Regime_sd[cv_rep_index]=sd(test_y[Opt_Regime_index],na.rm=T)
  Opt_Regime_count[cv_rep_index]=sum(Opt_Regime_index,na.rm=T)
}

boxplot(list(all.SER=All_SER_mean,all.BUP=All_BUP_mean,decision=Opt_Regime_mean),
        horizontal=T, xlab='mean outcome', ylab='treatment plan')



# ####CV using bhat
# 
# refit_newX_AX=refit_newX[,ncol(refit_X_intercept):(2*ncol(refit_X_intercept)-1)]
# refit_Etrt=readinmatrix[,2]
# receive_BUP=(refit_Etrt==1)
# receive_SER=(refit_Etrt==2)
# 
# 
# 
# cv_reps_num=200
# 
# All_SER_mean=rep(NA,cv_reps_num)
# All_SER_sd=rep(NA,cv_reps_num)
# All_SER_count=rep(NA,cv_reps_num)
# 
# All_BUP_mean=rep(NA,cv_reps_num)
# All_BUP_sd=rep(NA,cv_reps_num)
# All_BUP_count=rep(NA,cv_reps_num)
# 
# Opt_Regime_mean=rep(NA,cv_reps_num)
# Opt_Regime_sd=rep(NA,cv_reps_num)
# Opt_Regime_count=rep(NA,cv_reps_num)
# 
# sample_size_org=nrow(newX)
# 
# set.seed(1234)
# for(cv_rep_index in 1:cv_reps_num){
#   
#   rand_draw_index_train=sample(1:sample_size_org,200,replace=F)
#   rand_draw_index_test=(1:sample_size_org)[-rand_draw_index_train]
#   train_newX=refit_newX_AX[rand_draw_index_train,]
#   test_newX=refit_newX_AX[rand_draw_index_test,]
#   train_y=refit_y[rand_draw_index_train]
#   test_y=refit_y[rand_draw_index_test]
#   
#   refit_modelfit_train=pvalue_de(refit_newX[rand_draw_index_train,],
#                                  refit_y[rand_draw_index_train])
# 
#   decision_SER_test=(refit_X_intercept%*%
#                  (refit_modelfit_train$bhat)<0)[rand_draw_index_test]
#   decision_BUP_test=(refit_X_intercept%*%
#                  (refit_modelfit_train$bhat)>0)[rand_draw_index_test]
#   
#   receive_SER_test=receive_SER[rand_draw_index_test]
#   receive_BUP_test=receive_BUP[rand_draw_index_test]
#   
#   All_SER_index=receive_SER_test
#   All_SER_mean[cv_rep_index]=mean(test_y[All_SER_index],na.rm=T)
#   All_SER_sd[cv_rep_index]=sd(test_y[All_SER_index],na.rm=T)
#   All_SER_count[cv_rep_index]=sum(All_SER_index,na.rm=T)
#   
#   All_BUP_index=receive_BUP_test
#   All_BUP_mean[cv_rep_index]=mean(test_y[All_BUP_index],na.rm=T)
#   All_BUP_sd[cv_rep_index]=sd(test_y[All_BUP_index],na.rm=T)
#   All_BUP_count[cv_rep_index]=sum(All_BUP_index,na.rm=T)
#   
#   Opt_Regime_index=((decision_SER_test*receive_SER_test)+
#                       (decision_BUP_test*receive_BUP_test)==1)
#   Opt_Regime_mean[cv_rep_index]=mean(test_y[Opt_Regime_index],na.rm=T)
#   Opt_Regime_sd[cv_rep_index]=sd(test_y[Opt_Regime_index],na.rm=T)
#   Opt_Regime_count[cv_rep_index]=sum(Opt_Regime_index,na.rm=T)
# }
# 
# boxplot(list(all.SER=All_SER_mean,all.BUP=All_BUP_mean,decision=Opt_Regime_mean),
#         horizontal=T, xlab='mean outcome', ylab='treatment plan')
