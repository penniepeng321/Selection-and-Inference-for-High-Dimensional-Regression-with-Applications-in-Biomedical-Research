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
refit_X=scaled.dat
refit_X_intercept=cbind(rep(1,n),refit_X)
refit_trt=readinmatrix[,2]-1

pro=0.5

refit_y=readinmatrix[,3]

receive_BUP=(refit_trt==0)
receive_SER=(refit_trt==1)

#prohat=sum(receive_SER)/n
prohat=0.5
refit_newX=cbind(refit_X,(refit_trt-prohat)*refit_X_intercept)


#set.seed(7654)
model_fit=pvalue_de(refit_newX,refit_y)

selected_model=model_fit$sequence[(model_fit$pvalue<=0.05)]


ls_fit=lm(refit_y~refit_newX[,c(selected_model,p+1,p+selected_model)])
ls_fit_coef=ls_fit$coefficients
ls_fit_coef=as.numeric(ls_fit_coef)[(length(selected_model)+2):
                                      (length(selected_model)*2+2)]

ls_fit_coef

cv_reps_num=1000

All_SER_mean=rep(NA,cv_reps_num)
All_SER_sd=rep(NA,cv_reps_num)
All_SER_count=rep(NA,cv_reps_num)

All_BUP_mean=rep(NA,cv_reps_num)
All_BUP_sd=rep(NA,cv_reps_num)
All_BUP_count=rep(NA,cv_reps_num)

Opt_Regime_mean=rep(NA,cv_reps_num)
Opt_Regime_sd=rep(NA,cv_reps_num)
Opt_Regime_count=rep(NA,cv_reps_num)

sample_size_org=nrow(X)

set.seed(1234)
for(cv_rep_index in 1:cv_reps_num){
  
  rand_draw_index_boot=sample(1:sample_size_org,n,replace=T)
  
  ls_fit=lm(refit_y[rand_draw_index_boot]~
              refit_newX[rand_draw_index_boot,c(selected_model,p+1,p+selected_model)])
  ls_fit_coef=ls_fit$coefficients
  ls_fit_coef=as.numeric(ls_fit_coef)[(length(selected_model)+2):
                                        (length(selected_model)*2+2)]
  
  temp=refit_X_intercept[rand_draw_index_boot,c(1,1+selected_model)]%*%ls_fit_coef
  
  refit_y_boot=refit_y[rand_draw_index_boot]
  
  decision_SER=(temp>0)
  decision_BUP=(temp<0)
  
  All_SER_index=receive_SER[rand_draw_index_boot]
  
  All_SER_mean[cv_rep_index]=mean(refit_y_boot[All_SER_index],na.rm=T)
  All_SER_sd[cv_rep_index]=sd(refit_y_boot[All_SER_index],na.rm=T)
  All_SER_count[cv_rep_index]=sum(All_SER_index,na.rm=T)
  
  All_BUP_index=receive_BUP[rand_draw_index_boot]
  
  All_BUP_mean[cv_rep_index]=mean(refit_y_boot[All_BUP_index],na.rm=T)
  All_BUP_sd[cv_rep_index]=sd(refit_y_boot[All_BUP_index],na.rm=T)
  All_BUP_count[cv_rep_index]=sum(All_BUP_index,na.rm=T)
  
  Opt_Regime_index=((decision_SER*All_SER_index)+
                      (decision_BUP*All_BUP_index)==1)
  Opt_Regime_mean[cv_rep_index]=mean(refit_y_boot[Opt_Regime_index],na.rm=T)
  Opt_Regime_sd[cv_rep_index]=sd(refit_y_boot[Opt_Regime_index],na.rm=T)
  Opt_Regime_count[cv_rep_index]=sum(Opt_Regime_index,na.rm=T)
  
  
}

boxplot(list(all.SER=All_SER_mean,all.BUP=All_BUP_mean,decision=Opt_Regime_mean),
        horizontal=T, xlab='mean outcome', ylab='treatment plan')

quantile(Opt_Regime_mean-All_SER_mean,c(0.025,0.975),na.rm=T)
quantile(Opt_Regime_mean-All_BUP_mean,c(0.025,0.975),na.rm=T)

quantile(Opt_Regime_mean-All_SER_mean,c(0.05,0.95),na.rm=T)
quantile(Opt_Regime_mean-All_BUP_mean,c(0.05,0.95),na.rm=T)

