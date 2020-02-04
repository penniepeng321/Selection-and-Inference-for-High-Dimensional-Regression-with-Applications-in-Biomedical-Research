rm(list = ls())

library(MASS)
library(hdi) 

readin=read.table('stard.csv',sep=',',header=T)
str(readin)
colnames(readin)

readinmatrix=matrix(c(unlist(readin)),319,308)
sum(is.na(readinmatrix))
n=319
p=305
X=readinmatrix[,4:308]

refit_X=scale(X)
refit_X_intercept=cbind(rep(1,n),refit_X)
refit_trt=readinmatrix[,2]-1

refit_Y=-readinmatrix[,3]

receive_BUP=(refit_trt==0)
receive_SER=(refit_trt==1)

prohat=sum(receive_SER)/n
#prohat=0.5
refit_newX=cbind(refit_X,(refit_trt-prohat)*refit_X_intercept)

debias_stard = lasso.proj(refit_newX, refit_Y, betainit = "scaled lasso")

selected_ind=which(debias_stard$pval[(p+1): (2*p+1)]<0.05) #index based on refit_X_intercept

selected_coef = debias_stard$bhat[p+selected_ind]

selected_se = debias_stard$se[p+selected_ind]

selected_pval = debias_stard$pval[p+selected_ind]

selected_sort_p = sort(selected_pval, index=T)

selected_sort_ind = selected_ind[selected_sort_p$ix]
readin[1, selected_sort_ind-1+3] #minus intercept plus first 3 colums of readin
colnames(readin[,selected_sort_ind-1+3]) #minus intercept plus first 3 colums of readin
selected_sort_p$x
selected_coef[selected_sort_p$ix]
selected_se[selected_sort_p$ix]

#h1=colnames(readin[,selected_sort_ind-1+3])
h1=1:14
h2=selected_coef[selected_sort_p$ix]
h3=selected_coef[selected_sort_p$ix]+qnorm(0.975)*selected_se[selected_sort_p$ix]
h4=selected_coef[selected_sort_p$ix]-qnorm(0.975)*selected_se[selected_sort_p$ix]

require(plotrix)
plotCI(h1, 
       h2, 
       ui=h3, 
       li=h4,xlab='index of selected interaction',
       ylab='confidence interval of interaction coefficient')

### select 14 covariates when set pval threshold at 0.05
#1.qccur_r_rate(QCCURRATE): QIDS-C score changing rates at Level (chengchun)
#2.URNONE: no symptoms in patients' urination category (chengchun)
#3.NVTRM: CNS: Tremors (chengchun)
#4.IMPWR: indicting whether patients thought they have special powers (chengchun)
#5.HWL: HRS Weight loss
#7.DSMTD: Recurrent thoughts of death, recurrent suicidal ideation, or suicide attempt
#9.HMNIN: HRS Middle insomnia
#14.EMSTU: Did you worry a lot that you might do something to make people think that you were stupid or foolish?


### Bootstrap to compare distributions of response using different trt regime

#### refit OLS and get coefficients
#nonzero_ind = which(debias_stard$pval<0.05) #index based on newX
#nonzero_ind = which(debias_stard$pval<0.01) #index based on newX
#refit = lm(refit_Y~refit_newX[,nonzero_ind])
#refit_coef = refit$coefficients


bootstrap_num=1000

SER_mean = rep(NA, bootstrap_num)
SER_count = rep(NA, bootstrap_num)


BUP_mean = rep(NA, bootstrap_num)
BUP_count = rep(NA, bootstrap_num)


Opt_mean = rep(NA, bootstrap_num)
Opt_count = rep(NA, bootstrap_num)


set.seed(1234)
for(t in 1:bootstrap_num){
  
  draw_ind=sample(1:n,n,replace=T)
  
  temp=refit_X_intercept[draw_ind,selected_ind]%*%selected_coef
  #temp=refit_X_intercept[draw_ind,selected_ind]%*%refit_coef[28:41]
  #temp=refit_X_intercept[draw_ind,selected_ind]%*%refit_coef[14:16]
  
  refit_Y_draw = refit_Y[draw_ind]
  
  decision_SER=(temp>0)
  decision_BUP=(temp<0)
  
  SER_ind=receive_SER[draw_ind]
  SER_mean[t]=mean(refit_Y_draw[SER_ind])
  SER_count[t]=sum(SER_ind,na.rm=T)

  
  BUP_ind=receive_BUP[draw_ind]
  BUP_mean[t]=mean(refit_Y_draw[BUP_ind])
  BUP_count[t]=sum(BUP_ind,na.rm=T)

  
  Opt_ind=((decision_SER*SER_ind)+(decision_BUP*BUP_ind)==1)
  Opt_mean[t]=mean(refit_Y_draw[Opt_ind])
  Opt_count[t]=sum(Opt_ind,na.rm=T)

}


boxplot(list(all.SER=SER_mean,all.BUP=BUP_mean,optimal=Opt_mean),
        horizontal=T, xlab='mean outcome', ylab='treatment regime')






cv_num=1000

SER_mean = rep(NA, cv_num)
SER_count = rep(NA, cv_num)


BUP_mean = rep(NA, cv_num)
BUP_count = rep(NA, cv_num)


Opt_mean = rep(NA, cv_num)
Opt_count = rep(NA, cv_num)


set.seed(1234)
for(t in 1:cv_num){
  
  train_ind=sample(1:n,ceiling(n/2),replace=F)
  test_ind=setdiff(1:n,train_ind)
  
  
  
  refit_Y_train = refit_Y[train_ind]
  refit_X_train=refit_newX[train_ind,selected_sort_ind+p]
  test_dec=refit_newX[test_ind,selected_sort_ind-1]%*%(as.numeric(lm(refit_Y_train~refit_X_train)$coef)[-1])
  
  refit_Y_test = refit_Y[test_ind]
  
  decision_SER=(test_dec>0)
  decision_BUP=(test_dec<0)
  
  SER_ind=receive_SER[test_ind]
  SER_mean[t]=mean(refit_Y_test[SER_ind])
  SER_count[t]=sum(SER_ind,na.rm=T)
  
  
  BUP_ind=receive_BUP[test_ind]
  BUP_mean[t]=mean(refit_Y_test[BUP_ind])
  BUP_count[t]=sum(BUP_ind,na.rm=T)
  
  
  Opt_ind=((decision_SER*SER_ind)+(decision_BUP*BUP_ind)==1)
  Opt_mean[t]=mean(refit_Y_test[Opt_ind])
  Opt_count[t]=sum(Opt_ind,na.rm=T)
  
}


boxplot(list(all.SER=SER_mean,all.BUP=BUP_mean,optimal=Opt_mean),
        horizontal=T, xlab='mean outcome', ylab='treatment regime')






