#simulation 
rm(list=ls())

#setwd("/mnt/home/hpeng2")

library(covTest)
library(lars)
library(MASS)
library(hdi)





CEU_raw <- read.table("/home/huimin/Documents/CEU.raw", header = T, quote="\"")
CEU_name <- CEU_raw$IID
CEU_data <- CEU_raw[, -(1 : 6)]
CEU_filter <- CEU_data[, is.na(colSums(CEU_raw[, -(1 : 6)])) == F]

CEU_filter <- CEU_filter[, colSums(CEU_filter) >= 0.05 * 2 * nrow(CEU_filter)]


true_included_CEU=c("rs2831459_C","rs7278456_G",'rs965951_A',"rs2832253_T",
                    "rs2832332_G","rs13046799_A")
filter_index=match(true_included_CEU[1],colnames(CEU_filter))
for( hjk in 2:length(true_included_CEU)){
  filter_index=append(filter_index,match(true_included_CEU[hjk],colnames(CEU_filter)))
}
temp=sample(setdiff(colnames(CEU_filter),true_included_CEU),
            3000-length(true_included_CEU))
for(hjk in 1:length(temp)){
  filter_index=append(filter_index,match(temp[hjk],colnames(CEU_filter)))
}
length(filter_index)

CEU_filter <- CEU_filter[filter_index]


namelist <- c("NA06985","NA06993","NA07022","NA07034","NA07055","NA07056","NA07345","NA07357",
              "NA11829","NA11830","NA11831","NA11832","NA11839","NA11840","NA11881","NA11882",
              "NA11992","NA11993","NA11994","NA11995","NA12003","NA12004","NA12056","NA12145",
              "NA12146","NA12154","NA12155","NA12156","NA12236","NA12248","NA12249","NA12264",
              "NA12716","NA12717","NA12750","NA12751","NA12760","NA12761","NA12762","NA12763",
              "NA12812","NA12813","NA12814","NA12815","NA12872","NA12873","NA12874","NA12875",
              "NA12891","NA12892","NA12044","NA12144","NA06994","NA12006","NA12005","NA12057",
              "NA12239","NA12043","NA07000","NA12234")

## Gene ID GI_6005726
y <- c(13.2702080357681, 13.6910317319355, 13.4860788127517, 13.4193482302501, 13.7434130352314,
       13.3479201413790, 13.8043261933933, 13.4317826888239, 13.5310309677652, 13.5932472503920,
       13.5820508920561, 13.7114512904573, 13.7570567579120, 13.8270536759757, 13.607526756319,
       13.7402904908830, 13.9716885519597, 13.4148895323929, 13.8383053398298, 13.6052340534644,
       13.8237916703524, 14.0116833558187, 13.3458216017790, 13.715732383364, 13.5604478792499,
       13.7312797780573, 13.4655309466930, 13.8271088312680, 13.3701661495819, 13.7026521697938,
       13.4329169274760, 13.3046986150768, 13.4171695765824, 13.7927435628403, 13.8265899460629,
       13.5202904154497, 13.8190593626588, 13.4016998181281, 13.3981288335911, 13.3830642376141,
       13.4815716443848, 13.6502264734442, 13.6178822577041, 13.2765829043134, 13.930695227279,
       13.3173366087084, 13.6450971315783, 13.4570184856483, 13.7443947660240, 13.5504467280272,
       13.4179725582229, 13.8810276213521, 13.5535749446040, 14.2258928216667, 13.5972258523391,
       13.7082415009722, 14.0046410157686, 13.7533311543939, 13.3306562912772, 13.5721632074088)

y_order <- y[order(namelist)]




############################################################

n=length(y_order)
p=ncol(CEU_filter)
n;p

X=matrix(as.numeric(unlist(CEU_filter)),n,p,byrow=F)
dim(X)

y=as.numeric(unlist(y_order))
length(y)

n=nrow(X)
p=ncol(X)
X=X-matrix(colMeans(X),n,p,byrow=T)
Xcol=sqrt(colSums(X^2))
for(h in 1:p){
  X[,h]=X[,h]/Xcol[h]
}
y=y-rep(mean(y),n)
for(h in 1:n){
  y[h]=y[h]/sqrt(sum(y^2))
}
obj1=lars(X,y,type="lasso",use.Gram=F)

res1=covTest(obj1,X,y,sigma.est=1)$results
res1
act1=res1[,1]
act1
TS1=res1[,2]
TS1
#drop NAs
chos=!is.na(TS1)
TS1=TS1[chos]
TS1
act1=act1[chos]
act1

qts1=rep(0,length(TS1))
for(ql in 1:length(TS1)){  
  qts1[ql]=exp(-sum(TS1[ql:length(TS1)]))
}

qts1


lpm=length(qts1)

#bounding sequence
m=lpm
sim=10000
Vm=rep(0,sim)
set.seed(1234)
tmep=matrix(rnorm(m*sim,0,1),m,sim)
for(hj in 1:sim){
  ret=c(tmep[,hj])
  ret=abs(ret)
  ret=sort(ret,decreasing=T)
  TT=rep(0,m-1)
  for(j in 1:(m-1)){
    TT[j]=ret[j]*(ret[j]-ret[j+1])
  }
  qts=rep(0,m-1)
  
  for(ql in 1:(m-1)){  
    qts[ql]=exp(-sum(TT[ql:(m-1)]))
  }
  Um=qts
  Vm[hj]=max((((1:(m-1))/(m-1)-Um)/sqrt(Um*(1-Um)))[1:round(m/2)])
}
Vm=sort(Vm)


ppp=1-1/sqrt((log(lpm)))
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(CEU_filter)[abs(act1)]


(colnames(CEU_filter)[abs(act1)])[1:NPE]



lpm=length(qts1)
BICc=rep(0,lpm)
for(lp in 1:lpm){
  temp=lm(y~X[,abs(act1[1:lp])])
  BICc[lp]=sum((temp$residuals)^2)+lp*log(p)
}
BIC=which.min(BICc)
BIC
BICc

###############################lasso cross validation
selected_set=lasso.cv(X,y)
LCV=length(selected_set)
LCV

###################################### Bon cut
BON=sum(((qts1-0.05/lpm)<0)==T)
BON

###################################### FDR cut
FDR=sum(((qts1-0.05*(1:lpm)/lpm)<0)==T)
FDR

(colnames(CEU_filter)[abs(act1)])[1:NPE]
(colnames(CEU_filter)[abs(act1)])[1:LCV]
(colnames(CEU_filter)[abs(act1)])[1:BIC]
(colnames(CEU_filter)[abs(act1)])[1:BON]
(colnames(CEU_filter)[abs(act1)])[1:FDR]


match("rs2831459_A",colnames(CEU_filter))
match("rs2831459_T",colnames(CEU_filter))
match("rs2831459_C",colnames(CEU_filter))
match("rs2831459_G",colnames(CEU_filter))
match("rs7277536_A",colnames(CEU_filter))
match("rs7277536_T",colnames(CEU_filter))
match("rs7277536_C",colnames(CEU_filter))
match("rs7277536_G",colnames(CEU_filter))
match("rs7278456_A",colnames(CEU_filter))
match("rs7278456_T",colnames(CEU_filter))
match("rs7278456_C",colnames(CEU_filter))
match("rs7278456_G",colnames(CEU_filter))
match("rs2248610_A",colnames(CEU_filter))
match("rs2248610_T",colnames(CEU_filter))
match("rs2248610_C",colnames(CEU_filter))
match("rs2248610_G",colnames(CEU_filter))
match('rs965951_A',colnames(CEU_filter))
match('rs965951_T',colnames(CEU_filter))
match('rs965951_C',colnames(CEU_filter))
match('rs965951_G',colnames(CEU_filter))
match("rs3787662_A",colnames(CEU_filter))
match("rs3787662_T",colnames(CEU_filter))
match("rs3787662_C",colnames(CEU_filter))
match("rs3787662_G",colnames(CEU_filter))
match("rs2832253_A",colnames(CEU_filter))
match("rs2832253_T",colnames(CEU_filter))
match("rs2832253_C",colnames(CEU_filter))
match("rs2832253_G",colnames(CEU_filter))
match("rs2832332_A",colnames(CEU_filter))
match("rs2832332_T",colnames(CEU_filter))
match("rs2832332_C",colnames(CEU_filter))
match("rs2832332_G",colnames(CEU_filter))
match("rs13046799_A",colnames(CEU_filter))
match("rs13046799_T",colnames(CEU_filter))
match("rs13046799_C",colnames(CEU_filter))
match("rs13046799_G",colnames(CEU_filter))

true=c("rs2831459","rs7277536","rs13046799","rs2832332","rs2832253","rs3787662",
       'rs965951',"rs2248610","rs7278456")
count=0
for(l in 1:9){
  count=count+4-sum(is.na(c(match(paste0(true[l],"_A"),colnames(CEU_filter)[abs(act1)]),
                            match(paste0(true[l],"_T"),colnames(CEU_filter)[abs(act1)]),
                            match(paste0(true[l],"_C"),colnames(CEU_filter)[abs(act1)]),
                            match(paste0(true[l],"_G"),colnames(CEU_filter)[abs(act1)]))))
}
count


true_included_CEU=c("rs2831459_C","rs7278456_G",'rs965951_A',"rs2832253_T",
                    "rs2832332_G","rs13046799_A")


match("rs2831459_A",colnames(CEU_filter)[abs(act1)])
match("rs2831459_T",colnames(CEU_filter)[abs(act1)])
match("rs2831459_C",colnames(CEU_filter)[abs(act1)])
match("rs2831459_G",colnames(CEU_filter)[abs(act1)])
match("rs7277536_A",colnames(CEU_filter)[abs(act1)])
match("rs7277536_T",colnames(CEU_filter)[abs(act1)])
match("rs7277536_C",colnames(CEU_filter)[abs(act1)])
match("rs7277536_G",colnames(CEU_filter)[abs(act1)])
match("rs7278456_A",colnames(CEU_filter)[abs(act1)])
match("rs7278456_T",colnames(CEU_filter)[abs(act1)])
match("rs7278456_C",colnames(CEU_filter)[abs(act1)])
match("rs7278456_G",colnames(CEU_filter)[abs(act1)])
match("rs2248610_A",colnames(CEU_filter)[abs(act1)])
match("rs2248610_T",colnames(CEU_filter)[abs(act1)])
match("rs2248610_C",colnames(CEU_filter)[abs(act1)])
match("rs2248610_G",colnames(CEU_filter)[abs(act1)])
match('rs965951_A',colnames(CEU_filter)[abs(act1)])
match('rs965951_T',colnames(CEU_filter)[abs(act1)])
match('rs965951_C',colnames(CEU_filter)[abs(act1)])
match('rs965951_G',colnames(CEU_filter)[abs(act1)])
match("rs3787662_A",colnames(CEU_filter)[abs(act1)])
match("rs3787662_T",colnames(CEU_filter)[abs(act1)])
match("rs3787662_C",colnames(CEU_filter)[abs(act1)])
match("rs3787662_G",colnames(CEU_filter)[abs(act1)])
match("rs2832253_A",colnames(CEU_filter)[abs(act1)])
match("rs2832253_T",colnames(CEU_filter)[abs(act1)])
match("rs2832253_C",colnames(CEU_filter)[abs(act1)])
match("rs2832253_G",colnames(CEU_filter)[abs(act1)])
match("rs2832332_A",colnames(CEU_filter)[abs(act1)])
match("rs2832332_T",colnames(CEU_filter)[abs(act1)])
match("rs2832332_C",colnames(CEU_filter)[abs(act1)])
match("rs2832332_G",colnames(CEU_filter)[abs(act1)])
match("rs13046799_A",colnames(CEU_filter)[abs(act1)])
match("rs13046799_T",colnames(CEU_filter)[abs(act1)])
match("rs13046799_C",colnames(CEU_filter)[abs(act1)])
match("rs13046799_G",colnames(CEU_filter)[abs(act1)])
