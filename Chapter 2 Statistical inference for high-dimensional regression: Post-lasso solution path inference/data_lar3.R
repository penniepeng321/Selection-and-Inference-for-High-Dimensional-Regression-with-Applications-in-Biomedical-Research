#simulation 
rm(list=ls())

#setwd("/mnt/home/hpeng2")

library(covTest)
library(lars)
library(MASS)
library(hdi)


YRI_raw <- read.table("/home/huimin/Documents/YRI.raw", header = T, quote="\"")
YRI_name <- YRI_raw$IID
YRI_data <- YRI_raw[, -(1 : 6)]
YRI_filter <- YRI_data[, is.na(colSums(YRI_raw[, -(1 : 6)])) == F]

true_included_YRI=c("rs9982023_G","rs1236427_C","rs2831972_G","rs2091966_T",
                    "rs2832024_A" ,"rs2205413_T", "rs2832053_A" ,
                    "rs16983288_T", "rs16983303_A" ,"rs8134601_A" ,
                    "rs1006903_G" ,"rs7277685_C" ,"rs9982426_C", "rs2832115_G" ,
                    "rs2243503_G", "rs2243552_T" ,"rs2247809_A", "rs878797_A" ,
                    "rs8128844_T" ,"rs965951_A", "rs2832159_A", "rs2832178_T",
                    "rs2832186_G" ,"rs2832190_A" ,"rs7275293_T" ,"rs2251381_T",
                    "rs2251517_G", "rs2832225_A")

filter_index=match(true_included_YRI[1],colnames(YRI_filter))
for( hjk in 2:length(true_included_YRI)){
  filter_index=append(filter_index,match(true_included_YRI[hjk],colnames(YRI_filter)))
}
set.seed(555)
temp=sample(setdiff(colnames(YRI_filter),true_included_YRI),
                                   2000-length(true_included_YRI))
for(hjk in 1:length(temp)){
  filter_index=append(filter_index,match(temp[hjk],colnames(YRI_filter)))
}
length(filter_index)



YRI_filter <- YRI_filter[filter_index]


namelist <- c("NA18501","NA18502","NA18504","NA18505","NA18507","NA18508","NA18516",
              "NA18517","NA18522","NA18523","NA18852","NA18853","NA18855","NA18856",
              "NA18858","NA18859","NA18861","NA18862","NA18870","NA18871","NA18912",
              "NA18913","NA19092","NA19093","NA19098","NA19099","NA19101","NA19102",
              "NA19116","NA19119","NA19127","NA19130","NA19131","NA19137","NA19138",
              "NA19141","NA19143","NA19144","NA19152","NA19153","NA19159","NA19160",
              "NA19171","NA19172","NA19192","NA19193","NA19200","NA19201","NA19203",
              "NA19204","NA19206","NA19207","NA19209","NA19210","NA19222","NA19223",
              "NA19238","NA19239","NA19128","NA19140")

y <- as.numeric(c("13.2375466794262","13.8666316064246","13.4655762561699","13.4374899740909",
                  "13.1982215235080","13.4790073827739","13.6118354910471","13.6530547845626",
                  "13.3191710024573","13.4377736437223","13.2000824646128","13.2937598492530",
                  "13.1657753216251","13.9501070585572","13.6190447208569","13.0371567957961",
                  "13.2159988193338","13.8584159738624","13.2066214537857","13.7205135542957",
                  "13.6222374129441","13.3799665723633","13.2520444888282","13.1976767854297",
                  "13.2193245716207","13.1271384560751","13.6966409864918","13.4324185775067",
                  "13.8218242517983","13.3577530080554","13.4153545941554","13.7028607491570",
                  "13.6494878664414","13.0311549336935","13.9243452058029","13.5544190903444",
                  "13.6725585098348","13.6888593534800","13.3842738929415","13.4001842245995",
                  "13.4696023133097","13.3271346036410","13.5566339811415","13.7541549684152",
                  "13.7133025198463","13.6383335549412","13.6556216547411","13.4113206801535",
                  "13.5824214051641","13.4899802520369","13.8085968045668","13.1940579362263",
                  "13.0626421123032","13.1358052671756","13.3464927205494","13.3601911432172",
                  "13.2548005050623","13.4645532431137","13.3716651291550","13.2611088888499"))

y_order <- y[order(namelist)]




############################################################

n=length(y_order)
p=ncol(YRI_filter)
n;p

X=matrix(as.numeric(unlist(YRI_filter)),n,p,byrow=F)
dim(X)

y=as.numeric(unlist(y_order))
length(y)



#covariance test p values
pvalue3=function(X,y){
  n=nrow(X)
  p=ncol(X)
  X=X-matrix(colSums(X),n,p,byrow=T)
  Xcol=sqrt(colMeans(X^2))
  for(h in 1:p){
    X[,h]=X[,h]/Xcol[h]
  }
  y=y-rep(mean(y),n)
  for(h in 1:n){
    y[h]=y[h]/sqrt(mean(y^2))
  }
  l1=lars(X,y,type="lar",use.Gram=F)
  #l1=glmnet(X,y,family="gaussian")
  #sequential index of added varibles
  LARpath=as.numeric(unlist(l1$actions))
  #knots
  knot=l1$lambda
  
  temp=(LARpath>0)
  newlength=sum(temp)
  LARpath=LARpath[temp]
  #p+1 row p col entire rows are A_0,A_1,...,A_p
  #l1$beta displays a matrix of parameter estimates along LAR
  #the size of the display is (length(knot)+1)*p
  #each row is a vector of estimates at a certain knot
  #betahat=matrix(as.numeric(l1$beta),nrow=(length(knot)+1),ncol=p)
  #num of ts= m in the paper
  num_ts=length(knot)-1
  #store test statistic
  covTS=rep(0,num_ts)
  covtscomp=function(i){
    knot[i]*(knot[i]-knot[i+1]) 
  }
  covTS=sapply(1:num_ts,covtscomp)
  #at the last knot, the covariance ts is not well defined
  # and is thus not computed
  #q^* to ensure monotonicity
  qstar=rep(0,num_ts)
  qstarcomp=function(i){
    exp(-sum(covTS[i:num_ts]))
  }
  qstar=sapply(1:num_ts,qstarcomp)
  #apply two-cut method to qstar
  #because qstar are like ordered independent p values
  pvalue=list(model=l1,LAR=LARpath,q_value=qstar,covT=covTS,knot=knot)
  return(pvalue)
}

Lk=pvalue3(X,y)
qts1=Lk$q_value
act1=Lk$LAR

Lk$knot
Lk$covT

tempcovt=Lk$covT
covtspvalue=rep(NA,20)
for(plm in 1:20){
  covtspvalue[plm]=pexp(tempcovt[plm],rate=plm,lower.tail = F)
}
covtspvalue

qts1

lpm=length(qts1)

#bounding sequence
m=lpm
sim=1000
Vm=rep(0,sim)
set.seed(4567)
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


ppp=1-1/sqrt(log(lpm))
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=ceiling(max(((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1)[1:round(lpm/2)])*lpm)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]




(colnames(YRI_filter)[abs(act1)])[1:NPE]


ppp=1-0.1
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]


(colnames(YRI_filter)[abs(act1)])[1:NPE]


ppp=1-0.3
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]


(colnames(YRI_filter)[abs(act1)])[1:NPE]


ppp=1-0.5
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]


(colnames(YRI_filter)[abs(act1)])[1:NPE]


ppp=1-0.7
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]


(colnames(YRI_filter)[abs(act1)])[1:NPE]


ppp=1-0.9
Vm=as.numeric(unlist(Vm))
cbd=as.numeric(quantile(Vm,max(ppp,0)))
NPE=max(ceiling(max(((((1:lpm)/lpm-cbd*sqrt(qts1*(1-qts1))-qts1))/(1-qts1))[1:round(lpm/2)])*m)
        ,0)

ppp;cbd;NPE

act1

colnames(YRI_filter)[abs(act1)]


(colnames(YRI_filter)[abs(act1)])[1:NPE]


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


(colnames(YRI_filter)[abs(act1)])[1:NPE]
(colnames(YRI_filter)[abs(act1)])[1:LCV]
(colnames(YRI_filter)[abs(act1)])[1:BIC]
(colnames(YRI_filter)[abs(act1)])[1:BON]
(colnames(YRI_filter)[abs(act1)])[1:FDR]


true=c("rs9982023" ,"rs1236427" ,"rs2831972", "rs2091966" ,"rs2832010",
       "rs2832024" ,"rs2205413", "rs2832042", "rs2832053" ,"rs8130766",
       "rs16983288", "rs16983303" ,"rs8134601" ,"rs7276141" ,"rs7281691",
       "rs1006903" ,"rs7277685" ,"rs9982426", "rs2832115" ,"rs11910981",
       "rs2243503", "rs2243552" ,"rs2247809", "rs878797" ,"rs6516887",
       "rs8128844" ,"rs965951", "rs2070610" ,"rs2832159", "rs2832178",
       "rs2832186" ,"rs2832190" ,"rs7275293" ,"rs16983792", "rs2251381",
       "rs2251517", "rs2832225" ,"rs7283854")

true_included_YRI=c("rs9982023_G","rs1236427_C","rs2831972_G","rs2091966_T",
                    "rs2832024_A" ,"rs2205413_T", "rs2832053_A" ,
                    "rs16983288_T", "rs16983303_A" ,"rs8134601_A" ,
                    "rs1006903_G" ,"rs7277685_C" ,"rs9982426_C", "rs2832115_G" ,
                    "rs2243503_G", "rs2243552_T" ,"rs2247809_A", "rs878797_A" ,
                    "rs8128844_T" ,"rs965951_A", "rs2832159_A", "rs2832178_T",
                    "rs2832186_G" ,"rs2832190_A" ,"rs7275293_T" ,"rs2251381_T",
                    "rs2251517_G", "rs2832225_A")



wip=length(true)
for(l in 1:5){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
  match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
  match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
  match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 6:10){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 11:15){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 16:20){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 21:25){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 26:30){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 31:35){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

for(l in 36:38){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)])))
}

count=0
for(l in 1:38){
  count=count+4-sum(is.na(c(match(paste0(true[l],"_A"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_T"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_C"),colnames(YRI_filter)[abs(act1)]),
          match(paste0(true[l],"_G"),colnames(YRI_filter)[abs(act1)]))))
}
count


for(l in 1:5){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 6:10){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 11:15){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 16:20){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 21:25){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 26:30){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 31:35){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}

for(l in 36:38){
  print(c(match(paste0(true[l],"_A"),colnames(YRI_filter)),
          match(paste0(true[l],"_T"),colnames(YRI_filter)),
          match(paste0(true[l],"_C"),colnames(YRI_filter)),
          match(paste0(true[l],"_G"),colnames(YRI_filter))))
}



