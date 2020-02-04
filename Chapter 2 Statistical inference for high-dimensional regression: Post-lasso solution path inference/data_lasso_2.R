#simulation 
rm(list=ls())

#setwd("/mnt/home/hpeng2")

library(covTest)
library(lars)
library(MASS)
library(hdi)


JPT_CHB_raw <- read.table("/home/huimin/Documents/JPT_CHB.raw", header = T, quote="\"")
JPT_CHB_name <- JPT_CHB_raw$IID
JPT_CHB_data <- JPT_CHB_raw[, -(1 : 6)]
JPT_CHB_filter <- JPT_CHB_data[, is.na(colSums(JPT_CHB_raw[, -(1 : 6)])) == F]

JPT_CHB_filter <- JPT_CHB_filter[, colSums(JPT_CHB_filter) >= 
                                   0.05*2  * nrow(JPT_CHB_filter)]


true_included_JPT_CHB=c("rs2832159_A","rs2245431_T","rs7282280_T")
filter_index=match(true_included_JPT_CHB[1],colnames(JPT_CHB_filter))
for( hjk in 2:length(true_included_JPT_CHB)){
  filter_index=append(filter_index,match(true_included_JPT_CHB[hjk],colnames(JPT_CHB_filter)))
}
temp=sample(setdiff(colnames(JPT_CHB_filter),true_included_JPT_CHB),
            5000-length(true_included_JPT_CHB))
for(hjk in 1:length(temp)){
  filter_index=append(filter_index,match(temp[hjk],colnames(JPT_CHB_filter)))
}
length(filter_index)

JPT_CHB_filter <- JPT_CHB_filter[filter_index]



namelist <- c("NA18940","NA18942","NA18943","NA18944","NA18945","NA18947","NA18948","NA18949",
              "NA18951","NA18952","NA18953","NA18956","NA18959","NA18960","NA18961","NA18964",
              "NA18965","NA18966","NA18967","NA18968","NA18969","NA18970","NA18971","NA18972",
              "NA18973","NA18975","NA18976","NA18978","NA18980","NA18981","NA18987","NA18990",
              "NA18991","NA18992","NA18994","NA18995","NA18997","NA18998","NA18999","NA19000",
              "NA19003","NA19005","NA19012","NA19007","NA18974","NA18524","NA18526","NA18529",
              "NA18537","NA18540","NA18542","NA18545","NA18547","NA18550","NA18552","NA18555",
              "NA18558","NA18561","NA18562","NA18563","NA18564","NA18566","NA18570","NA18571",
              "NA18572","NA18573","NA18576","NA18577","NA18579","NA18582","NA18592","NA18593",
              "NA18594","NA18603","NA18605","NA18608","NA18609","NA18611","NA18612","NA18620",
              "NA18621","NA18622","NA18623","NA18624","NA18632","NA18633","NA18635","NA18636",
              "NA18637","NA18532")

y <- as.numeric(c("13.3138935876198","13.4765768436313","13.4744110218750","13.3414682897444",
                  "13.6498970554683","13.2358992446565","13.7677520557194","13.2568844047249",
                  "13.4657487093833","13.2722417936417","13.4158245078951","13.3001165272561",
                  "12.8285402357208","13.2217094405126","13.2309554600056","13.2866547466829",
                  "13.4110391139191","13.3739491269624","13.0959716623119","13.5547789464538",
                  "13.39289466652","13.5483557111725","13.2494389397586","13.3760414161565",
                  "13.1095933999577","13.3604075485222","13.1418812789614","13.4671730617591",
                  "13.151081102449","13.5447545885485","13.7327201859884","13.4519370617366",
                  "13.3972972472058","13.4445921476462","13.3471690423879","13.4418473600311",
                  "13.4663437437659","13.5129819051744","13.6723795824062","13.2140909018436",
                  "13.4419548589422","13.1344690599592","13.3334431273999","13.9235074065227",
                  "13.3948121395240", "13.3251682878942","13.2795866114057","13.2251679202102",
                  "13.3528183560039","13.4786497515687","13.4780351365868","13.5203929348286",
                  "13.5246294236922","13.3357624938171","13.3959068747341","13.4630555610632",
                  "13.4428880707915","13.4014388374885","13.3334273633403","13.2502918739290",
                  "13.535508756536","13.2659841559861","13.3821410402617","13.4079710701753",
                  "13.4939433030891","13.5520663322600","13.769645123503","13.2822270944276",
                  "13.3568533211745","13.2201876945441","13.0983387091758","13.5459730562942",
                  "13.5940777576623","13.3668571880270","13.6217630710201","13.5292323189050",
                  "13.7327643612330","13.3151492461755","13.6426912055455","13.6584679151410",
                  "13.1607350904759","13.8394877953868","13.7448664206536","14.0739814291627",
                  "13.3539474256382","13.7118890058913","13.6875624510660","13.7448386231938",
                  "13.7918800278628","13.7738122104445"))

y_order <- y[order(namelist)]




############################################################

n=length(y_order)
p=ncol(JPT_CHB_filter)
n;p

X=matrix(as.numeric(unlist(JPT_CHB_filter)),n,p,byrow=F)
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

colnames(JPT_CHB_filter)[abs(act1)]




(colnames(JPT_CHB_filter)[abs(act1)])[1:NPE]



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


(colnames(JPT_CHB_filter)[abs(act1)])[1:NPE]
(colnames(JPT_CHB_filter)[abs(act1)])[1:LCV]
(colnames(JPT_CHB_filter)[abs(act1)])[1:BIC]
(colnames(JPT_CHB_filter)[abs(act1)])[1:BON]
(colnames(JPT_CHB_filter)[abs(act1)])[1:FDR]

match("rs16981663_A",colnames(JPT_CHB_filter))
match("rs16981663_T",colnames(JPT_CHB_filter))
match("rs16981663_C",colnames(JPT_CHB_filter))
match("rs16981663_G",colnames(JPT_CHB_filter))
match("rs9981984_A",colnames(JPT_CHB_filter))
match("rs9981984_T",colnames(JPT_CHB_filter))
match("rs9981984_C",colnames(JPT_CHB_filter))
match("rs9981984_G",colnames(JPT_CHB_filter))
match("rs7282280_A",colnames(JPT_CHB_filter))
match("rs7282280_T",colnames(JPT_CHB_filter))
match("rs7282280_C",colnames(JPT_CHB_filter))
match("rs7282280_G",colnames(JPT_CHB_filter))
match("rs2245431_A",colnames(JPT_CHB_filter))
match("rs2245431_T",colnames(JPT_CHB_filter))
match("rs2245431_C",colnames(JPT_CHB_filter))
match("rs2245431_G",colnames(JPT_CHB_filter))
match("rs2832159_A",colnames(JPT_CHB_filter))
match("rs2832159_T",colnames(JPT_CHB_filter))
match("rs2832159_C",colnames(JPT_CHB_filter))
match("rs2832159_G",colnames(JPT_CHB_filter))
match("rs1999321_A",colnames(JPT_CHB_filter))
match("rs1999321_T",colnames(JPT_CHB_filter))
match("rs1999321_C",colnames(JPT_CHB_filter))
match("rs1999321_G",colnames(JPT_CHB_filter))
match("rs2832224_A",colnames(JPT_CHB_filter))
match("rs2832224_T",colnames(JPT_CHB_filter))
match("rs2832224_C",colnames(JPT_CHB_filter))
match("rs2832224_G",colnames(JPT_CHB_filter))

true=c("rs1999321","rs2832224","rs16981663","rs9981984","rs7282280",
       "rs2245431","rs2832159")
count=0
for(l in 1:7){
  count=count+4-sum(is.na(c(match(paste0(true[l],"_A"),colnames(JPT_CHB_filter)[abs(act1)]),
                            match(paste0(true[l],"_T"),colnames(JPT_CHB_filter)[abs(act1)]),
                            match(paste0(true[l],"_C"),colnames(JPT_CHB_filter)[abs(act1)]),
                            match(paste0(true[l],"_G"),colnames(JPT_CHB_filter)[abs(act1)]))))
}
count

true_included_JPT_CHB=c("rs2832159_A","rs2245431_T","rs7282280_T")


match("rs16981663_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs16981663_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs16981663_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs16981663_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs9981984_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs9981984_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs9981984_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs9981984_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs7282280_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs7282280_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs7282280_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs7282280_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2245431_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2245431_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2245431_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2245431_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832159_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832159_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832159_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832159_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs1999321_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs1999321_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs1999321_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs1999321_G",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832224_A",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832224_T",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832224_C",colnames(JPT_CHB_filter)[abs(act1)])
match("rs2832224_G",colnames(JPT_CHB_filter)[abs(act1)])
