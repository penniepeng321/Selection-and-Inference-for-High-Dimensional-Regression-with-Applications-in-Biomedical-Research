#par(mfrow=c(2,4), mar=c(4.5,4.5,2,0.5), oma=c(0,0,0,0))


par(mfrow=c(2,4), mar=c(4.5,4.5,2,0.5), oma=c(0,0,0,0))
#lars+QVS under high-dim setting
#plots 17
n=200;p=100;s=5
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)
l1=c(0.6070,0.9280,0.9962,1,1#,1,1
)
l2=c(0.2472,0.8468,0.9926,1,1#,1,1
)
l3=c(0.1788,0.6294,0.9340,1,1#,1,1
)
l4=c(0.2194,0.7370,0.9716,1,1#,1,1
)
l5=c(0.5328,0.9022,0.9950,1,1#,1,1
)
h1=c(0.0255,0.0224,0.0190,0.0177,0.0173#,0.0172,0.0171
)
h2=c(0.0094,0.0250,0.0297,0.0300,0.0298#,0.0302,0.0300
)
h3=c(0.0014,0.0034,0.0043,0.0045,0.0045#,0.0045,0.0045
)
h4=c(0.0025,0.0067,0.0084,0.0085,0.0085#,0.0085,0.0085
)
h5=c(0.0172,0.0209,0.0230,0.0234,0.0234#,0.0233,0.0227
)

gms1=sqrt(l1*h1)
gms2=sqrt(l2*h2)
gms3=sqrt(l3*h3)
gms4=sqrt(l4*h4)
gms5=sqrt(l5*h5)

plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
legend("bottomright", 
       c("BIC","LCV","BON","FDR","QVS"), 
       pch=c("I","C","","",""),
       lty=c(2,3,4,5,1),
       lwd=c(1,1,2,2,2),cex=1)



plot(bs,h1,pch="I",ylim=c(0,0.04),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)




# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))
# 
# plot(1:10,1:10)
# text(x=1,y=1,"T",cex=2)
# mtext("T",  WEST <0, at=9, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=130, cex=1.5, col="blue")
# 
# plot(1:20,1:20)
# text(x=10,y=10,"T",cex=2)
# mtext("F",  WEST <0, at=9, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=120, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=120, line=547, cex=1.5, col="blue")

n=200;p=100;s=5
m=min(n-1,p)
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)

B1=read.table('B1.csv',header=T)
B2=read.table('B2.csv',header=T)
B3=read.table('B3.csv',header=T)
B4=read.table('B4.csv',header=T)
B5=read.table('B5.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.2,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)


FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0,0.6),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)




#lars+QVS under high-dim setting
#plots 17
n=200;p=100;s=10
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)
l1=c(0.5699,0.8985,0.9871,1,1#,1,1
)
l2=c(0.4827,0.9453,0.9966,1,1#,1,1
)
l3=c(0.1838,0.5909,0.8892,0.9833,1#,1,1
)
l4=c(0.2663,0.7493,0.9613,1,1#,1,1
)
l5=c(0.5217,0.8690,0.9838,1,1#,1,1
)
h1=c(0.0326,0.0342,0.0274,0.0207,0.0176#,0.0159,0.0149
)
h2=c(0.0395,0.0812,0.0887,0.0893,0.0892#,0.0896,0.0898
)
h3=c(0.0037,0.0095,0.0137,0.0149,0.0151#,0.0152,0.0152
)
h4=c(0.0076,0.0204,0.0249,0.0258,0.0260#,0.0260,0.0260
)
h5=c(0.0293,0.0365,0.0391,0.0402,0.0403#,0.0403,0.0403
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.10),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# mtext("T",  WEST <0, at=230, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=230, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=350, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=350, line=547, cex=1.5, col="blue")
n=200;p=100;s=10
m=min(n-1,p)
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)
B1=read.table('B11.csv',header=T)
B2=read.table('B21.csv',header=T)
B3=read.table('B31.csv',header=T)
B4=read.table('B41.csv',header=T)
B5=read.table('B51.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.2,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
#dev.off()


FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0,0.6),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)



#lars+QVS under high-dim setting
#plots 17
n=200;p=100;s=20
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)
l1=c(0.5145,0.8448,0.9703,0.9969,1#,1,1
)
l2=c(0.7119,0.9659,0.9979,1,1#,1,1
)
l3=c(0.2285,0.5702,0.8225,0.9474,0.9905#,0.9983,1
)
l4=c(0.3678,0.7577,0.9397,0.9912,1#,1,1
)
l5=c(0.5287,0.8249,0.9588,0.9948,1#,1,1
)
h1=c(0.0584,0.0814,0.0840,0.0738,0.0630#,0.0556,0.0514
)
h2=c(0.1410,0.1979,0.2086,0.2086,0.2085#,0.2088,0.2091
)
h3=c(0.0179,0.0415,0.0570,0.0646,0.0679#,0.0691,0.0697
)
h4=c(0.0377,0.0715,0.0865,0.0907,0.0925#,0.0933,0.0940
)
h5=c(0.0688,0.0880,0.0957,0.0991,0.1013#,0.1022,0.1028
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.25),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# 
# mtext("T",  WEST <0, at=455, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=455, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=580, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=580, line=547, cex=1.5, col="blue")
n=200;p=100;s=20
m=min(n-1,p)
bs=c(0.2,0.3,0.4,0.5,0.6#,0.7,0.8
)

B1=read.table('B12.csv',header=T)
B2=read.table('B22.csv',header=T)
B3=read.table('B32.csv',header=T)
B4=read.table('B42.csv',header=T)
B5=read.table('B52.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.2,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0,0.6),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)


par(mfrow=c(2,4), mar=c(4.5,4.5,2,0.5), oma=c(0,0,0,0))

#lars+QVS under high-dim setting
#plots 17
n=200;p=400;s=10
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
l1=c(0.8069,0.9661,1,1,1#,1,1
)
l2=c(0.8516,0.9858,1,1,1#,1,1
)
l3=c(0.5332,0.8459,0.9699,1,1#,1,1
)
l4=c(0.6925,0.9384,0.9947,1,1#,1,1
)
l5=c(0.7943,0.9600,1,1,1#,1,1
)
h1=c(0.0205,0.0176,0.0121,0.0091,0.0079#,0.0072,0.0069
)
h2=c(0.0509,0.0654,0.0673,0.0673,0.0672#,0.0675,0.0672
)
h3=c(0.0053,0.0083,0.0094,0.0096,0.0097#,0.0097,0.0097
)
h4=c(0.0131,0.0181,0.0192,0.0193,0.0193#,0.0193,0.0193
)
h5=c(0.0247,0.0258,0.0256,0.0256,0.0256#,0.0256,0.0256
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
legend("bottomright", 
       c("BIC","LCV","BON","FDR","QVS"), 
       pch=c("I","C","","",""),
       lty=c(2,3,4,5,1),
       lwd=c(1,1,2,2,2),cex=1)

plot(bs,h1,pch="I",ylim=c(0,0.08),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# mtext("T",  WEST <0, at=9, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=9, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=120, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=120, line=547, cex=1.5, col="blue")

n=200;p=400;s=10
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
m=min(n-1,p)

B1=read.table('C1.csv',header=T)
B2=read.table('C2.csv',header=T)
B3=read.table('C3.csv',header=T)
B4=read.table('C4.csv',header=T)
B5=read.table('C5.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.6,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0,0.6),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)

#lars+QVS under high-dim setting
#plots 17
n=200;p=400;s=20
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
l1=c(0.6870,0.8885,0.9732,1,1#,1,1
)
l2=c(0.8614,0.9789,1,1,1#,1,1
)
l3=c(0.4858,0.7426,0.8934,0.9698,1#,1,1
)
l4=c(0.6676,0.8795,0.9730,1,1#,1,1
)
l5=c(0.7311,0.8965,0.9752,1,1#,1,1
)
h1=c(0.0414,0.0509,0.0506,0.0429,0.0359#,0.0309,0.0276
)
h2=c(0.1338,0.1712,0.1804,0.1819,0.1821#,0.1821,0.1819
)
h3=c(0.0194,0.0311,0.0390,0.0427,0.0441#,0.0447,0.0450
)
h4=c(0.0436,0.0599,0.0664,0.0681,0.0686#,0.0690,0.0693
)
h5=c(0.0583,0.0669,0.0689,0.0694,0.0696#,0.0699,0.0702
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.20),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# mtext("T",  WEST <0, at=230, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=230, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=350, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=350, line=547, cex=1.5, col="blue")

n=200;p=400;s=20
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
m=min(n-1,p)

B1=read.table('C11.csv',header=T)
B2=read.table('C21.csv',header=T)
B3=read.table('C31.csv',header=T)
B4=read.table('C41.csv',header=T)
B5=read.table('C51.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.6,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0,0.6),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)
#dev.off()

#lars+QVS under high-dim setting
#plots 17
n=200;p=400;s=40
bs=c(0.3,0.5,0.7,0.9,1.1#,1.3,1.5
)
l1=c(0.5520,0.8575,0.9764,1,1#,1,1
)
l2=c(0.8042,0.9813,1,1,1#,1,1
)
l3=c(0.4230,0.7181,0.8825,0.9598,0.9884#,1,1
)
l4=c(0.5980,0.8634,0.9678,1,1#,1,1
)
l5=c(0.6309,0.8575,0.9572,0.9901,1#,1,1
)
h1=c(0.0957,0.1746,0.2096,0.2061,0.1942#,0.1861,0.1804
)
h2=c(0.2682,0.3951,0.4248,0.4312,0.4337#,0.4349,0.4344
)
h3=c(0.0626,0.1266,0.1667,0.1880,0.1982#,0.2034,0.2062
)
h4=c(0.1190,0.1931,0.2241,0.2351,0.2401#,0.2437,0.2461
)
h5=c(0.1336,0.1901,0.2136,0.2226,0.2270#,0.2304,0.2326
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.45),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# mtext("T",  WEST <0, at=455, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=455, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=580, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=580, line=547, cex=1.5, col="blue")
n=200;p=400;s=40
bs=c(0.3,0.5,0.7,0.9,1.1#,1.3,1.5
)
m=min(n-1,p)

B1=read.table('C12.csv',header=T)
B2=read.table('C22.csv',header=T)
B3=read.table('C32.csv',header=T)
B4=read.table('C42.csv',header=T)
B5=read.table('C52.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.6,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)

FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.2,0.8),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)


par(mfrow=c(2,4), mar=c(4.5,4.5,2,0.5), oma=c(0,0,0,0))

#lars+QVS under high-dim setting
#plots 17
n=200;p=2000;s=10
bs=c(0.3,0.4,0.5,0.6#,0.7#,0.8,0.9
)
l1=c(0.6438,0.8795,0.9764,0.9976#,0.9998#,1,1
)
l2=c(0.5809,0.9220,0.9914,0.9997#,1.0000#,1,1
)
l3=c(0.5333,0.8289,0.9665,0.9961#,0.9998#,1,1
)
l4=c(0.7312,0.9412,0.9935,0.9997#,1.0000#,1,1
)
l5=c(0.7597,0.9373,0.9920,0.9994#,1.0000#,1,1
)
h1=c(0.0254,0.0243,0.0183,0.0124#,0.0092#,0.0073,0.0065
)
h2=c(0.0493,0.0902,0.1052,0.1078#,0.1081#,0.1080,0.1083
)
h3=c(0.0156,0.0235,0.0269,0.0276#,0.0277#,0.0278,0.0278
)
h4=c(0.0529,0.0679,0.0711,0.0714#,0.0715#,0.0715,0.0715
)
h5=c(0.0642,0.0679,0.0667,0.0666#,0.0666#,0.0666,0.0666
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
legend("bottomright", 
       c("BIC","LCV","BON","FDR","QVS"), 
       pch=c("I","C","","",""),
       lty=c(2,3,4,5,1),
       lwd=c(1,1,2,2,2),cex=1)

plot(bs,h1,pch="I",ylim=c(0,0.15),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))
n=200;p=2000;s=10
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
m=min(n-1,p)

B1=read.table('A38.csv',header=T)
B2=read.table('A39.csv',header=T)
B3=read.table('A40.csv',header=T)
B4=read.table('A41.csv',header=T)
B5=read.table('A42.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.3,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.1,0.7),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)

#lars+QVS under high-dim setting
#plots 17
n=200;p=2000;s=20
bs=c(0.3,0.4,0.5,0.6,0.7,0.8#,0.9,1.1
)
l1=c(0.5025,0.6791,0.8251,0.9164,0.9677,0.9887#,0.9965,1
)
l2=c(0.4984,0.8078,0.9458,0.9864,0.9966,0.9992#,0.9999,1
)
l3=c(0.5023,0.7105,0.8551,0.9403,0.9788,0.9936#,0.9981,1
)
l4=c(0.7080,0.8670,0.9506,0.9856,0.9965,0.9993#,0.9999,1
)
l5=c(0.7063,0.8467,0.9314,0.9774,0.9933,0.9980#,0.9996,1
)
h1=c(0.0492,0.0658,0.0815,0.0887,0.0906,0.0854#,0.0790,0.0666
)
h2=c(0.0901,0.1896,0.2633,0.3011,0.3168,0.3238#,0.3262,0.3282
)
h3=c(0.0535,0.0881,0.1131,0.1275,0.1344,0.1373#,0.1386,0.1401
)
h4=c(0.1464,0.1941,0.2183,0.2280,0.2315,0.2327#,0.2333,0.2344
)
h5=c(0.1442,0.1738,0.1876,0.1922,0.1941,0.1947#,0.1952,0.1969
)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.35),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))
n=200;p=2000;s=20
bs=c(0.3,0.4,0.5,0.6,0.7,0.8#,0.9,1.1
)
m=min(n-1,p)

B1=read.table('A45.csv',header=T)
B2=read.table('A46.csv',header=T)
B3=read.table('A47.csv',header=T)
B4=read.table('A48.csv',header=T)
B5=read.table('A49.csv',header=T)
B6=read.table('A50.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
B6=as.matrix(B6)
Bmat=cbind(B1,B2,B3,B4,B5,B6)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.3,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.2,0.8),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)

#lars+QVS under high-dim setting
#plots 17
n=200;p=2000;s=40
bs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9)
l1=c(0.4390,0.5489,0.6301,0.6947,0.7441,0.7833,0.8118,0.8600,0.8929,0.9148,0.9327,0.9427)
l2=c(0.3014,0.4470,0.5388,0.5982,0.6350,0.6631,0.6811,0.7030,0.7182,0.7250,0.7334,0.7372)
l3=c(0.5152,0.6443,0.7286,0.7852,0.8254,0.8541,0.8773,0.9109,0.9328,0.9469,0.9574,0.9642)
l4=c(0.7022,0.7987,0.8511,0.8886,0.9115,0.9289,0.9428,0.9594,0.9700,0.9768,0.9829,0.9847)
l5=c(0.6729,0.7572,0.7954,0.8156,0.8232,0.8250,0.8262,0.8268,0.8272,0.8273,0.8282,0.8273)
h1=c(0.1069,0.1605,0.2129,0.2666,0.3135,0.3581,0.3950,0.4655,0.5214,0.5665,0.6050,0.6366)
h2=c(0.0935,0.1576,0.2091,0.2506,0.2783,0.3008,0.3166,0.3375,0.3503,0.3591,0.3654,0.3692)
h3=c(0.1512,0.2356,0.3115,0.3747,0.4283,0.4735,0.5143,0.5785,0.6278,0.6681,0.6998,0.7258)
h4=c(0.3085,0.4112,0.4901,0.5524,0.6014,0.6417,0.6756,0.7277,0.7646,0.7940,0.8164,0.8348)
h5=c(0.2797,0.3542,0.3993,0.4194,0.4274,0.4291,0.4283,0.4266,0.4256,0.4249,0.4243,0.4242)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.85),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))


# mtext("T",  WEST <0, at=9, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=9, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=120, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=120, line=547, cex=1.5, col="blue")
# 
# mtext("T",  WEST <0, at=230, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=230, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=340, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=340, line=547, cex=1.5, col="blue")
# 
# 
# mtext("T",  WEST <0, at=455, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=455, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=580, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=580, line=547, cex=1.5, col="blue")
n=200;p=2000;s=40
bs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9)
m=min(n-1,p)

B1=read.table('A52.csv',header=T)
B2=read.table('A53.csv',header=T)
B3=read.table('A54.csv',header=T)
B4=read.table('A55.csv',header=T)
B5=read.table('A56.csv',header=T)
B6=read.table('A57.csv',header=T)
B7=read.table('A58.csv',header=T)
B8=read.table('D1.csv',header=T)
B9=read.table('D2.csv',header=T)
B10=read.table('D3.csv',header=T)
B11=read.table('D4.csv',header=T)
B12=read.table('D5.csv',header=T)


B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
B6=as.matrix(B6)
B7=as.matrix(B7)
B8=as.matrix(B8)
B9=as.matrix(B9)
B10=as.matrix(B10)
B11=as.matrix(B11)
B12=as.matrix(B12)
Bmat=cbind(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.3,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.4,1),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)


par(mfrow=c(2,4), mar=c(4.5,4.5,2,0.5), oma=c(0,0,0,0))

#lars+QVS under high-dim setting
#plots 17
n=200;p=10000;s=10
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
l1=c(0.4987,0.7258,0.8924,0.9700,0.9945#,0.9992,1
)
l2=c(0.2828,0.7043,0.9278,0.9881,0.9993#,0.9999,1
)
l3=c(0.5671,0.8188,0.9503,0.9904,0.9990#,0.9998,1
)
l4=c(0.7918,0.9336,0.9869,0.9984,0.9999#,1.0000,1
)
l5=c(0.7086,0.8882,0.9727,0.9954,0.9997#,0.9999,1
)
h1=c(0.0288,0.0307,0.0286,0.0228,0.0166#,0.0114,0.0083
)
h2=c(0.0267,0.0739,0.1187,0.1393,0.1455#,0.1470,0.1466
)
h3=c(0.0488,0.0730,0.0841,0.0875,0.0881#,0.0882,0.0882
)
h4=c(0.1780,0.2133,0.2238,0.2259,0.2261#,0.2262,0.2262
)
h5=c(0.1143,0.1342,0.1391,0.1400,0.1400#,0.1400,0.1400
)

plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
legend("bottomright", 
       c("BIC","LCV","BON","FDR","QVS"), 
       pch=c("I","C","","",""),
       lty=c(2,3,4,5,1),
       lwd=c(1,1,2,2,2),cex=1)

plot(bs,h1,pch="I",ylim=c(0,0.25),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))
n=200;p=10000;s=10
bs=c(0.3,0.4,0.5,0.6,0.7#,0.8,0.9
)
m=min(n-1,p)

B1=read.table('A80.csv',header=T)
B2=read.table('A81.csv',header=T)
B3=read.table('A82.csv',header=T)
B4=read.table('A83.csv',header=T)
B5=read.table('A84.csv',header=T)

B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
Bmat=cbind(B1,B2,B3,B4,B5)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.1,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.2,1),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)
#lars+QVS under high-dim setting
#plots 17
n=200;p=10000;s=20
bs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3)
l1=c(0.4136,0.5225,0.6131,0.6917,0.7516,0.8108,0.8587,0.9261,0.9601)
l2=c(0.1747,0.3617,0.5467,0.6677,0.7540,0.8180,0.8607,0.9102,0.9341)
l3=c(0.5739,0.7248,0.8179,0.8748,0.9173,0.9471,0.9663,0.9856,0.9935)
l4=c(0.7730,0.8642,0.9173,0.9472,0.9690,0.9807,0.9885,0.9957,0.9978)
l5=c(0.6880,0.7962,0.8616,0.9020,0.9301,0.9513,0.9631,0.9773,0.9828)
h1=c(0.0556,0.0756,0.0948,0.1178,0.1382,0.1586,0.1777,0.2030,0.2160)
h2=c(0.0287,0.0703,0.1256,0.1807,0.2313,0.2744,0.3072,0.3558,0.3821)
h3=c(0.1334,0.2093,0.2723,0.3209,0.3567,0.3812,0.3988,0.4186,0.4292)
h4=c(0.3264,0.4132,0.4751,0.5161,0.5424,0.5594,0.5708,0.5830,0.5891)
h5=c(0.2258,0.2939,0.3437,0.3745,0.3925,0.4017,0.4061,0.4087,0.4093)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,0.65),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))
n=200;p=10000;s=20
bs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3)
m=min(n-1,p)

B1=read.table('A87.csv',header=T)
B2=read.table('A88.csv',header=T)
B3=read.table('A89.csv',header=T)
B4=read.table('A90.csv',header=T)
B5=read.table('A91.csv',header=T)
B6=read.table('A92.csv',header=T)
B7=read.table('A93.csv',header=T)
B8=read.table('D13.csv',header=T)
B9=read.table('D14.csv',header=T)


B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
B6=as.matrix(B6)
B7=as.matrix(B7)
B8=as.matrix(B8)
B9=as.matrix(B9)
Bmat=cbind(B1,B2,B3,B4,B5,B6,B7,B8,B9)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.1,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)

FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.4,1),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)

#lars+QVS under high-dim setting
#plots 17
n=200;p=10000;s=40
bs=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9)
l1=c(0.3952,0.4884,0.5691,0.6333,0.6814,0.7196,0.7532,0.7995,0.8321,0.8577,0.8787,0.8917)
l2=c(0.0708,0.1049,0.1329,0.1575,0.1708,0.1783,0.1862,0.1957,0.2001,0.2059,0.2092,0.2103)
l3=c(0.6318,0.7409,0.8038,0.8484,0.8764,0.8997,0.9155,0.9336,0.9470,0.9598,0.9675,0.9690)
l4=c(0.7962,0.8576,0.8975,0.9239,0.9386,0.9494,0.9586,0.9675,0.9750,0.9814,0.9853,0.9854)
l5=c(0.7015,0.7640,0.7850,0.7921,0.7941,0.7945,0.7963,0.7967,0.7939,0.7971,0.7997,0.7978)
h1=c(0.1069,0.1560,0.2083,0.2618,0.3097,0.3553,0.3982,0.4711,0.5326,0.5804,0.6229,0.6561)
h2=c(0.0196,0.0280,0.0359,0.0438,0.0490,0.0514,0.0541,0.0572,0.0599,0.0612,0.0626,0.0635)
h3=c(0.2863,0.4046,0.4975,0.5703,0.6288,0.6750,0.7133,0.7706,0.8125,0.8436,0.8668,0.8851)
h4=c(0.5033,0.6115,0.6873,0.7426,0.7860,0.8175,0.8431,0.8791,0.9031,0.9211,0.9350,0.9438)
h5=c(0.3694,0.4391,0.4637,0.4684,0.4676,0.4670,0.4672,0.4658,0.4656,0.4652,0.4648,0.4648)
plot(bs,l1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="TPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,l1,lty=2,lwd=1)
points(bs,l2,pch="C")
lines(bs,l2,lty=3,lwd=1)
points(bs,l3,pch="")
lines(bs,l3,lty=4,lwd=2)
points(bs,l4,pch="")
lines(bs,l4,lty=5,lwd=2)
lines(bs,l5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

plot(bs,h1,pch="I",ylim=c(0,1),xlab="
     intensity",ylab="FPP",main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,h1,lty=2,lwd=1)
points(bs,h2,pch="C")
lines(bs,h2,lty=3,lwd=1)
points(bs,h3,pch="")
lines(bs,h3,lty=4,lwd=2)
points(bs,h4,pch="")
lines(bs,h4,lty=5,lwd=2)
lines(bs,h5,lty=1,lwd=2)
# legend("topright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,3))

# mtext("T",  WEST <0, at=9, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=9, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=9, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=120, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=120, line=547, cex=1.5, col="blue")
# 
# mtext("T",  WEST <0, at=230, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=230, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=230, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=340, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=340, line=547, cex=1.5, col="blue")
# 
# 
# mtext("T",  WEST <0, at=455, line=100, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=115, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=130, cex=1.5, col="blue")
# 
# 
# mtext("F",  WEST <0, at=455, line=380, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=395, cex=1.5, col="blue")
# mtext("P",  WEST <0, at=455, line=410, cex=1.5, col="blue")
# 
# 
# mtext("intensity",  SOUTH <0, at=580, line=275, cex=1.5, col="blue")
# mtext("intensity",  SOUTH <0, at=580, line=547, cex=1.5, col="blue")

n=200;p=10000;s=40
bs=c(0.3,0.4,0.5,0.6,0.7,0.8, 1.1,1.3,1.5,1.7,1.9
)
m=min(n-1,p)

B1=read.table('A94.csv',header=T)
B2=read.table('A95.csv',header=T)
B3=read.table('A96.csv',header=T)
B4=read.table('A97.csv',header=T)
B5=read.table('A98.csv',header=T)
B6=read.table('A99.csv',header=T)

B8=read.table('D6.csv',header=T)
B9=read.table('D7.csv',header=T)
B10=read.table('D8.csv',header=T)
B11=read.table('D9.csv',header=T)
B12=read.table('D10.csv',header=T)


B1=as.matrix(B1)
B2=as.matrix(B2)
B3=as.matrix(B3)
B4=as.matrix(B4)
B5=as.matrix(B5)
B6=as.matrix(B6)

B8=as.matrix(B8)
B9=as.matrix(B9)
B10=as.matrix(B10)
B11=as.matrix(B11)
B12=as.matrix(B12)
Bmat=cbind(B1,B2,B3,B4,B5,B6  ,B8,B9,B10,B11,B12
)

gmesmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    t1=B1[(B1[,16]!=0),2*k]/B1[(B1[,16]!=0),16]
    t2=1-(B1[(B1[,16]!=0),(2*k-1)]-B1[(B1[,16]!=0),(2*k)])/(m-B1[(B1[,16]!=0),16])
    gmesmat[k,j]=mean(sqrt(t1*t2))
  }
}

gmes1=gmesmat[1,]
gmes2=gmesmat[2,]
gmes3=gmesmat[3,]
gmes4=gmesmat[4,]
gmes5=gmesmat[5,]

plot(bs,gmes1,pch="I",ylim=c(0.1,1),xlab="
     intensity",ylab='g-measure',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.5, cex.axis = 1.3)
lines(bs,gmes1,lty=2,lwd=1)
points(bs,gmes2,pch="C")
lines(bs,gmes2,lty=3,lwd=1)
points(bs,gmes3,pch="")
lines(bs,gmes3,lty=4,lwd=2)
points(bs,gmes4,pch="")
lines(bs,gmes4,lty=5,lwd=2)
lines(bs,gmes5,lty=1,lwd=2)
# legend("bottomright", 
#        c("BIC","LCV","BON","FDR","QVS"), 
#        pch=c("I","C","","",""),
#        lty=c(2,3,4,5,1),
#        lwd=c(1,1,2,2,2),cex=1.5)
FDPmat=matrix(0,5,length(bs))

for(j in 1:length(bs)){
  B1=Bmat[,(16*(j-1)+1):(16*j)]
  for(k in 1:5){
    FDPmat[k,j]=mean(1-B1[(B1[,(2*k-1)]!=0),(2*k)]/B1[(B1[,(2*k-1)]!=0),(2*k-1)])
  }
}

FDP1=FDPmat[1,]
FDP2=FDPmat[2,]
FDP3=FDPmat[3,]
FDP4=FDPmat[4,]
FDP5=FDPmat[5,]

plot(bs,FDP1,pch="I",ylim=c(0.4,1),xlab="
     intensity",ylab='FDP',main=
       paste("n=",n,"/p=",p,"/s=",s,"/AR"),cex.lab = 1.8, cex.axis = 1.5)
lines(bs,FDP1,lty=2,lwd=1)
points(bs,FDP2,pch="C")
lines(bs,FDP2,lty=3,lwd=1)
points(bs,FDP3,pch="")
lines(bs,FDP3,lty=4,lwd=2)
points(bs,FDP4,pch="")
lines(bs,FDP4,lty=5,lwd=2)
lines(bs,FDP5,lty=1,lwd=2)
#dev.off()
# 
dev.off()