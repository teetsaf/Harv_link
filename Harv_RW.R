##########Climate

library(treeclim)
library(dplR)
library(dplyr)
library(ggplot2)
library(tidyr)

harv_cored2=read.csv('F:/NAU/R/Harv_cored3.csv')
head(harv_cored2)
colnames(harv_cored2)

###############################DATA

rw=read.delim('F:/NAU/Harv_cores/LyfordAllCombinedOrder.txt',header=FALSE)
head(rw)
#rw1=gather(rw,coreID,RW,)

acru=read.rwl('F:/NAU/Harv_cores/LF_ACRU.rwl')
head(acru)
colnames(acru)
acru$YR=rownames(acru)
acru1=gather(acru,coreID,RW,LF103s:LF349n,na.rm=TRUE)
head(acru1)
acru1$TreeID=substr(acru1$coreID,1,5)
acru1$corloc=substr(acru1$coreID,6,6)
head(acru1)


beal=read.rwl('F:/NAU/Harv_cores/LF_BEAL.rwl')
head(beal)
beal$YR=rownames(beal)
beal1=gather(beal,coreID,RW,LF107n:LF320s,na.rm=TRUE)
head(beal1)
beal1$TreeID=substr(beal1$coreID,1,5)
beal1$corloc=substr(beal1$coreID,6,6)
head(beal1)


fagr=read.rwl('F:/NAU/Harv_cores/LF_FAGR.rwl')
head(fagr)
fagr$YR=rownames(fagr)
fagr1=gather(fagr,coreID,RW,LF109s:LF354s,na.rm=TRUE)
head(fagr1)
fagr1$TreeID=substr(fagr1$coreID,1,5)
fagr1$corloc=substr(fagr1$coreID,6,6)
head(fagr1)


quru=read.rwl('F:/NAU/Harv_cores/LF_QURU.rwl')
head(quru)
quru$YR=rownames(quru)
quru1=gather(quru,coreID,RW,LF102n:LF372w,na.rm=TRUE)
head(quru1)
quru1$TreeID=substr(quru1$coreID,1,5)
quru1$corloc=substr(quru1$coreID,6,6)
head(quru1)


tsca=read.rwl('F:/NAU/Harv_cores/LF_TSCA.rwl')
head(tsca)
tsca$YR=rownames(tsca)
tsca1=gather(tsca,coreID,RW,LF101a:LF128w,na.rm=TRUE)
head(tsca1)
tsca1$TreeID=substr(tsca1$coreID,1,5)
tsca1$corloc=substr(tsca1$coreID,6,6)
head(tsca1)

com1=rbind(acru1,beal1)
head(com1)
com2=rbind(com1,fagr1)
com3=rbind(com2,quru1)
head(com3)
comb=rbind(com3,tsca1)


#com1=combine.rwl(ACRU,BEAL)
#com2=combine.rwl(com1,FAGR)
#com3=combine.rwl(com2,QURU)
#com=combine.rwl(com3,TSCA)


rwdb=merge(comb,harv_cored2,by="TreeID",all=TRUE)
head(rwdb)
rwdb2=rwdb[complete.cases(rwdb[,c(2,4)]),]


setwd=('F:/NAU/R')
#write.csv(rwdb2,'Harv_RW_DB.csv')


#reconstruct diameters

RW=rwdb2
str(RW)
RW$YR=as.numeric(RW$YR)
RW=subset(RW,YR>1979)

RW <- RW[order(RW$coreID,-RW$YR),] 
head(RW)

RW$SPP=RW$Species
RW$bfactor=1
a=1


###BARK factor function (Dixon and Keyser 2011) 
bfactor=function(SPP){
  if(SPP=='ABBA'){a=0.9349}
  else if(SPP=='PIRU'){a=0.9324}
  else if(SPP=='ACPE'){a=0.90}
  else if(SPP=='ACRU'){a=0.95}
  else if(SPP=='ACSA'){a=0.92}
  else if(SPP=='BEAL'){a=0.948}
  else if(SPP=='BEPA'){a=0.948}
  else if(SPP=='BEPO'){a=0.948}
  else if(SPP=='FAGR'){a=0.95}
  else if(SPP=='FRAM'){a=0.9 }
  else if(SPP=='FRNI'){a=0.9 }
  else if(SPP=='OSVI'){a=0.9 }
  else if(SPP=='PIGL'){a=0.956}
  else if(SPP=='PIRE'){a=0.92}
  else if(SPP=='PIST'){a=0.92}
  else if(SPP=='POGR'){a=0.9 }
  else if(SPP=='POTR'){a=0.9 }
  else if(SPP=='PRSE'){a=0.94}
  else if(SPP=='QURU'){a=0.9 }
  else if(SPP=='THOC'){a=0.95}
  else if(SPP=='TIAM'){a=0.9 }
  else if(SPP=='TSCA'){a=0.934}
  else if(SPP=='ULAM'){a=0.9}
  bfactor=(1-a)
  return(bfactor=bfactor)}

RW$bfactor=mapply(bfactor,RW$SPP)

head(RW)

#bark width at year cored
RW$bark.w15=(RW$DBH*RW$bfactor)

#diameter inside bark at year cored
RW$DIB15<-(RW$DBH-RW$bark.w15)

#convert ring widths from mm to cm
RW$RWcm=RW$RW/10
RW$TAG=RW$coreID
TAG<-as.factor(RW$TAG)  


#####  Yearly diameter inside bark; DIB=DIB15-cumulative sum of rings from 2015 to year 'length'
head(RW)
RW$DIB=0
for(tag in unique(TAG)){
  RW[RW$TAG==tag,]$DIB=RW[RW$TAG==tag,"DIB15"]-c(0,cumsum(2*RW[RW$TAG==tag,"RWcm"])[-length(RW[RW$TAG==tag,"RWcm"])])}
RW$DOB=RW$DIB/(1-RW$bfactor)
RW$bark.w=RW$DOB-RW$DIB


#Basal area inside bark (m^2)
RW$BA=((RW$DIB/2)^2*pi)/10000


#Basal area outside bark (m^2)
RW$BA.ob=((RW$DOB/2)^2*pi)/10000

head(RW)

#Basal area increment outside bark 
RW$BAI.ob=c(abs(diff(RW$BA.ob)),0)
for(tag in unique(TAG)){
  RW[RW$TAG==tag,"BAI.ob"][length(RW[RW$TAG==tag,"BA.ob"])]=0}


RW$cBAI=0
for(tag in unique(TAG)){
  RW[RW$TAG==tag,]$cBAI=sum(RW[RW$TAG==tag,]$BA.ob)[-length(RW[RW$TAG==tag,"BA.ob"])]}

str(RW)
head(RW)

#write.csv(RW,"Harvard_2013_Diameters_reconstructed")
#View(RW)

RW$YR=as.numeric(RW$YR)
RW=subset(RW,YR>1979)

RW=read.csv("C:/R/Harvard_2013_Diameters_reconstructed.csv")


###### use young's equations

########Complete tree

###### ******QURU is not available in Youngs equations

compmass=function(SPP,DOB){
  if(SPP=='PIRU'){comp.b0=1.10651;comp.b1=2.298388}
  else if (SPP=='ABBA'){comp.b0=0.8161535;comp.b1=2.413978}
  else if(SPP=='PIST'){comp.b0=0.5724865;comp.b1=2.467798}
  else if(SPP=='TSCA'){comp.b0=0.8645225;comp.b1=2.38591}
  else if(SPP=='THOC'){comp.b0=1.32942;comp.b1=1.919051}
  else if(SPP=='BEAL'){comp.b0=1.345053;comp.b1=2.335473}
  else if(SPP=='BEPA'){comp.b0=0.74343;comp.b1=2.639837}
  else if(SPP=='ACRU'){comp.b0=1.187596;comp.b1=2.37025}
  else if(SPP=='ACSA'){comp.b0=1.477753;comp.b1=2.321565}
  else if(SPP=='FAGR'){comp.b0=1.55455;comp.b1=2.291087}
  else if(SPP=='QURU'){comp.b0=1.55455;comp.b1=2.291087}
  else {comp.b0=0;comp.b1=0}
  comp=exp(comp.b0+comp.b1*(log(DOB/2.54)))*0.4536
  return(comp)}

RW$comp=0
RW$comp=mapply(compmass,SPP=RW$SPP,DOB=RW$DOB)

##sort
RW=RW[order(RW$TAG,RW$YR),]
head(RW)

RW$TAG=as.factor(RW$TAG)

df1=data.frame(RW %>% 
                 group_by(TAG)%>%
                 mutate(compinc=c(NA,diff(comp))))
head(df1)
RW=df1

head(RW)
summary(unique(RW$TreeID))
summary(unique(RW$coreID))

RW2=RW[complete.cases(RW[,32]),]
head(RW2)
distil=aggregate(RW2$compinc,by=list(RW2$YR,RW2$Site,RW2$TreeID,RW2$SPP,RW2$DBH,RW2$Distance),FUN=mean)
colnames(distil)=c('YR','Site','TreeID','SPP','DBH','Distance','compinc')
head(distil)


summary(distil$Distance)
summary(distil$DBH)
distil$compinc=distil$compinc*1000

#plot areas
sma=pi*13^2
bia=pi*20^2



smplots=subset(distil,DBH<20.0)

colnames(smplots)
smplots=smplots[complete.cases(smplots[,7]),]
#smplots=subset(smplots,YR<2013)

ggplot(data=smplots, aes(x=YR,y=compinc,col=TreeID))+
  geom_line(aes(),size=1)+
  geom_point(aes(), shape=21,fill='black',size=2.8)+
  #geom_errorbar(aes(ymin=WBI-compSE, ymax=WBI+compSE), colour="black", width=.3,alpha=0.4) +
  theme_classic(base_size=14)+
  xlab('Year')+
  ylab(expression("Carbon mass increment  " ~ (g~C~m^{-2})))+
  #ylim(120,190)+
  #ggtitle("Standard error\n25x25 m plots") +
  #scale_x_continuous(breaks=seq(1990,2015,by=1),labels=c(1990,"",1992,"",1994,"",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"",2014,""))+
  #scale_y_continuous(breaks=seq(120,190,by=10),labels=c(120,130,140,150,160,170,180,190))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



smplotinc=aggregate(smplots$compinc,by=list(smplots$Site,smplots$YR),FUN=sum)
head(smplotinc)
colnames(smplotinc)=c('Site','YR','smcompinc')
smplotinc$incm2=smplotinc$smcompinc/sma

head(smplotinc)


biplots=subset(distil,DBH>19.9)
biplots=biplots[complete.cases(biplots[,7]),]


biplotinc=aggregate(biplots$compinc,by=list(biplots$Site,biplots$YR),FUN=sum)
head(biplotinc)
colnames(biplotinc)=c('Site','YR','bicompinc')
biplotinc$incm2=biplotinc$bicompinc/bia

head(biplotinc)

plotinc=merge(smplotinc,biplotinc,by=c("Site","YR"))
?merge

pinc=subset(plotinc,plotinc$YR>1989)
head(pinc)

pinc$incm2=(pinc$incm2.x+pinc$incm2.y)

pinc2=aggregate(pinc$incm2,by=list(pinc$Site,pinc$YR),FUN=mean)
head(pinc2)
colnames(pinc2)=c("Site","YR","compinc")

pinc2$compinc=pinc2$compinc/2


ggplot(data=pinc2, aes(x=YR,y=compinc,col=Site))+
  geom_line(aes(),size=1)+
  geom_point(aes(), shape=21,fill='black',size=2.8)+
  #geom_errorbar(aes(ymin=WBI-compSE, ymax=WBI+compSE), colour="black", width=.3,alpha=0.4) +
  theme_classic(base_size=14)+
  xlab('Year')+
  ylab(expression("Carbon mass increment  " ~ (g~C~m^{-2})))+
  #ylim(120,190)+
  #ggtitle("Standard error\n25x25 m plots") +
  #scale_x_continuous(breaks=seq(1990,2015,by=1),labels=c(1990,"",1992,"",1994,"",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"",2014,""))+
  #scale_y_continuous(breaks=seq(120,190,by=10),labels=c(120,130,140,150,160,170,180,190))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

?spread
pinc2b=spread(pinc2,Site,compinc)
head(pinc2b)

pinc3=aggregate(pinc2$compinc,by=list(pinc2$YR),FUN=mean)
pinc4=aggregate(pinc2$compinc,by=list(pinc2$YR),FUN=sd)
pinc5=cbind(pinc3,pinc4[,2])
head(pinc5)
colnames(pinc5)=c('YR','compinc','sd')

ggplot(data=pinc5, aes(x=YR,y=compinc))+
  geom_line(aes(),size=1)+
  geom_point(aes(), shape=21,fill='black',size=2.8)+
  geom_errorbar(aes(ymin=compinc-sd, ymax=compinc+sd), colour="black", width=.3,alpha=0.4) +
  theme_classic(base_size=14)+
  xlab('Year')+
  ylab(expression("Carbon mass increment  " ~ (g~C~m^{-2})))+
  #ylim(120,190)+
  #ggtitle("Standard error\n25x25 m plots") +
  #scale_x_continuous(breaks=seq(1990,2015,by=1),labels=c(1990,"",1992,"",1994,"",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"",2014,""))+
  #scale_y_continuous(breaks=seq(120,190,by=10),labels=c(120,130,140,150,160,170,180,190))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


pinc6=cbind(pinc2b,pinc5[,c(2,3)])
head(pinc6)
write.csv(pinc6,"Harvard_2013_plotlevelCMI.csv",row.names = FALSE)





ALLC=cbind(ACRU,BEAL,FAGR,QURU,TSCA)

##############################dplR

#Automatically detrend all
PIRUd=detrend(PIRU,method='Spline')

#Detrend one by one (recommended)
PIRUd=i.detrend(PIRU)
head(PIRUd)
PIRUd$YR=row.names(PIRUd)

###Here I saved the detrended ring widths for a different analysis, probably not needed for you
PIRUd2=gather(PIRUd,treeID,stdRW,-YR,na.rm=T)
head(PIRUd2)
write.csv(PIRUd2,'piruD.csv')

###Use detrended series to make chronology
PIRUc=chron(PIRUd)
head(PIRUc)
summary(PIRUc)
plot(PIRUc)

###Save chronology for later use
write.csv(PIRUc,'chron_piru_RW_R.csv')

###fix formatting when you read in the .csv for later use
PIRUc=read.csv('chron_piru_RW_R.csv')
rownames(PIRUc)=PIRUc$X
PIRUc=PIRUc[,-1]






######################################treeclim 

#spruce vs. tmax from previous april to current october, can choose correlation or response function in 'method'
PIRUtcc=dcc(PIRUc,tmax,selection = -4:10,method='correlation')
plot(PIRUtcc)

#spruce vs. precip
PIRUpcc=dcc(PIRUc,precip,selection = -4:10,method='correlation')
plot(PIRUpcc)






#######################################ggplot2


####Nicer looking graphs from treeclim output

PIRUctcc=dcc(PIRUc,tmax,selection = -4:10,method='correlation')
PIRUcbp=data.frame(PIRUctcc$coef$month,PIRUctcc$coef$coef,PIRUctcc$coef$significant)
colnames(PIRUcbp)=c('month','coef','significant')
PIRUcbp$mon=1:19
head(PIRUcbp)
gpiru1<-ggplot(data=PIRUcbp, aes(x=mon,y=coef,fill=significant))+
  geom_bar(stat='identity')+
  theme_classic(base_size=12)+
  #xlab('Month')+
  ylab('Coefficients\n')+
  ggtitle('Red spruce : Average daily maximum temperature')+
  scale_fill_manual(values=c("lightsteelblue","red3"))+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(breaks=c(1:19),labels=PIRUcbp$month)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gpiru1


PIRUcpcc=dcc(PIRUc,precip,selection = -4:10,method='correlation')
PIRUcbp=data.frame(PIRUcpcc$coef$month,PIRUcpcc$coef$coef,PIRUcpcc$coef$significant)
colnames(PIRUcbp)=c('month','coef','significant')
PIRUcbp$mon=1:19
head(PIRUcbp)
gpiru4<-ggplot(data=PIRUcbp, aes(x=mon,y=coef,fill=significant))+
  geom_bar(stat='identity')+
  theme_classic(base_size=12)+
  #xlab('Month')+
  ylab('Coefficients\n')+
  ggtitle('Red Spruce : Monthly total precipiation')+
  scale_fill_manual(values=c("lightsteelblue","red3"))+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(breaks=c(1:19),labels=PIRUcbp$month)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gpiru4




