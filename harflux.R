#### CO2 Flux adjustments
library(ggplot2)
library(lubridate)
library(date)
library(Hmisc)
setwd("C:/R")
library(plotrix)


#####################Data import
##Main Tower
harvf <- read.delim("HF_9215_filled.txt")
head(harvf)

harf=harvf[,c(1,2,3,4,11)]
head(harf)
colnames(harf)=c('YR','month','doy','hr','NEE')
flux=harf
str(flux$NEE)
YR<-as.factor(flux$YR)
flux$NEE=as.numeric(flux$NEE)
YRflux=aggregate(flux$NEE, by=list(flux$YR), FUN=sum)
head(YRflux)
colnames(YRflux)=c("YEAR","FLUX")
YRflux$FLUX=(abs(YRflux$FLUX)/10^6)*12.01*1800
year<-as.numeric(YRflux$YEAR)
plot(YRflux$FLUX~year,type="b")
head(YRflux)
colnames(YRflux)=c('YR','calendarflux')

###find optimal

head(flux)
flux=subset(flux,flux$YR>1991)


flux$f2YR<-ifelse(flux$doy>161,(flux$YR+1),(flux$YR))
head(flux)
summary(flux$YR)
f2YRflux=aggregate(flux$NEE, by=list(flux$f2YR), FUN=sum)
colnames(f2YRflux)=c("YR","f2FLUX")
f2YRflux$f2Cflux=(abs(f2YRflux$f2FLUX)/10^6)*12.01*1800
head(f2YRflux)
plot(f2YRflux$f2Cflux~f2YRflux$YR,type="b")

flux2=merge(f2YRflux,YRflux,by='YR')

inc=read.csv('Harvard_2013_plotlevelCMI.csv')
head(inc)
yrINC=inc
plot(yrINC$compinc~yrINC$YR,type='b')


d5=merge(yrINC,flux2,by='YR')
d5=subset(d5,d5$YR<2013)
head(d5)
d5=d5[-1,]
cor(d5$calendarflux,d5$compinc)

cor(d5$f2Cflux,d5$compinc)


####plot the relationship before finding optimal 

d6=d5
head(d6)
## Plot first set of data and draw its axis
plot(d6$YR, d6$compinc, pch=21, xlab="Year", ylab=expression("Woody biomass increment  " ~ (g~C~m^{-2})),type="o",
     col="black",cex.lab=1.25,cex.axis=1.1,cex=1.2,xaxt = "n",bg="black")
labels=c(1990,"",1992,"",1994,"",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"")
#text(1996:2015, par("usr")[3], srt = 45,adj=c(1.7,0.3),labels = labels, xpd = TRUE)
## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(d6$YR, d6$calendarflux, pch=24,bg='red',  xlab="",ylab="", 
     axes=FALSE, type="b", cex=1.2,col="red",cex.lab=1.25,cex.axis=1.25,lty=1,lwd=1.5)
## a little farther out (line=4) to make room for labels
mtext(expression("Net ecosystem exchange  " ~ (g~C~m^{-2})),side=4,col="red",line=4,cex=1.2) 
axis(4, ylim=c(1380,3550), col="red",col.axis="red",cex=1.25)
#abline(v=2007.5,lty=2)
x=1990:2013
axis(side = 1, at = x, labels = labels, tck = -0.02)
#arrows(c(2007,2008),c(355.0,355.0),c(2002,2013),c(355.0,355.0), lwd=2)
#text(2006.5,344,'(a)',cex=1.25)
#text(2008.5,344,'(b)',cex=1.25)

## Add Legend
legend("topleft",legend=c("Woody biomass ",expression("NEE")),
       text.col=c("black","red"),pch=c(19,17),col=c("black","red"),bg=c("black","red"),cex=1.2,bty='n')


head(flux)

## Plot first set of data and draw its axis
plot(d6$YR, d6$compinc, pch=21, xlab="Year", ylab=expression("Woody biomass increment  " ~ (g~C~m^{-2})),type="o",
     col="black",cex.lab=1.25,cex.axis=1.1,cex=1.2,xaxt = "n",bg="black")
labels=c(1990,"",1992,"",1994,"",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"")
#text(1996:2015, par("usr")[3], srt = 45,adj=c(1.7,0.3),labels = labels, xpd = TRUE)
## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(d6$YR, d6$f2Cflux, pch=24,bg='red',  xlab="",ylab="", 
     axes=FALSE, type="b", cex=1.2,col="red",cex.lab=1.25,cex.axis=1.25,lty=1,lwd=1.5)
## a little farther out (line=4) to make room for labels
mtext(expression("Net ecosystem exchange  " ~ (g~C~m^{-2})),side=4,col="red",line=4,cex=1.2) 
axis(4, ylim=c(1380,3550), col="red",col.axis="red",cex=1.25)
#abline(v=2007.5,lty=2)
x=1990:2013
axis(side = 1, at = x, labels = labels, tck = -0.02)
#arrows(c(2007,2008),c(355.0,355.0),c(2002,2013),c(355.0,355.0), lwd=2)
#text(2006.5,344,'(a)',cex=1.25)
#text(2008.5,344,'(b)',cex=1.25)

## Add Legend
legend("topleft",legend=c("Woody biomass ",expression("NEE")),
       text.col=c("black","red"),pch=c(19,17),col=c("black","red"),bg=c("black","red"),cex=1.2,bty='n')






m1=data.frame(doy=1:365)

for (i in 1:365){
  flux$f2YR<-ifelse(flux$doy>i,(flux$YR+1),(flux$YR))
  f2YRflux=aggregate(flux$NEE, by=list(flux$f2YR), FUN=sum)
  colnames(f2YRflux)=c("YR","f2FLUX")
  f2YRflux$f2Cflux=(abs(f2YRflux$f2FLUX)/10^6)*12.01*1800
  d5=merge(yrINC,f2YRflux,by='YR')
  d5=subset(d5,d5$YR>1990)
  d5=subset(d5,d5$YR<2013)
  m1$a[i]<-cor(d5$f2Cflux,d5$compinc)
}


head(m1)
m1$shft_days=m1$doy-366
plot(m1$a~m1$shft_days,type='l')


m1
#write.csv(m1,'Bartlett_shift_20180917.csv')


####falls apart here


doy=read.csv('Bartlett_shift_20180917.csv')
head(doy)
shift=read.csv('Bartlett_shift_20180917.csv')
head(shift)
shift2=subset(shift,shift$shft_days>-365)
shift2$length='day'
shift3=subset(shift2,shift2$jday %in% c('365','335','305','274','244','213','182','152','121','91','60'))
shift3$length='month'
shift2=rbind(shift2,shift3)
head(shift2)

pl3=ggplot(data=shift2,aes(x=shft_days,y=value))+
  geom_line(data=subset(shift2,length=='day'),aes(y=corr,color='DAY'),size=0.8,color='black')+
  geom_vline(xintercept=-121,linetype=2)+
  geom_point(data=subset(shift2,length=='month'),aes(y=corr,fill='Beginning of month'),shape=18,size=6,alpha=0.5)+
  theme_classic(base_size=18)+
  ylim(0.35,0.735)+
  xlab('Start of annual flux summary')+
  ylab('Correlation (r)\n')+
  #geom_vline(xintercept=-131,linetype=2)+
  #geom_vline(xintercept=-160,linetype=2)+
  scale_x_continuous(breaks=doy$shift,labels=doy$name)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=13))+
  theme(legend.justification=c(-0.05,1), legend.position=c(0,.99))+
  #guides(col = guide_legend(ncol = 3))+
  theme(legend.title=element_blank(),legend.key = element_rect(fill = "transparent", colour = "transparent"))+
  #guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

pl3



comp=d5
head(comp)
comp$YEAR=comp$YR
comp$opt=comp$f2Cflux

plot(comp$opt~comp$compinc)

par(mar=c(5, 5, 2, 5) + 0.1)

head(comp)
## Plot first set of data and draw its axis
plot(comp$YEAR, comp$compinc/2.96/10, pch=21, xlab="Year", ylab=expression("Woody biomass increment  " ~ (g~C~m^{-2})),type="o",
     col="black",cex.lab=1.25,cex.axis=1.1,cex=1.2,xaxt = "n",bg="black")
#labels=c("",1996,"",1998,"",2000,"",2002,"",2004,"",2006,"",2008,"",2010,"",2012,"",2014,"")
#text(1996:2015, par("usr")[3], srt = 45,adj=c(1.7,0.3),labels = labels, xpd = TRUE)
## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(comp$YEAR, comp$opt, pch=24,bg='red',  xlab="",ylab="", 
     axes=FALSE, type="b", cex=1.2,col="red",cex.lab=1.25,cex.axis=1.25,lty=1,lwd=1.5)
## a little farther out (line=4) to make room for labels
mtext(expression("Net ecosystem exchange  " ~ (g~C~m^{-2})),side=4,col="red",line=4,cex=1.2) 
axis(4, ylim=c(1380,3550), col="red",col.axis="red",cex=1.25)
#abline(v=2007.5,lty=2)
x=1992:2013
axis(side = 1, at = x, labels = labels, tck = -0.02)
#arrows(c(2007,2008),c(355.0,355.0),c(2002,2013),c(355.0,355.0), lwd=2)
#text(2006.5,344,'(a)',cex=1.25)
#text(2008.5,344,'(b)',cex=1.25)

## Add Legend
legend("topleft",legend=c("Woody biomass ",expression("NEE")),
       text.col=c("black","red"),pch=c(19,17),col=c("black","red"),bg=c("black","red"),cex=1.2,bty='n')





