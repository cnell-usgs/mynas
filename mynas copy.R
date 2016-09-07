##myna point count data
mynas<-read.csv("/Users/colleennell/Documents/R/mynas/mynas.csv")
View(mynas)
counts<-read.csv("/Users/colleennell/Documents/R/mynas/mynas_counts.csv")
View(counts)

library(unmarked)
library(dplyr)

##make values for all COMY, natives, JUMY?
counts$COMY.10[is.na(counts$COMY.10)] <- 0
counts$X25[is.na(counts$X25)] <- 0
counts$X50[is.na(counts$X50)] <- 0
counts$all_COMY<-counts$COMY.10+counts$X25+counts$X50

envdata<- left_join(counts,mynas,by=c("Transect","Point"))
View(envdata)
#estimate
#point rtansect surveys
comy<-subset(envdata,select=c(Date,Time,Transect,Point,Obs,COMY.10,X25,X50))
#no counts of COMY were made after line 676 
comy<-subset(comy,Date!="")
#distance intervals
breaks<-c(0,10,25,50)
n<-length(breaks)-1

y<-comy[,6:8]
y<-as.matrix(y)



####unmarked
obvars<-c("dist.forest","habitat","Develop","Tree.cover")
covars<-c("Transect","Point")
birds<-c("all_COMY","total_native")
str(envdata)
umdf<-select_(envdata,"all_COMY","total.native.1","dist.forest","area.1","Develop","Tree.cover","Transect","Point","Obs")
View(umdf)
umdm<-as.matrix(umdf)

countcov<- unmarkedFramePCount(
  y= umdm[,1:2],
  siteCovs = umdm[,7:9],
  obsCovs = list(enviro = umdm[,3:6])
)
