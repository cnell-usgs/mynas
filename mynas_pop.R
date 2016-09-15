##estimating myna poplation density and identifying mediating factors 
#using unmarked package in R

#http://www.montana.edu/screel/Webpages/conservation%20biology/BIOE%20440%20distance%20sampling.html

#data
transects<-read.csv("/Users/colleennell/Documents/R/mynas/mynas.csv") ##env data on transect points
counts<-read.csv("/Users/colleennell/Documents/R/mynas/mynas_counts.csv") #point-transect data counts for mynas annd al species within 10, 25, and 50m of point
View(transects)
View(counts)

#packages
library(unmarked)
library(dplyr)
library(reshape2)

##data prep
counts$IDENT<-paste0(counts$Transect,"0",counts$Point) #create unique ID for each transect point that corresponds to GPS pts
myna.df<-counts[1:675,c("Date","Time","Transect","Point","Obs","COMY.10","X25","X50")]
myna.df[is.na(myna.df)]<-0
View(myna.df) #3 separate observations of mynas at distance intervals for each pt

##distance sampling models correct density estimates for undetected animals by using the distribution of distances
#distsamp
#allows estimation of covariates for detection and density
#controllling for fx of detection is necessary for accurate estimation of density
#covariate of density is the primary interest for interpretation

#set up dataframe for unmarked
#data in long format to include all species, cast by distances, then subset for analyses
counts.obs<-counts[1:675,1:40]#ust the bird obs
counts.obs$IDENT<-paste0(counts.obs$Transect,"0",counts.obs$Point)
count.melt<-melt(counts.obs, id.vars=c("Date","Time","IDENT","Transect","Point","Obs"),,variable.name="distance", na.rm=FALSE, value.name="abundance")
count.melt$value[is.na(count.melt$value)]<-0 #make NA's 0
count.melt<-count.melt[!grepl("total",count.melt$variable),]##remove 'total' variables
##add a column for species ID, and distances
count.melt$species<-as.factor(ifelse(count.melt$variable == "COMY.10" | count.melt$variable == "X25" | count.melt$variable == "X50", "COMY",
                           ifelse(count.melt$variable == "JUMY" | count.melt$variable == "X25.1" | count.melt$variable == "X50.1", "JUMY",
                                  ifelse(count.melt$variable == "RVBU.10" | count.melt$variable == "X25.2" | count.melt$variable == "X50.2", "RVBU",
                                         ifelse(count.melt$variable == "WAHO" | count.melt$variable == "X25.3" | count.melt$variable == "X50.3", "WAHO",
                                                ifelse(count.melt$variable == "COLK" | count.melt$variable == "X25.4" | count.melt$variable == "X50.4", "COLK",
                                                       ifelse(count.melt$variable == "SAST" | count.melt$variable == "X25.5" | count.melt$variable == "X50.5", "SAST",
                                                              ifelse(count.melt$variable == "POST" | count.melt$variable == "X25.6" | count.melt$variable == "X50.6", "POST",
                                                                     ifelse(count.melt$variable == "CAHO" | count.melt$variable == "X25.7" | count.melt$variable == "X50.7", "CAHO",
                                                                            ifelse(count.melt$variable == "PCFD" | count.melt$variable == "X25.8" | count.melt$variable == "X50.8", "PCFD",NA
                                                                            ))))))))))
                           
count.melt$dist<-as.factor(ifelse(count.melt$variable == "COMY.10" | count.melt$variable == "JUMY" | count.melt$variable == "RVBU.10" | count.melt$variable == "WAHO" | count.melt$variable == "COLK" | count.melt$variable == "SAST" | count.melt$variable == "POST" | count.melt$variable == "CAHO" | count.melt$variable == "PCFD", "10",
                        ifelse(count.melt$variable == "X25" | count.melt$variable == "X25.1" | count.melt$variable == "X25.2" | count.melt$variable == "X25.3" | count.melt$variable == "X25.4" | count.melt$variable == "X25.5" | count.melt$variable == "X25.6" | count.melt$variable == "X25.7" | count.melt$variable == "X25.8", "25",
                               ifelse(count.melt$variable == "X50" | count.melt$variable == "X50.1" | count.melt$variable == "X50.2" | count.melt$variable == "X50.3" | count.melt$variable == "X50.4" | count.melt$variable == "X50.5" | count.melt$variable == "X50.6" | count.melt$variable == "X50.7" | count.melt$variable == "X50.8", "50", 
                                      NA)))) ##distance variable

##recast so distances are columns
drop<-c("variable")
count.melt<-count.melt[,!(names(count.melt) %in% drop)]#get rid of 'variable' now that it has been split into 2 vars
count.cast<-dcast(count.melt,Date+Time+IDENT+Transect+Point+Obs+species~dist,value.var="value")
View(count.cast) #now a beautiful df that cant be filtered by species 
write.csv(count.cast,file="~/mynas_dist_df.csv")

mynas<-filter(count.cast, species == "COMY")
str(mynas)
mynas$Point<-as.ordered(mynas$Point)

##unmakrd df
count.melt$dist<-as.character(count.melt$dist)
count.melt$dist<-as.numeric(count.melt$dist)
View(count.melt)
yDat<-formatDistData(count.melt,distCol="dist",transectNameCol="Transect",dist.breaks=c(0,10,25,50)) #this requires dist is numeric which it isnt

##cont
yDat<-mynas[,c(8,9,10)]
covs<-mynas[,c(1,2,4,5,6)]
umf<-unmarkedFrameDS(y=as.matrix(yDat),siteCovs=covs,survey="point",dist.breaks=c(0,10,25,50),unitsIn="m")

hist(umf, xlab="distance(m)", main=" ", cex.lab=0.8, cex.axis=0.8)#freq dist

##fitting distance models
m.half<-distsamp(~1 ~1, umf, keyfun="halfnorm", output="density", unitsOut="ha")
m.half #overall detection and density estimate proabilities
##fit appropriate distribution to distance histogram- halfnorm, hazard, uniform
m.haz<-distsamp(~1 ~1, umf, keyfun="hazard", output="density", unitsOut="ha")##lowest AIC
m.uni<-distsamp(~1 ~1, umf, keyfun="uniform", output="density", unitsOut="ha")
m.haz
m.uni

##fit a priori model set of covariates
#goodness of fit test

##backtransform
backTransform(m.haz,type="state")
backTransform(m.haz,type="det")
