# Course Bayesian thinking - UZH 2021
# Author: Gesa von Hirschheydt


#  Add working path  ####
#-----------------------#

# This code needs to be adjusted

# set working directory
setwd("H:/GitHub/Git_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/0_Code")

# define path to data folder
folder <- "H:/GitHub/Git_repositories_to_pull_from/Bayesian_thinking_UZH/Group_projects/SDM_swiss_birds/1_Data/"



#  Load data  ####
#----------------#

bb <- read.table(paste(folder,"EBBA2subset.csv",sep=""), header=TRUE, sep =";")
str(bb)
head(bb)

# load packages
library(reshape2)  # for function reshape()



#  SURVEYS  ####
#--------------#

# table which contains all visit-specific data: site, date, survey number, duration of visit (minutes)

bb$Date.of.survey <- as.Date(bb$Date.of.survey, format="%d.%m.%Y")

surveys <- unique(bb[,c("Square.code","Date.of.survey","Survey.number","Duration.of.survey..minutes.")])
(t <- table(surveys$Square.code))
colnames(surveys) <- c("site","date","survey","time")

# export cleaned data
# write.table(surveys, file=paste(folder,"surveys.csv",sep=""), sep=";", row.names=FALSE)



#  SITES  ####
#------------#

# table which contains all site-specific data: site, coordinates, environmental variables

sites.help <- unique(bb[,c("Square.code","x_lv95","y_lv95","kmx_lv03","kmy_lv03","Latitude","Longitude")])
colnames(sites.help)[1] <- "site"
str(sites.help)
sites.help$id <- paste(sites.help$x_lv95-2000000, "_", sites.help$y_lv95-999999, sep=""); length(unique(sites.help$id))

# combine with Umweltdaten
umwelt <- read.delim(paste(folder,"Umweltdaten_use.txt",sep=""), header=TRUE)
umwelt$id <- paste(umwelt$x, "_", umwelt$y, sep="")
length(unique(umwelt$id)) # = all km2 in Switzerland
sum(unique(sites.help$id) %in% unique(umwelt$id)) # for check

sites <- merge(sites.help, umwelt, by="id", all.x=TRUE)
str(sites); head(sites)
paste(colnames(sites), collapse="','")
sites <- sites[,c("site","buildings","elevation","farmland","forest","glacier",
                  "grassland","lakes","northness","rivers","roads","rocks",
                  "shoreline","slope","structures","wetlands","x_lv95","y_lv95",
                  "kmx_lv03","kmy_lv03","Latitude","Longitude","x","y","id")]
head(sites)

# export cleaned data
# write.table(sites, file=paste(folder,"sites.csv",sep=""), sep=";", row.names=FALSE)



#  KITE  ####
#-----------#

kite.help1 <- subset(bb, EBBA2.scientific_name=="Milvus migrans")
str(kite.help1)
colnames(kite.help1)[1] <- "survey"
kite.help1$det <- 1
kite.help1 <- kite.help1[,c("survey","det")]
head(kite.help1)

# merge with surveys to fill in non-detections
kite.help2 <- merge(kite.help1, surveys, by="survey", all=TRUE)
kite.help2 <- kite.help2[order(kite.help2$site),c("site","date","det","time","survey")]
kite.help2[is.na(kite.help2$det),"det"] <- 0
head(kite.help2); str(kite.help2)

# Number visits through from 1:max.number.of.visits.per.site
kite.help2$visit <- NA
for(i in 1:nlevels(as.factor(kite.help2$site))){
  help <- which(kite.help2$site==levels(as.factor(kite.help2$site))[i])
  kite.help2[help,"visit"] <- 1:length(help)
}

# merge with sites to fill in environmental variables
kite <- merge(kite.help2, sites, by="site", all=TRUE)

str(kite)
head(kite[,1:21])

# export cleaned data
# write.table(kite, file=paste(folder,"kite.csv",sep=""), sep=";", row.names=FALSE)



#  KITE_FOR_JAGS  ####
#--------------------#

str(kite); head(kite)
help1 <- kite[order(kite$visit),c("det","visit","site")]
help2 <- reshape(help1,
                v.names="det",
                timevar="visit",
                idvar="site",
                direction="wide")
head(help2,10)
rownames(help2) <- help2[,1]; help2 <- help2[,-1]; colnames(help2) <- paste("visit",1:9,sep="")
kite_for_jags <- help2; dim(kite_for_jags); head(kite_for_jags)


# export cleaned data
# write.table(kite_for_jags, file=paste(folder,"kite_for_jags.csv",sep=""), sep=";", row.names=FALSE)
