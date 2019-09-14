#Contact: Vaibhav Menon (vmeno002@ucr.edu)
#The purpose of this script is to generate a code that quantifies feeding events and licking events. The
#skeleton of the workhorse function here, feedslicksplots, is a heavily adapted version of Edric Tam's
#makePlotsandTables.R large function, used in the Elife paper Joseph et al, 2017. 
#The basis for separation of feeding from licking is through an intensity threshold (feeding > 100 A.U. 
#above baseline, tasting < 100 A.U. above baseline) that is consistent with other papers using FLIC assays.
#This code is clunky and can definitely be implemented more efficiently, but extracted the information we were interested in.

library(dplyr)
library(gtools)
library(MESS)
library(openxlsx)
#modified function makeplotsandtables
df3 <- read.csv("DFM_3.csv")
#This function is called feedslicksplots because it generates a list of candidate feed/lick events and peak plots for each well.
feedslicksplots<-function(df) { 
  baseline.offset.for.contact.threshold <- 20
  begin.time = 1500
  end.time = df$Sample[length(df$Sample)-1] #refers to sample column in raw data
  # create the raw data trace plots 
  pdf('peak_plots.pdf')
  durationList = list()
  # loop over each of the 12 wells
  for (i in 1:12) {
    # +4 due to the df having time, samples and other columns (totalling 4 columns) before the 12 actual data columns
    datum<-df[begin.time:end.time,i+4]
    # set threshold for sip/contact/licks. mode is baseline
    data.baseline <- find.mode(datum)
    contactThreshold <- data.baseline + baseline.offset.for.contact.threshold #mode + 20 = is threshold for contact
    # define licks as data that is above contactThreshold
    if(length(which(datum>contactThreshold)) >= 0) {
      lickintsfile <- as.character(paste("lickints",i,".csv",sep=""))
      feeds <- which((datum>=contactThreshold)) #feeds contains andidate Events (not yet classified as feeding or tasting)
      feedints <- data.frame(feeds, datum[feeds]) 
      lickints <- data.frame(feeds, datum[feeds]) #create variables to store events when tehy are classified
      colnames(lickints) <- c("licks","datum.licks.")
      feedaucs <- auc(feedints[,1],feedints[,2]) #area under curve
      lickaucs <- auc(lickints[,1],lickints[,2])
      feedaucsglob <- as.character(paste("feedaucs",i,sep=""))
      lickaucsglob <- as.character(paste("lickaucs",i,sep=""))
      modes <- as.character(paste("mode",i,sep=""))
      assign(feedaucsglob,feedaucs,envir=globalenv()) #creating global variables for AUC
      assign(lickaucsglob,lickaucs,envir=globalenv())
      assign(modes, data.baseline, envir=globalenv())
      feedintsfile <- as.character(paste("feedints",i,".csv",sep=""))
      write.csv(lickints, lickintsfile) 
      write.csv(feedints, feedintsfile) #write these 2col dataframes to current directory
    } else {
      licks <-c(0)
    }
   #Generates peak plots
    plot(seq_along(as.vector(datum)), as.vector(datum), type="l",main=colnames(df)[i+4], xlab="Time Indices (5/s)", ylab="Peak Intensity")#ylim=c(min(peakIntensity),max(peakIntensity)))
    
  }
  dev.off()
}

makePlotsAndTables(df3)
#source('run.R')
source('flic loading organization.R')
#create function to group events based on time

#read lick and feed csvs
for (i in 1:12){
  lickfile <- as.character(paste("lickints", i, ".csv", sep=""))
  lickreadfile <- as.character(paste("licks",i,sep=""))
  assign(lickreadfile, read.csv(lickfile))
  feedfile <- as.character(paste("feedints",i,".csv",sep=""))
  feedreadfile <- as.character(paste("feeds",i,sep=""))
  assign(feedreadfile, read.csv(feedfile))}

#annotating from description text file. NOTE: description is from flic loading organization.R, the numbers after it in 
#str_sub are specific to the directory which contains information I want to annotate the well
#if it's another directory, adjust numbers
runinfo <- str_sub(description, 57,66)
runinfo <- paste(runinfo, " DFM3",sep="")
#WELL 1 CHUNK
#grouping Feeding events- Well 1
j <- 1
intervalF1 <- NULL
if (length(feeds1$feeds)>0){
feeds1$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds1$feeds)-1)){
  if (length(feeds1$feeds)==0){
    break
  }
  else if ((feeds1[i+1,2])-(feeds1[i,2])<5){
    feeds1$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds1$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
}
aucfeed1 <- split(feeds1,feeds1$EventNumber)
aucfeed1 <- aucfeed1[sapply(aucfeed1, function(d) max(d[[3]]) > mode1+100)] #max intensity of feeding event >100
feeds1 <- do.call("rbind",aucfeed1)
feeds1 <- feeds1[order(feeds1$feeds),]
j <- 1
intervalF1 <- NULL
if (length(feeds1$feeds)>0){
  feeds1$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds1$feeds)-1)){
  if (length(feeds1$feeds)==0){
    break
  }
  else if ((feeds1[i+1,2])-(feeds1[i,2])<5){
    feeds1$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds1$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds1$feeds[i+1]-feeds1$feeds[i]
    intervalF1 <- rbind(intervalF1,interval)
  }
}


##feeding events may start below feeding threshold (100 above baseline), and cross threshold later in event

if (length(feeds1$feeds)>0){
feeds1 <- group_by(feeds1, EventNumber) %>% filter(n()>1) #get rid of events only 1 second long 
#finding max intensity in each feeding event
aveF1 <-aggregate(feeds1[, 3], list(feeds1$EventNumber), max)
aveF1 <-cbind(aveF1, aggregate(feeds1[, 3], list(feeds1$EventNumber), length)) 
aveF1 <- aveF1[,-3]
w1f1 <-feeds1[feeds1$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF1 <- cbind(aveF1, mean(aveF1[,3]), mean(feeds1$datum.feeds.),sum(aveF1[,3]),mean(intervalF1),
               auc((1:length(feeds1$feeds)),feeds1$datum.feeds.),genotype1,taste5)
colnames(aveF1) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
aveF1 <- aveF1[mixedorder(aveF1$`Event Number`),]}
if (length(feeds1$feeds)==0){
  aveF1 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF1) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)","Genotype", "Tastant")
}
#grouping Licking events- Well 1
j <- 1
intervalL1 <- NULL
if (length(licks1$licks)>0){
licks1$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks1$licks)-1)){
  if (length(licks1$licks)==0){
    break
  }
  else if ((licks1[i+1,2])-(licks1[i,2])<5){ #if two consecutive tracked lick events separated by 5 units of time aka 1s (1 unit =.2s), break into 2 events 
    licks1$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks1$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
}
auclick1 <- split(licks1,licks1$EventNumber)
auclick1 <- auclick1[sapply(auclick1, function(d) max(d[[3]]) < mode1+100)] #max intensity of lick events <100 over baseline
licks1 <- do.call("rbind",auclick1)
licks1 <- licks1[order(licks1$licks),]

j <- 1
intervalL1 <- NULL
if (length(licks1$licks)>0){
  licks1$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks1$licks)-1)){
  if (length(licks1$licks)==0){
    break
  }
  else if ((licks1[i+1,2])-(licks1[i,2])<5){
    licks1$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks1$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks1$licks[i+1]-licks1$licks[i]
    intervalL1 <- rbind(intervalL1,interval)
  }
}
#get rid of events only 1 time unit long
if(length(licks1$licks)>0){
licks1 <- group_by(licks1, EventNumber) %>% filter(n()>1)
#finding max intensity of each tasting event
aveL1 <-aggregate(licks1[, 3], list(licks1$EventNumber), max) 
aveL1 <-cbind(aveL1, aggregate(licks1[, 3], list(licks1$EventNumber), length)) 
aveL1 <- aveL1[,-3]
aveL1 <- cbind(aveL1, mean(aveL1[,3]), mean(licks1$datum.licks.),sum(aveL1[,3]),mean(intervalL1),
               auc((1:length(licks1$licks)),licks1$datum.licks.),genotype1, taste5)
colnames(aveL1) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype", "Tastant")
aveL1 <- aveL1[mixedorder(aveL1$`Event Number`),]}
if (length(licks1$licks)==0){
  aveL1 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL1) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events","Lick AUC (a.u^2)", "Genotype", "Tastant")
}
if (is.null(intervalF1)==TRUE){
  intervalF1 <- NA
}
if (length(feeds1)>0){ #create summary stats for this well (corresponds to Well 1 in Plate Stats.csv)
  aveW1 <- cbind(aveF1$`Mean Length of all Events`[1],aveL1$`Mean Length of all Events`[1],mean(aveF1$`Mean Intensity(per event)`),
                 mean(aveL1$`Mean Intensity(per event)`),aveF1$`Total time (All Events)`[1],aveL1$`Total time (All Events)`[1],
                 aveF1$`Mean interval between Events`[1],aveL1$`Mean interval between Events`[1], length(aveF1$`Event Number`),
                 length(aveL1$`Event Number`),
                 aveF1$`Feed AUC (a.u^2)`, aveL1$`Lick AUC (a.u^2)`,max(w1f1$datum.feeds.),length(w1f1$datum.feeds.), 
                 auc((1:length(w1f1$feeds)),w1f1$datum.feeds.),intervalF1[1], feeds1$feeds[1],
                 NA,genotype1, taste5, runinfo, 1)
  aucvec1 <- lapply(aucfeed1, function(d) (auc(d[[2]],d[[3]])))
  aucvec1 <- unlist(aucvec1)
  aucvec1 <- as.vector(aucvec1)
  deltaintw1 <- data.frame(aveF1$`Event Number`, ((aveF1$`Mean Intensity(per event)`
                                                   -mean(aveF1$`Mean Intensity(per event)`))/mean(aveF1$`Mean Intensity(per event)`)),((aveF1$`Length of Event`-mean(aveF1$`Length of Event`))/
                                                                                                                                         mean(aveF1$`Length of Event`)),(aucvec1-mean(aucvec1))/mean(aucvec1),
                           ((intervalF1[1]-mean(intervalF1))/mean(intervalF1)),((intervalF1[2]-mean(intervalF1))/mean(intervalF1)),((intervalF1[3]-mean(intervalF1))/mean(intervalF1)),
                           ((intervalF1[4]-mean(intervalF1))/mean(intervalF1)),((intervalF1[5]-mean(intervalF1))/mean(intervalF1)),
                           genotype1, taste5,paste(runinfo,"Well 1", sep=" "))}
if (length(feeds1)==0){
  aveW1 <- cbind(aveF1$`Mean Length of all Events`[1],aveL1$`Mean Length of all Events`[1],aveF1$`Mean Intensity (all Events)`[1],
                 aveL1$`Mean Intensity (all Events)`[1],aveF1$`Total time (All Events)`[1],aveL1$`Total time (All Events)`[1],
                 aveF1$`Mean interval between Events`[1],aveL1$`Mean interval between Events`[1], length(aveF1$`Event Number`),
                 length(aveL1$`Event Number`),
                 aveF1$`Feed AUC (a.u^2)`, aveL1$`Lick AUC (a.u^2)`,0,0, 
                 0,intervalF1[1], NA, NA,genotype1, taste5, runinfo, 1)
  deltaintw1 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype1,taste5,paste(runinfo,"Well 1", sep=" "))}
colnames(deltaintw1) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW1) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
#NOTE: Well 2-12 are replicates of the code for Well 1, save for appropriate tastant and genotype annotated information.
###WELL 2 CHUNK
#grouping feeding events- Well2
j <- 1
intervalF2 <- NULL
if (length(feeds2$feeds)>0){
feeds2$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds2$feeds)-1)){
  if (length(feeds2$feeds)==0){
    break
  }
  else if ((feeds2[i+1,2])-(feeds2[i,2])<5){
    feeds2$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds2$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds2$feeds[i+1]-feeds2$feeds[i]
    intervalF2 <- rbind(intervalF2,interval)
  }
}
aucfeed2 <- split(feeds2,feeds2$EventNumber)
aucfeed2 <- aucfeed2[sapply(aucfeed2, function(d) max(d[[3]]) >mode2+100)]
feeds2 <- do.call("rbind",aucfeed2)
feeds2 <- feeds2[order(feeds2$feeds),]
j <- 1
intervalF2 <- NULL
if (length(feeds2$feeds)>0){
  feeds2$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds2$feeds)-1)){
  if (length(feeds2$feeds)==0){
    break
  }
  else if ((feeds2[i+1,2])-(feeds2[i,2])<5){
    feeds2$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds2$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds2$feeds[i+1]-feeds2$feeds[i]
    intervalF2 <- rbind(intervalF2,interval)
  }
}
#get rid of events only 1 second long
if (length(feeds2$feeds)){
feeds2 <- group_by(feeds2, EventNumber) %>% filter(n()>1)
#calculating mean
aveF2 <-aggregate(feeds2[, 3], list(feeds2$EventNumber), max) 
aveF2 <-cbind(aveF2, aggregate(feeds2[, 3], list(feeds2$EventNumber), length)) 
aveF2 <- aveF2[,-3]
w2f1 <- feeds2[feeds2$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF2 <- cbind(aveF2, mean(aveF2[,3]), mean(feeds2$datum.feeds.),sum(aveF2[,3]),mean(intervalF2), 
               auc((1:length(feeds2$feeds)),feeds2$datum.feeds.),genotype1, taste6)
colnames(aveF2) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
aveF2 <- aveF2[mixedorder(aveF2$`Event Number`),]}
if (length(feeds2$feeds)==0){
  aveF2 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF2) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
}
#grouping Licking events- Well 2
j <- 1
intervalL2 <- NULL
if (length(licks2$licks)>0){
licks2$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks2$licks)-1)){
  if (length(licks2$licks)==0){
    break
  }
  else if ((licks2[i+1,2])-(licks2[i,2])<5){
    licks2$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks2$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks2$licks[i+1]-licks2$licks[i]
    intervalL2 <- rbind(intervalL2,interval)
  }
}
auclick2 <- split(licks2,licks2$EventNumber)
auclick2 <- auclick2[sapply(auclick2, function(d) max(d[[3]]) < mode2+100)]
licks2 <- do.call("rbind",auclick2)
licks2 <- licks2[order(licks2$licks),]
j <- 1
intervalL2 <- NULL
if (length(licks2$licks)>0){
  licks2$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks2$licks)-1)){
  if (length(licks2$licks)==0){
    break
  }
  else if ((licks2[i+1,2])-(licks2[i,2])<5){
    licks2$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks2$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks2$licks[i+1]-licks2$licks[i]
    intervalL2 <- rbind(intervalL2,interval)
  }
}
#get rid of events only 1 second long
if (length(licks2$licks)>0){
licks2 <- group_by(licks2, EventNumber) %>% filter(n()>1)
#calculating mean
aveL2 <-aggregate(licks2[, 3], list(licks2$EventNumber), max) 
aveL2 <-cbind(aveL2, aggregate(licks2[, 3], list(licks2$EventNumber), length)) 
aveL2 <- aveL2[,-3]
aveL2 <- cbind(aveL2, mean(aveL2[,3]), mean(licks2$datum.licks.),sum(aveL2[,3]),mean(intervalL2), 
               auc((1:length(licks2$licks)),licks2$datum.licks.),genotype1, taste6)
colnames(aveL2)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype", "Tastant")
aveL2 <- aveL2[mixedorder(aveL2$`Event Number`),]}
if (length(licks2$licks)==0){
  aveL2 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL2)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype", "Tastant")
}
if (is.null(intervalF2)==TRUE){
  intervalF2 <- NA
}
if (length(feeds2)>0){
  aveW2 <- cbind(aveF2$`Mean Length of all Events`[1],aveL2$`Mean Length of all Events`[1],mean(aveF2$`Mean Intensity(per event)`),
                 mean(aveL2$`Mean Intensity(per event)`),aveF2$`Total time (All Events)`[1],aveL2$`Total time (All Events)`[1],
                 aveF2$`Mean interval between Events`[1],aveL2$`Mean interval between Events`[1], length(aveF2$`Event Number`),
                 length(aveL2$`Event Number`),
                 aveF2$`Feed AUC (a.u^2)`, aveL2$`Lick AUC (a.u^2)`,max(w2f1$datum.feeds.),length(w2f1$datum.feeds.), 
                 auc((1:length(w2f1$feeds)),w2f1$datum.feeds.),intervalF2[1],feeds2$feeds[1],
                 sum((aggregate(licks2[, 2], list(licks2$EventNumber), max) < feeds2$feeds[1]),na.rm=TRUE),
                 genotype1,taste6, runinfo, 2)
  aucvec2 <- lapply(aucfeed2, function(d) (auc(d[[2]],d[[3]])))
  aucvec2 <- unlist(aucvec2)
  aucvec2 <- as.vector(aucvec2)
  deltaintw2 <- data.frame(aveF2$`Event Number`, ((aveF2$`Mean Intensity(per event)`
                                                   -mean(aveF2$`Mean Intensity(per event)`))/mean(aveF2$`Mean Intensity(per event)`)),((aveF2$`Length of Event`-mean(aveF2$`Length of Event`))/mean(aveF2$`Length of Event`))
                           ,(aucvec2-mean(aucvec2))/mean(aucvec2),
                           ((intervalF2[1]-mean(intervalF2))/mean(intervalF2)),((intervalF2[2]-mean(intervalF2))/mean(intervalF2)),((intervalF2[3]-mean(intervalF2))/mean(intervalF2)),
                           ((intervalF2[4]-mean(intervalF2))/mean(intervalF2)),((intervalF2[5]-mean(intervalF2))/mean(intervalF2)),
                           genotype1, taste6,paste(runinfo,"Well 2", sep=" "))
  }
if (length(feeds2)==0){
  aveW2 <- cbind(aveF2$`Mean Length of all Events`[1],aveL2$`Mean Length of all Events`[1],aveF2$`Mean Intensity (all Events)`[1],
                 aveL2$`Mean Intensity (all Events)`[1],aveF2$`Total time (All Events)`[1],aveL2$`Total time (All Events)`[1],
                 aveF2$`Mean interval between Events`[1],aveL2$`Mean interval between Events`[1], length(aveF2$`Event Number`),
                 length(aveL2$`Event Number`),
                 aveF2$`Feed AUC (a.u^2)`, aveL2$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype1,taste6, runinfo, 2)
  deltaintw2 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype1,taste6,paste(runinfo,"Well 2", sep=" "))}
colnames(deltaintw2) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW2) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")

###WELL 3 CHUNK
#grouping feeding events- Well3
j <- 1
intervalF3 <- NULL
if(length(feeds3$feeds)){
feeds3$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds3$feeds)-1)){
  if (length(feeds3$feeds)==0){
    break
  }
  else if ((feeds3[i+1,2])-(feeds3[i,2])<5){
    feeds3$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds3$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds3$feeds[i+1]-feeds3$feeds[i]
    intervalF3 <- rbind(intervalF3,interval)
  }
}

aucfeed3 <- split(feeds3,feeds3$EventNumber)
aucfeed3 <- aucfeed3[sapply(aucfeed3, function(d) max(d[[3]]) > mode3+100)]
feeds3 <- do.call("rbind",aucfeed3)
feeds3 <- feeds3[order(feeds3$feeds),]
j <- 1
intervalF3 <- NULL
if(length(feeds3$feeds)){
  feeds3$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds3$feeds)-1)){
  if (length(feeds3$feeds)==0){
    break
  }
  else if ((feeds3[i+1,2])-(feeds3[i,2])<5){
    feeds3$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds3$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds3$feeds[i+1]-feeds3$feeds[i]
    intervalF3 <- rbind(intervalF3,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds3$feeds)>0){
feeds3 <- group_by(feeds3, EventNumber) %>% filter(n()>1)
#calculating mean
aveF3 <-aggregate(feeds3[, 3], list(feeds3$EventNumber), max)
aveF3 <-cbind(aveF3, aggregate(feeds3[, 3], list(feeds3$EventNumber), length)) 
aveF3 <- aveF3[,-3]
w3f1 <- feeds3[feeds3$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF3 <- cbind(aveF3, mean(aveF3[,3]), mean(feeds3$datum.feeds.),sum(aveF3[,3]),mean(intervalF3), 
               auc((1:length(feeds3$feeds)),feeds3$datum.feeds.),genotype2, taste5)
colnames(aveF3) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
aveF3 <- aveF3[mixedorder(aveF3$`Event Number`),]}
if (length(feeds3$feeds)==0){
  aveF3 <- data.frame(matrix(ncol=10, nrow=1))
  colnames(aveF3) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
}
#grouping Licking events- Well 3
j <- 1
intervalL3 <- NULL
if (length(licks3$licks)>0){
licks3$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks3$licks)-1)){
  if (length(licks3$licks)==0){
    break
  }
  else if ((licks3[i+1,2])-(licks3[i,2])<5){
    licks3$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks3$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks3$licks[i+1]-licks3$licks[i]
    intervalL3 <- rbind(intervalL3,interval)
  }
}
auclick3 <- split(licks3,licks3$EventNumber)
auclick3 <- auclick3[sapply(auclick3, function(d) max(d[[3]]) < mode3+100)]
licks3 <- do.call("rbind",auclick3)
licks3 <- licks3[order(licks3$licks),]
j <- 1
intervalL3 <- NULL
if (length(licks3$licks)>0){
  licks3$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks3$licks)-1)){
  if (length(licks3$licks)==0){
    break
  }
  else if ((licks3[i+1,2])-(licks3[i,2])<5){
    licks3$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks3$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks3$licks[i+1]-licks3$licks[i]
    intervalL3 <- rbind(intervalL3,interval)
  }
}
#get rid of events only 1 second long
if (length(licks3$licks)>0){
licks3 <- group_by(licks3, EventNumber) %>% filter(n()>1)
#calculating mean
aveL3 <-aggregate(licks3[, 3], list(licks3$EventNumber), max)
aveL3 <-cbind(aveL3, aggregate(licks3[, 3], list(licks3$EventNumber), length)) 
aveL3 <- aveL3[,-3]
aveL3 <- cbind(aveL3, mean(aveL3[,3]), mean(licks3$datum.licks.),sum(aveL3[,3]),mean(intervalL3), 
               auc((1:length(licks3$licks)),licks3$datum.licks.),genotype2, taste5)
colnames(aveL3)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype", "Tastant")
aveL3 <- aveL3[mixedorder(aveL3$`Event Number`),]}
if (length(licks3$licks)==0){
  aveL3 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL3)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype", "Tastant")
}
if (is.null(intervalF3)==TRUE){
  intervalF3 <- NA}
if (length(feeds3)>0){
  aveW3 <- cbind(aveF3$`Mean Length of all Events`[1],aveL3$`Mean Length of all Events`[1],mean(aveF3$`Mean Intensity(per event)`),
                 mean(aveL3$`Mean Intensity(per event)`),aveF3$`Total time (All Events)`[1],aveL3$`Total time (All Events)`[1],
                 aveF3$`Mean interval between Events`[1],aveL3$`Mean interval between Events`[1], 
                 length(aveF3$`Event Number`),length(aveL3$`Event Number`),
                 aveF3$`Feed AUC (a.u^2)`, aveL3$`Lick AUC (a.u^2)`,max(w3f1$datum.feeds.),length(w3f1$datum.feeds.), 
                 auc((1:length(w3f1$feeds)),w3f1$datum.feeds.),intervalF3[1],feeds3$feeds[1],
                 sum((aggregate(licks3[, 2], list(licks3$EventNumber), max) < feeds3$feeds[1]),na.rm=TRUE),
                 genotype2, taste5,runinfo,3)
  aucvec3 <- lapply(aucfeed3, function(d) (auc(d[[2]],d[[3]])))
  aucvec3 <- unlist(aucvec3)
  aucvec3 <- as.vector(aucvec3)
  deltaintw3 <- data.frame(aveF3$`Event Number`, ((aveF3$`Mean Intensity(per event)`
                                                   -mean(aveF3$`Mean Intensity(per event)`))/mean(aveF3$`Mean Intensity(per event)`)),
                           ((aveF3$`Length of Event`-mean(aveF3$`Length of Event`))/mean(aveF3$`Length of Event`))
                           ,(aucvec3-mean(aucvec3))/mean(aucvec3),
                           ((intervalF3[1]-mean(intervalF3))/mean(intervalF3)),((intervalF3[2]-mean(intervalF3))/mean(intervalF3)),((intervalF3[3]-mean(intervalF3))/mean(intervalF3)),
                           ((intervalF3[4]-mean(intervalF3))/mean(intervalF3)),((intervalF3[5]-mean(intervalF3))/mean(intervalF3)),
                           genotype2, taste5,paste(runinfo,"Well 3", sep=" "))}
if (length(feeds3)==0){
  aveW3 <- cbind(aveF3$`Mean Length of all Events`[1],aveL3$`Mean Length of all Events`[1],aveF3$`Mean Intensity (all Events)`[1],
                 aveL3$`Mean Intensity (all Events)`[1],aveF3$`Total time (All Events)`[1],aveL3$`Total time (All Events)`[1],
                 aveF3$`Mean interval between Events`[1],aveL3$`Mean interval between Events`[1], 
                 length(aveF3$`Event Number`),length(aveL3$`Event Number`),
                 aveF3$`Feed AUC (a.u^2)`, aveL3$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype2, taste5,runinfo,3)  
  deltaintw3 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype2,taste5,paste(runinfo,"Well 3", sep=" "))
}
colnames(deltaintw3) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW3) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
###WELL 4 CHUNK
#grouping feeding events- Well4
j <- 1
intervalF4 <- NULL
if (length(feeds4$feeds)>0){
feeds4$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds4$feeds)-1)){
  if (length(feeds4$feeds)==0){
    break
  }
  else if ((feeds4[i+1,2])-(feeds4[i,2])<5){
    feeds4$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds4$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds4$feeds[i+1]-feeds4$feeds[i]
    intervalF4 <- rbind(intervalF4,interval)
  }
}
aucfeed4 <- split(feeds4,feeds4$EventNumber)
aucfeed4 <- aucfeed4[sapply(aucfeed4, function(d) max(d[[3]]) >mode4+100)]
feeds4 <- do.call("rbind",aucfeed4)
feeds4 <- feeds4[order(feeds4$feeds),]
j <- 1
intervalF4 <- NULL
if (length(feeds4$feeds)>0){
  feeds4$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds4$feeds)-1)){
  if (length(feeds4$feeds)==0){
    break
  }
  else if ((feeds4[i+1,2])-(feeds4[i,2])<5){
    feeds4$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds4$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds4$feeds[i+1]-feeds4$feeds[i]
    intervalF4 <- rbind(intervalF4,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds4$feeds)>0){
feeds4 <- group_by(feeds4, EventNumber) %>% filter(n()>1)
#calculating mean
aveF4 <-aggregate(feeds4[, 3], list(feeds4$EventNumber), max) 
aveF4 <-cbind(aveF4, aggregate(feeds4[, 3], list(feeds4$EventNumber), length)) 
aveF4 <- aveF4[,-3]
w4f1 <- feeds4[feeds4$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF4 <- cbind(aveF4, mean(aveF4[,3]), mean(feeds4$datum.feeds.),sum(aveF4[,3]),mean(intervalF4),
               auc((1:length(feeds4$feeds)),feeds4$datum.feeds.),genotype2,taste6)
colnames(aveF4) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype", "Tastant")
aveF4 <- aveF4[mixedorder(aveF4$`Event Number`),]}
if (length(feeds4$feeds)==0){
  aveF4 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF4) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
}
#grouping Licking events- Well 4
j <- 1
intervalL4 <- NULL
if (length(licks4$licks)>0){
licks4$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks4$licks)-1)){
  if (length(licks4$licks)==0){
    break
  }
  else if ((licks4[i+1,2])-(licks4[i,2])<5){
    licks4$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks4$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks4$licks[i+1]-licks4$licks[i]
    intervalL4 <- rbind(intervalL4,interval)
  }
}
auclick4 <- split(licks4,licks4$EventNumber)
auclick4 <- auclick4[sapply(auclick4, function(d) max(d[[3]]) < mode4+100)]
licks4 <- do.call("rbind",auclick4)
licks4 <- licks4[order(licks4$licks),]
j <- 1
intervalL4 <- NULL
if (length(licks4$licks)>0){
  licks4$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks4$licks)-1)){
  if (length(licks4$licks)==0){
    break
  }
  else if ((licks4[i+1,2])-(licks4[i,2])<5){
    licks4$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks4$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks4$licks[i+1]-licks4$licks[i]
    intervalL4 <- rbind(intervalL4,interval)
  }
}

#get rid of events only 1 second long
if (length(licks4$licks)>0){
licks4 <- group_by(licks4, EventNumber) %>% filter(n()>1)
#calculating mean
aveL4 <-aggregate(licks4[, 3], list(licks4$EventNumber), max) 
aveL4 <-cbind(aveL4, aggregate(licks4[, 3], list(licks4$EventNumber), length)) 
aveL4 <- aveL4[,-3]
aveL4 <- cbind(aveL4, mean(aveL4[,3]), mean(licks4$datum.licks.),sum(aveL4[,3]),mean(intervalL4), auc((1:length(licks4$licks)),licks4$datum.licks.),genotype2,taste6)
colnames(aveL4)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype", "Tastant")
aveL4 <- aveL4[mixedorder(aveL4$`Event Number`),]}
if (length(licks4$licks)==0){
  aveL4 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL4)<- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype", "Tastant")
}
if (is.null(intervalF4)==TRUE){
  intervalF4 <- NA
}
if (length(feeds4)>0){
  aveW4 <- cbind(aveF4$`Mean Length of all Events`[1],aveL4$`Mean Length of all Events`[1],mean(aveF4$`Mean Intensity(per event)`),
                 mean(aveL4$`Mean Intensity(per event)`),aveF4$`Total time (All Events)`[1],aveL4$`Total time (All Events)`[1],
                 aveF4$`Mean interval between Events`[1],aveL4$`Mean interval between Events`[1], 
                 length(aveF4$`Event Number`),length(aveL4$`Event Number`),
                 aveF4$`Feed AUC (a.u^2)`, aveL4$`Lick AUC (a.u^2)`,max(w4f1$datum.feeds.),length(w4f1$datum.feeds.), 
                 auc((1:length(w4f1$feeds)),w4f1$datum.feeds.),intervalF4[1],feeds4$feeds[1],
                 sum((aggregate(licks4[, 2], list(licks4$EventNumber), max) < feeds4$feeds[1]),na.rm=TRUE),
                 genotype2, taste6, runinfo, 4)
  aucvec4 <- lapply(aucfeed4, function(d) (auc(d[[2]],d[[3]])))
  aucvec4 <- unlist(aucvec4)
  aucvec4 <- as.vector(aucvec4)
  deltaintw4 <- data.frame(aveF4$`Event Number`, ((aveF4$`Mean Intensity(per event)`
                                                   -mean(aveF4$`Mean Intensity(per event)`))/mean(aveF4$`Mean Intensity(per event)`)),((aveF4$`Length of Event`-mean(aveF4$`Length of Event`))/
                                                                                                                                         mean(aveF4$`Length of Event`)),(aucvec4-mean(aucvec4))/mean(aucvec4),
                           ((intervalF4[1]-mean(intervalF4))/mean(intervalF4)),((intervalF4[2]-mean(intervalF4))/mean(intervalF4)),((intervalF4[3]-mean(intervalF4))/mean(intervalF4)),
                           ((intervalF4[4]-mean(intervalF4))/mean(intervalF4)),((intervalF4[5]-mean(intervalF4))/mean(intervalF4)),
                           genotype2, taste6,paste(runinfo,"Well 4", sep=" "))}
if (length(feeds4)==0){
  aveW4 <- cbind(aveF4$`Mean Length of all Events`[1],aveL4$`Mean Length of all Events`[1],mean(aveF4$`Mean Intensity(per event)`),
                 mean(aveL4$`Mean Intensity(per event)`),aveF4$`Total time (All Events)`[1],aveL4$`Total time (All Events)`[1],
                 aveF4$`Mean interval between Events`[1],aveL4$`Mean interval between Events`[1], 
                 length(aveF4$`Event Number`),length(aveL4$`Event Number`),
                 aveF4$`Feed AUC (a.u^2)`, aveL4$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype2, taste6, runinfo, 4)
  deltaintw4 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype2,taste6,paste(runinfo,"Well 4", sep=" "))}
colnames(deltaintw4) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW4) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")

###WELL 5 CHUNK
#grouping feeding events- Well5
j <- 1
intervalF5 <- NULL
if (length(feeds5$feeds)>0){
feeds5$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds5$feeds)-1)){
  if (length(feeds5$feeds)==0){
    break
  }
  else if ((feeds5[i+1,2])-(feeds5[i,2])<5){
    feeds5$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds5$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds5$feeds[i+1]-feeds5$feeds[i]
    intervalF5 <- rbind(intervalF5,interval)
  }
}
aucfeed5 <- split(feeds5,feeds5$EventNumber)
aucfeed5 <- aucfeed5[sapply(aucfeed5, function(d) max(d[[3]]) > mode5+ 100)]
feeds5 <- do.call("rbind",aucfeed5)
feeds5 <- feeds5[order(feeds5$feeds),]
j <- 1
intervalF5 <- NULL
if (length(feeds5$feeds)>0){
  feeds5$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds5$feeds)-1)){
  if (length(feeds5$feeds)==0){
    break
  }
  else if ((feeds5[i+1,2])-(feeds5[i,2])<5){
    feeds5$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds5$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds5$feeds[i+1]-feeds5$feeds[i]
    intervalF5 <- rbind(intervalF5,interval)
  }
}
#get rid of events only 1 second long
if (length(feeds5$feeds)>0){
feeds5 <- group_by(feeds5, EventNumber) %>% filter(n()>1)
#calculating mean
aveF5 <-aggregate(feeds5[, 3], list(feeds5$EventNumber), max) 
aveF5 <-cbind(aveF5, aggregate(feeds5[, 3], list(feeds5$EventNumber), length)) 
aveF5 <- aveF5[,-3]
w5f1 <- feeds5[feeds5$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF5 <- cbind(aveF5, mean(aveF5[,3]), mean(feeds5$datum.feeds.),sum(aveF5[,3]),mean(intervalF5), 
               auc((1:length(feeds5$feeds)),feeds5$datum.feeds.),genotype3,taste5)
colnames(aveF5) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
aveF5 <- aveF5[mixedorder(aveF5$`Event Number`),]}
if (length(feeds5$feeds)==0){
  aveF5 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF5) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)", "Genotype", "Tastant")
  }
#grouping Licking events- Well 5
j <- 1
intervalL5 <- NULL
if (length(licks5$licks)>0){
licks5$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks5$licks)-1)){
  if (length(licks5$licks)==0){
    break
  }
  else if ((licks5[i+1,2])-(licks5[i,2])<5){
    licks5$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks5$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks5$licks[i+1]-licks5$licks[i]
    intervalL5 <- rbind(intervalL5,interval)
  }
}
auclick5 <- split(licks5,licks5$EventNumber)
auclick5 <- auclick5[sapply(auclick5, function(d) max(d[[3]]) < mode5+100)]
licks5 <- do.call("rbind",auclick5)
licks5 <- licks5[order(licks5$licks),]
j <- 1
intervalL5 <- NULL
if (length(licks5$licks)>0){
  licks5$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks5$licks)-1)){
  if (length(licks5$licks)==0){
    break
  }
  else if ((licks5[i+1,2])-(licks5[i,2])<5){
    licks5$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks5$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks5$licks[i+1]-licks5$licks[i]
    intervalL5 <- rbind(intervalL5,interval)
  }
}

#get rid of events only 1 second long
if (length(licks5$licks)>0){
licks5 <- group_by(licks5, EventNumber) %>% filter(n()>1)
#calculating mean
aveL5 <-aggregate(licks5[, 3], list(licks5$EventNumber), max) 
aveL5 <-cbind(aveL5, aggregate(licks5[, 3], list(licks5$EventNumber), length)) 
aveL5 <- aveL5[,-3]
aveL5 <- cbind(aveL5, mean(aveL5[,3]), mean(licks5$datum.licks.),sum(aveL5[,3]),mean(intervalL5),
               auc((1:length(licks5$licks)),licks5$datum.licks.), genotype3,taste5)
colnames(aveL5) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
aveL5 <- aveL5[mixedorder(aveL5$`Event Number`),]}
if (length(licks5$licks)==0){
  aveL5 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL5) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
}
if (is.null(intervalF5)==TRUE){
  intervalF5 <- NA
}
if (length(feeds5)>0){
  aveW5 <- cbind(aveF5$`Mean Length of all Events`[1],aveL5$`Mean Length of all Events`[1],mean(aveF5$`Mean Intensity(per event)`),
                 mean(aveL5$`Mean Intensity(per event)`),aveF5$`Total time (All Events)`[1],aveL5$`Total time (All Events)`[1],
                 aveF5$`Mean interval between Events`[1],aveL5$`Mean interval between Events`[1], 
                 length(aveF5$`Event Number`),length(aveL5$`Event Number`),
                 aveF5$`Feed AUC (a.u^2)`, aveL5$`Lick AUC (a.u^2)`,max(w5f1$datum.feeds.),length(w5f1$datum.feeds.), 
                 auc((1:length(w5f1$feeds)),w5f1$datum.feeds.),intervalF5[1],feeds5$feeds[1],
                 sum((aggregate(licks5[, 2], list(licks5$EventNumber), max) < feeds5$feeds[1]),na.rm=TRUE),
                 genotype3,taste5, runinfo, 5)
  aucvec5 <- lapply(aucfeed5, function(d) (auc(d[[2]],d[[3]])))
  aucvec5 <- unlist(aucvec5)
  aucvec5 <- as.vector(aucvec5)
  deltaintw5 <- data.frame(aveF5$`Event Number`, ((aveF5$`Mean Intensity(per event)`
                                                   -mean(aveF5$`Mean Intensity(per event)`))/mean(aveF5$`Mean Intensity(per event)`)),((aveF5$`Length of Event`-mean(aveF5$`Length of Event`))/
                                                                                                                                         mean(aveF5$`Length of Event`)),(aucvec5-mean(aucvec5))/mean(aucvec5),
                           ((intervalF5[1]-mean(intervalF5))/mean(intervalF5)),((intervalF5[2]-mean(intervalF5))/mean(intervalF5)),((intervalF5[3]-mean(intervalF5))/mean(intervalF5)),
                           ((intervalF5[4]-mean(intervalF5))/mean(intervalF5)),((intervalF5[5]-mean(intervalF5))/mean(intervalF5)),
                           genotype3, taste5,paste(runinfo,"Well 5", sep=" "))}
if (length(feeds5)==0){
  aveW5 <- cbind(aveF5$`Mean Length of all Events`[1],aveL5$`Mean Length of all Events`[1],mean(aveF5$`Mean Intensity(per event)`),
                 mean(aveL5$`Mean Intensity(per event)`),aveF5$`Total time (All Events)`[1],aveL5$`Total time (All Events)`[1],
                 aveF5$`Mean interval between Events`[1],aveL5$`Mean interval between Events`[1], 
                 length(aveF5$`Event Number`),length(aveL5$`Event Number`),
                 aveF5$`Feed AUC (a.u^2)`, aveL5$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype3,taste5, runinfo, 5)
  deltaintw5 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype3,taste5,paste(runinfo,"Well 5", sep=" "))}
colnames(deltaintw5) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC',"Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5", 'Genotype','Tastant','Run.Well.ID')
colnames(aveW5) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
###WELL 6 CHUNK
#grouping feeding events- Well6
j <- 1
intervalF6 <- NULL
if (length(feeds6$feeds)>0){
feeds6$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds6$feeds)-1)){
  if (length(feeds6$feeds)==0){
    break
  }
  else if ((feeds6[i+1,2])-(feeds6[i,2])<5){
    feeds6$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds6$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds6$feeds[i+1]-feeds6$feeds[i]
    intervalF6 <- rbind(intervalF6,interval)
  }
}
aucfeed6 <- split(feeds6,feeds6$EventNumber)
aucfeed6 <- aucfeed6[sapply(aucfeed6, function(d) max(d[[3]]) > mode6+100)]
feeds6 <- do.call("rbind",aucfeed6)
feeds6 <- feeds6[order(feeds6$feeds),]
j <- 1
intervalF6 <- NULL
if (length(feeds6$feeds)>0){
  feeds6$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds6$feeds)-1)){
  if (length(feeds6$feeds)==0){
    break
  }
  else if ((feeds6[i+1,2])-(feeds6[i,2])<5){
    feeds6$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds6$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds6$feeds[i+1]-feeds6$feeds[i]
    intervalF6 <- rbind(intervalF6,interval)
  }
}
#get rid of events only 1 second long
if (length(feeds6$feeds)>0){
feeds6 <- group_by(feeds6, EventNumber) %>% filter(n()>1)
#calculating mean
aveF6 <-aggregate(feeds6[, 3], list(feeds6$EventNumber), max)
aveF6 <-cbind(aveF6, aggregate(feeds6[, 3], list(feeds6$EventNumber), length)) 
aveF6 <- aveF6[,-3]
w6f1 <- feeds6[feeds6$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF6 <- cbind(aveF6, mean(aveF6[,3]), mean(feeds6$datum.feeds.),sum(aveF6[,3]),mean(intervalF6),
               auc((1:length(feeds6$feeds)),feeds6$datum.feeds.),genotype3,taste6)
colnames(aveF6) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
aveF6 <- aveF6[mixedorder(aveF6$`Event Number`),]}
if (length(feeds6$feeds)==0){
  aveF6 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF6) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
}
#grouping Licking events- Well 6
j <- 1
intervalL6 <- NULL
if (length(licks6$licks)>0){
licks6$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks6$licks)-1)){
  if (length(licks6$licks)==0){
    break
  }
  else if ((licks6[i+1,2])-(licks6[i,2])<5){
    licks6$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks6$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks6$licks[i+1]-licks6$licks[i]
    intervalL6 <- rbind(intervalL6,interval)
  }
}
auclick6 <- split(licks6,licks6$EventNumber)
auclick6 <- auclick6[sapply(auclick6, function(d) max(d[[3]]) < mode6+100)]
licks6 <- do.call("rbind",auclick6)
licks6 <- licks6[order(licks6$licks),]
j <- 1
intervalL6 <- NULL
if (length(licks6$licks)>0){
  licks6$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks6$licks)-1)){
  if (length(licks6$licks)==0){
    break
  }
  else if ((licks6[i+1,2])-(licks6[i,2])<5){
    licks6$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks6$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks6$licks[i+1]-licks6$licks[i]
    intervalL6 <- rbind(intervalL6,interval)
  }
}

#get rid of events only 1 second long
if (length(licks6$licks)>0){
licks6 <- group_by(licks6, EventNumber) %>% filter(n()>1)
#calculating mean
aveL6 <-aggregate(licks6[, 3], list(licks6$EventNumber), max) 
aveL6 <-cbind(aveL6, aggregate(licks6[, 3], list(licks6$EventNumber), length)) 
aveL6 <- aveL6[,-3]
aveL6 <- cbind(aveL6, mean(aveL6[,3]), mean(licks6$datum.licks.),sum(aveL6[,3]),mean(intervalL6), 
               auc((1:length(licks6$licks)),licks6$datum.licks.),genotype3,taste6)
colnames(aveL6) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
aveL6 <- aveL6[mixedorder(aveL6$`Event Number`),]}
if (length(licks6$licks)==0){
  aveL6 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL6) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
}
if (is.null(intervalF6)==TRUE){
  intervalF6 <- NA
}
if (length(feeds6)>0){
  aveW6 <- cbind(aveF6$`Mean Length of all Events`[1],aveL6$`Mean Length of all Events`[1],mean(aveF6$`Mean Intensity(per event)`),
                 mean(aveL6$`Mean Intensity(per event)`),aveF6$`Total time (All Events)`[1],aveL6$`Total time (All Events)`[1],
                 aveF6$`Mean interval between Events`[1],aveL6$`Mean interval between Events`[1], length(aveF6$`Event Number`),
                 length(aveL6$`Event Number`),
                 aveF6$`Feed AUC (a.u^2)`, aveL6$`Lick AUC (a.u^2)`,max(w6f1$datum.feeds.),length(w6f1$datum.feeds.), 
                 auc((1:length(w6f1$feeds)),w6f1$datum.feeds.),intervalF6[1],feeds6$feeds[1],
                 sum((aggregate(licks6[, 2], list(licks6$EventNumber), max) < feeds6$feeds[1]),na.rm=TRUE),genotype3, taste6, runinfo, 6)
  aucvec6 <- lapply(aucfeed6, function(d) (auc(d[[2]],d[[3]])))
  aucvec6 <- unlist(aucvec6)
  aucvec6 <- as.vector(aucvec6)
  deltaintw6 <- data.frame(aveF6$`Event Number`, ((aveF6$`Mean Intensity(per event)`
                                                   -mean(aveF6$`Mean Intensity(per event)`))/mean(aveF6$`Mean Intensity(per event)`)),((aveF6$`Length of Event`-mean(aveF6$`Length of Event`))/
                                                                                                                                         mean(aveF6$`Length of Event`)),(aucvec6-mean(aucvec6))/mean(aucvec6),
                           ((intervalF6[1]-mean(intervalF6))/mean(intervalF6)),((intervalF6[2]-mean(intervalF6))/mean(intervalF6)),((intervalF6[3]-mean(intervalF6))/mean(intervalF6)),
                           ((intervalF6[4]-mean(intervalF6))/mean(intervalF6)),((intervalF6[5]-mean(intervalF6))/mean(intervalF6)),
                           genotype3, taste6,paste(runinfo,"Well 6", sep=" "))}
if (length(feeds6)==0){
  aveW6 <- cbind(aveF6$`Mean Length of all Events`[1],aveL6$`Mean Length of all Events`[1],mean(aveF6$`Mean Intensity(per event)`),
                 mean(aveL6$`Mean Intensity(per event)`),aveF6$`Total time (All Events)`[1],aveL6$`Total time (All Events)`[1],
                 aveF6$`Mean interval between Events`[1],aveL6$`Mean interval between Events`[1], length(aveF6$`Event Number`),
                 length(aveL6$`Event Number`),
                 aveF6$`Feed AUC (a.u^2)`, aveL6$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype3, taste6, runinfo, 6)  
  deltaintw6 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype3,taste6,paste(runinfo,"Well 6", sep=" "))}
colnames(deltaintw6) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW6) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
###WELL 7 CHUNK
#grouping feeding events- Well7
j <- 1
intervalF7 <- NULL
if (length(feeds7$feeds)>0){
  feeds7$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds7$feeds)-1)){
  if (length(feeds7$feeds)==0){
    break
  }
  else if (((feeds7[i+1,2])-(feeds7[i,2])<5)&(length(feeds7$feeds)>0)){
    feeds7$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds7$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds7$feeds[i+1]-feeds7$feeds[i]
    intervalF7 <- rbind(intervalF7,interval)
  }
}
aucfeed7 <- split(feeds7,feeds7$EventNumber)
aucfeed7 <- aucfeed7[sapply(aucfeed7, function(d) max(d[[3]]) >mode7+ 100)]
feeds7 <- do.call("rbind",aucfeed7)
feeds7 <- feeds7[order(feeds7$feeds),]
j <- 1
intervalF7 <- NULL
if (length(feeds7$feeds)>0){
  feeds7$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds7$feeds)-1)){
  if (length(feeds7$feeds)==0){
    break
  }
  else if (((feeds7[i+1,2])-(feeds7[i,2])<5)&(length(feeds7$feeds)>0)){
    feeds7$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds7$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds7$feeds[i+1]-feeds7$feeds[i]
    intervalF7 <- rbind(intervalF7,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds7$feeds)>0){
feeds7 <- group_by(feeds7, EventNumber) %>% filter(n()>1)
#calculating mean
aveF7 <-aggregate(feeds7[, 3], list(feeds7$EventNumber), max) 
aveF7 <-cbind(aveF7, aggregate(feeds7[, 3], list(feeds7$EventNumber), length)) 
aveF7 <- aveF7[,-3]
w7f1 <- feeds7[feeds7$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF7 <- cbind(aveF7, mean(aveF7[,3]), mean(feeds7$datum.feeds.),sum(aveF7[,3]),mean(intervalF7), 
               auc((1:length(feeds7$feeds)),feeds7$datum.feeds.),genotype4,taste5)
colnames(aveF7) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
aveF7 <- aveF7[mixedorder(aveF7$`Event Number`),]}
if (length(feeds7$feeds)==0){
  aveF7 <- data.frame(matrix(ncol=10, nrow=1))
  colnames(aveF7) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
}

#grouping Licking events- Well 7
j <- 1
intervalL7 <- NULL
if (length(licks7$licks)>0){
licks7$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks7$licks)-1)){
  if (length(licks7$licks)==0){
    break
  }
  else if ((licks7[i+1,2])-(licks7[i,2])<5){
    licks7$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks7$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks7$licks[i+1]-licks7$licks[i]
    intervalL7 <- rbind(intervalL7,interval)
  }
}
auclick7 <- split(licks7,licks7$EventNumber)
auclick7 <- auclick7[sapply(auclick7, function(d) max(d[[3]]) < mode7+ 100)]
licks7 <- do.call("rbind",auclick7)
licks7 <- licks7[order(licks7$licks),]
j <- 1
intervalL7 <- NULL
if (length(licks7$licks)>0){
  licks7$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks7$licks)-1)){
  if (length(licks7$licks)==0){
    break
  }
  else if ((licks7[i+1,2])-(licks7[i,2])<5){
    licks7$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks7$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks7$licks[i+1]-licks7$licks[i]
    intervalL7 <- rbind(intervalL7,interval)
  }
}
#get rid of events only 1 second long
if (length(licks7$licks)>0){
licks7<- group_by(licks7, EventNumber) %>% filter(n()>1)
#calculating mean
aveL7 <-aggregate(licks7[, 3], list(licks7$EventNumber), max) 
aveL7 <-cbind(aveL7, aggregate(licks7[, 3], list(licks7$EventNumber), length)) 
aveL7 <- aveL7[,-3]
aveL7 <- cbind(aveL7, mean(aveL7[,3]), mean(licks7$datum.licks.),sum(aveL7[,3]),mean(intervalL7), 
               auc((1:length(licks7$licks)),licks7$datum.licks.),genotype4,taste5)
colnames(aveL7) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
aveL7 <- aveL7[mixedorder(aveL7$`Event Number`),]}
if (length(licks7$licks)==0){
  aveL7 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL7) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
}
if (is.null(intervalF7)==TRUE){
  intervalF7 <- NA
}
if (length(feeds7)>0){
  aveW7 <- cbind(aveF7$`Mean Length of all Events`[1],aveL7$`Mean Length of all Events`[1],mean(aveF7$`Mean Intensity(per event)`),
                 mean(aveL7$`Mean Intensity(per event)`),aveF7$`Total time (All Events)`[1],aveL7$`Total time (All Events)`[1],
                 aveF7$`Mean interval between Events`[1],aveL7$`Mean interval between Events`[1], 
                 length(aveF7$`Event Number`),length(aveL7$`Event Number`),
                 aveF7$`Feed AUC (a.u^2)`, aveL7$`Lick AUC (a.u^2)`, max(w7f1$datum.feeds.),length(w7f1$datum.feeds.), 
                 auc((1:length(w7f1$feeds)),w7f1$datum.feeds.),intervalF7[1],feeds7$feeds[1],
                 sum((aggregate(licks7[, 2], list(licks7$EventNumber), max) < feeds7$feeds[1]),na.rm=TRUE),
                 genotype4, taste5, runinfo, 7)
  aucvec7 <- lapply(aucfeed7, function(d) (auc(d[[2]],d[[3]])))
  aucvec7 <- unlist(aucvec7)
  aucvec7 <- as.vector(aucvec7)
  deltaintw7 <- data.frame(aveF7$`Event Number`, ((aveF7$`Mean Intensity(per event)`
                                                   -mean(aveF7$`Mean Intensity(per event)`))/mean(aveF7$`Mean Intensity(per event)`)),((aveF7$`Length of Event`-mean(aveF7$`Length of Event`))/
                                                                                                                                         mean(aveF7$`Length of Event`)),(aucvec7-mean(aucvec7))/mean(aucvec7),
                           ((intervalF7[1]-mean(intervalF7))/mean(intervalF7)),((intervalF7[2]-mean(intervalF7))/mean(intervalF7)),((intervalF7[3]-mean(intervalF7))/mean(intervalF7)),
                           ((intervalF7[4]-mean(intervalF7))/mean(intervalF7)),((intervalF7[5]-mean(intervalF7))/mean(intervalF7)),
                           genotype4, taste5,paste(runinfo,"Well 7", sep=" "))}
if (length(feeds7)==0){
  aveW7 <- cbind(aveF7$`Mean Length of all Events`[1],aveL7$`Mean Length of all Events`[1],mean(aveF7$`Mean Intensity(per event)`),
                 mean(aveL7$`Mean Intensity(per event)`),aveF7$`Total time (All Events)`[1],aveL7$`Total time (All Events)`[1],
                 aveF7$`Mean interval between Events`[1],aveL7$`Mean interval between Events`[1], 
                 length(aveF7$`Event Number`),length(aveL7$`Event Number`),
                 aveF7$`Feed AUC (a.u^2)`, aveL7$`Lick AUC (a.u^2)`, 0,0, 
                 0,0,NA,NA,genotype4, taste5, runinfo, 7)
  deltaintw7 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype4,taste5,paste(runinfo,"Well 7", sep=" "))}
colnames(deltaintw7) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC',"Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5", 'Genotype','Tastant','Run.Well.ID')
colnames(aveW7) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
###WELL 8 CHUNK
#grouping feeding events- Well8
j <- 1
intervalF8 <- NULL
if (length(feeds8$feeds)>0){
feeds8$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds8$feeds)-1)){
  if (length(feeds8$feeds)==0){
    break
  }
  else if ((feeds8[i+1,2])-(feeds8[i,2])<5){
    feeds8$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds8$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds8$feeds[i+1]-feeds8$feeds[i]
    intervalF8 <- rbind(intervalF8,interval)
  }
}
aucfeed8 <- split(feeds8,feeds8$EventNumber)
aucfeed8 <- aucfeed8[sapply(aucfeed8, function(d) max(d[[3]]) > mode8+ 100)]
feeds8 <- do.call("rbind",aucfeed8)
feeds8 <- feeds8[order(feeds8$feeds),]
j <- 1
intervalF8 <- NULL
if (length(feeds8$feeds)>0){
  feeds8$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds8$feeds)-1)){
  if (length(feeds8$feeds)==0){
    break
  }
  else if ((feeds8[i+1,2])-(feeds8[i,2])<5){
    feeds8$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds8$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds8$feeds[i+1]-feeds8$feeds[i]
    intervalF8 <- rbind(intervalF8,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds8$feeds)>0){
feeds8 <- group_by(feeds8, EventNumber) %>% filter(n()>1)
#calculating mean
aveF8 <-aggregate(feeds8[, 3], list(feeds8$EventNumber), max)
aveF8 <-cbind(aveF8, aggregate(feeds8[, 3], list(feeds8$EventNumber), length)) 
aveF8 <- aveF8[,-3]
w8f1 <- feeds8[feeds8$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF8 <- cbind(aveF8, mean(aveF8[,3]), mean(feeds8$datum.feeds.),sum(aveF8[,3]),mean(intervalF8), 
               auc((1:length(feeds8$feeds)),feeds8$datum.feeds.),genotype4,taste6)
colnames(aveF8) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
aveF8 <- aveF8[mixedorder(aveF8$`Event Number`),]}
if (length(feeds8$feeds)==0){
  aveF8 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF8) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
}

#grouping Licking events- Well 8
j <- 1
intervalL8 <- NULL
if (length(licks8$licks)>0){
licks8$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks8$licks)-1)){
  if (length(licks8$licks)==0){
    break
  }
  else if ((licks8[i+1,2])-(licks8[i,2])<5){
    licks8$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks8$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks8$licks[i+1]-licks8$licks[i]
    intervalL8 <- rbind(intervalL8,interval)
  }
}
auclick8 <- split(licks8,licks8$EventNumber)
auclick8 <- auclick8[sapply(auclick8, function(d) max(d[[3]]) < mode8+100)]  
licks8 <- do.call("rbind",auclick8)
licks8 <- licks8[order(licks8$licks),]
j <- 1
intervalL8 <- NULL
if (length(licks8$licks)>0){
  licks8$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks8$licks)-1)){
  if (length(licks8$licks)==0){
    break
  }
  else if ((licks8[i+1,2])-(licks8[i,2])<5){
    licks8$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks8$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks8$licks[i+1]-licks8$licks[i]
    intervalL8 <- rbind(intervalL8,interval)
  }
}

#get rid of events only 1 second long
if (length(licks8$licks)>0){
licks8 <- group_by(licks8, EventNumber) %>% filter(n()>1)
#calculating mean
aveL8 <-aggregate(licks8[, 3], list(licks8$EventNumber), max) 
aveL8 <-cbind(aveL8, aggregate(licks8[, 3], list(licks8$EventNumber), length)) 
aveL8 <- aveL8[,-3]
aveL8 <- cbind(aveL8, mean(aveL8[,3]), mean(licks8$datum.licks.),sum(aveL8[,3]),mean(intervalL8), 
               auc((1:length(licks8$licks)),licks8$datum.licks.),genotype4,taste6)
colnames(aveL8) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
aveL8 <- aveL8[mixedorder(aveL8$`Event Number`),]}
if (length(licks8$licks)==0){
  aveL8 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL8) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events", "Lick AUC (a.u^2)","Genotype", "Tastant")
}
if (is.null(intervalF8)==TRUE){
  intervalF8 <- NA
}
if (length(feeds8)>0){
  aveW8 <- cbind(aveF8$`Mean Length of all Events`[1],aveL8$`Mean Length of all Events`[1],mean(aveF8$`Mean Intensity(per event)`),
                 mean(aveL8$`Mean Intensity(per event)`),aveF8$`Total time (All Events)`[1],aveL8$`Total time (All Events)`[1],
                 aveF8$`Mean interval between Events`[1],aveL8$`Mean interval between Events`[1], 
                 length(aveF8$`Event Number`),length(aveL8$`Event Number`),
                 aveF8$`Feed AUC (a.u^2)`, aveL8$`Lick AUC (a.u^2)`,max(w8f1$datum.feeds.),length(w8f1$datum.feeds.), 
                 auc((1:length(w8f1$feeds)),w8f1$datum.feeds.),intervalF8[1],feeds8$feeds[1],
                 sum((aggregate(licks8[, 2], list(licks8$EventNumber), max) < feeds8$feeds[1]),na.rm=TRUE),
                 genotype4, taste6,runinfo,8)
  aucvec8 <- lapply(aucfeed8, function(d) (auc(d[[2]],d[[3]])))
  aucvec8 <- unlist(aucvec8)
  aucvec8 <- as.vector(aucvec8)
  deltaintw8 <- data.frame(aveF8$`Event Number`, ((aveF8$`Mean Intensity(per event)`
                                                   -mean(aveF8$`Mean Intensity(per event)`))/mean(aveF8$`Mean Intensity(per event)`)),((aveF8$`Length of Event`-mean(aveF8$`Length of Event`))/
                                                                                                                                         mean(aveF8$`Length of Event`)),(aucvec8-mean(aucvec8))/mean(aucvec8),
                           ((intervalF8[1]-mean(intervalF8))/mean(intervalF8)),((intervalF8[2]-mean(intervalF8))/mean(intervalF8)),((intervalF8[3]-mean(intervalF8))/mean(intervalF8)),
                           ((intervalF8[4]-mean(intervalF8))/mean(intervalF8)),((intervalF8[5]-mean(intervalF8))/mean(intervalF8)),
                           genotype4, taste6,paste(runinfo,"Well 8", sep=" "))}
if (length(feeds8)==0){
  aveW8 <- cbind(aveF8$`Mean Length of all Events`[1],aveL8$`Mean Length of all Events`[1],mean(aveF8$`Mean Intensity(per event)`),
                 mean(aveL8$`Mean Intensity(per event)`),aveF8$`Total time (All Events)`[1],aveL8$`Total time (All Events)`[1],
                 aveF8$`Mean interval between Events`[1],aveL8$`Mean interval between Events`[1], 
                 length(aveF8$`Event Number`),length(aveL8$`Event Number`),
                 aveF8$`Feed AUC (a.u^2)`, aveL8$`Lick AUC (a.u^2)`,0,0, 
                 0,0,NA,NA,genotype4, taste6,runinfo,8)
  deltaintw8 <- data.frame(NA,NA,NA,NA,genotype4,taste6,paste(runinfo,"Well 8", sep=" "))}
colnames(deltaintw8) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW8) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
###WELL 9 CHUNK
#grouping feeding events- Well9
j <- 1
intervalF9 <- NULL
if (length(feeds9$feeds)){
feeds9$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds9$feeds)-1)){
  if (length(feeds9$feeds)==0){
    break
  }
  else if ((feeds9[i+1,2])-(feeds9[i,2])<5){
    feeds9$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds9$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds9$feeds[i+1]-feeds9$feeds[i]
    intervalF9 <- rbind(intervalF9,interval)
  }
}
aucfeed9 <- split(feeds9,feeds9$EventNumber)
aucfeed9 <- aucfeed9[sapply(aucfeed9, function(d) max(d[[3]]) > mode9+100)]
feeds9 <- do.call("rbind",aucfeed9)
feeds9 <- feeds9[order(feeds9$feeds),]
j <- 1
intervalF9 <- NULL
if (length(feeds9$feeds)){
  feeds9$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds9$feeds)-1)){
  if (length(feeds9$feeds)==0){
    break
  }
  else if ((feeds9[i+1,2])-(feeds9[i,2])<5){
    feeds9$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds9$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds9$feeds[i+1]-feeds9$feeds[i]
    intervalF9 <- rbind(intervalF9,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds9$feeds)>0){
feeds9 <- group_by(feeds9, EventNumber) %>% filter(n()>1)
#calculating mean
aveF9 <-aggregate(feeds9[, 3], list(feeds9$EventNumber), max) 
aveF9 <-cbind(aveF9, aggregate(feeds9[, 3], list(feeds9$EventNumber), length)) 
aveF9 <- aveF9[,-3]
w9f1 <- feeds9[feeds9$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF9 <- cbind(aveF9, mean(aveF9[,3]), mean(feeds9$datum.feeds.),sum(aveF9[,3]),mean(intervalF9), 
               auc((1:length(feeds9$feeds)),feeds9$datum.feeds.),genotype5,taste5)
colnames(aveF9) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Feed AUC (a.u^2)","Genotype","Tastant")
aveF9 <- aveF9[mixedorder(aveF9$`Event Number`),]}
if (length(feeds9$feeds)==0){
  aveF9 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF9) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events","Feed AUC (a.u^2)","Genotype", "Tastant")

  }

#grouping Licking events- Well 9
j <- 1
intervalL9 <- NULL
if (length(licks9$licks)>0){
licks9$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks9$licks)-1)){
  if (length(licks9$licks)==0){
    break
  }
  else if ((licks9[i+1,2])-(licks9[i,2])<5){
    licks9$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks9$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks9$licks[i+1]-licks9$licks[i]
    intervalL9 <- rbind(intervalL9,interval)
  }
}
auclick9 <- split(licks9,licks9$EventNumber)
auclick9 <- auclick9[sapply(auclick9, function(d) max(d[[3]]) < mode9+100)]
licks9 <- do.call("rbind",auclick9)
licks9 <- licks9[order(licks9$licks),]
j <- 1
intervalL9 <- NULL
if (length(licks9$licks)>0){
  licks9$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks9$licks)-1)){
  if (length(licks9$licks)==0){
    break
  }
  else if ((licks9[i+1,2])-(licks9[i,2])<5){
    licks9$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks9$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks9$licks[i+1]-licks9$licks[i]
    intervalL9 <- rbind(intervalL9,interval)
  }
}

#get rid of events only 1 second long
if (length(licks9$licks)>0){
licks9 <- group_by(licks9, EventNumber) %>% filter(n()>1)
#calculating mean
aveL9 <-aggregate(licks9[, 3], list(licks9$EventNumber), max)
aveL9 <-cbind(aveL9, aggregate(licks9[, 3], list(licks9$EventNumber), length)) 
aveL9 <- aveL9[,-3]
aveL9 <- cbind(aveL9, mean(aveL9[,3]), mean(licks9$datum.licks.),sum(aveL9[,3]), mean(intervalL9),
               auc((1:length(licks9$licks)),licks9$datum.licks.),genotype5,taste5)
colnames(aveL9) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
aveL9 <- aveL9[mixedorder(aveL9$`Event Number`),]}
if (length(licks9$licks)==0){
  aveL9 <- data.frame(matrix(ncol=10, nrow=1))
  colnames(aveL9) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                       "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                       "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
}
if (is.null(intervalF9)==TRUE){
  intervalF9 <- NA
}
if (length(feeds9)>0){
  aveW9 <- cbind(aveF9$`Mean Length of all Events`[1],aveL9$`Mean Length of all Events`[1],mean(aveF9$`Mean Intensity(per event)`),
                 mean(aveL9$`Mean Intensity(per event)`),aveF9$`Total time (All Events)`[1],aveL9$`Total time (All Events)`[1],
                 aveF9$`Mean interval between Events`[1],aveL9$`Mean interval between Events`[1], 
                 length(aveF9$`Event Number`),length(aveL9$`Event Number`),
                 aveF9$`Feed AUC (a.u^2)`, aveL9$`Lick AUC (a.u^2)`,max(w9f1$datum.feeds.),length(w9f1$datum.feeds.), 
                 auc((1:length(w9f1$feeds)),w9f1$datum.feeds.),intervalF9[1], feeds9$feeds[1],
                 sum((aggregate(licks9[, 2], list(licks9$EventNumber), max) < feeds9$feeds[1]),na.rm=TRUE),
                 genotype5, taste5, runinfo,9)
  aucvec9 <- lapply(aucfeed9, function(d) (auc(d[[2]],d[[3]])))
  aucvec9 <- unlist(aucvec9)
  aucvec9 <- as.vector(aucvec9)
  deltaintw9 <- data.frame(aveF9$`Event Number`, ((aveF9$`Mean Intensity(per event)`
                                                   -mean(aveF9$`Mean Intensity(per event)`))/mean(aveF9$`Mean Intensity(per event)`)),((aveF9$`Length of Event`-mean(aveF9$`Length of Event`))/
                                                                                                                                         mean(aveF9$`Length of Event`)),(aucvec9-mean(aucvec9))/mean(aucvec9),
                           ((intervalF9[1]-mean(intervalF9))/mean(intervalF9)),((intervalF9[2]-mean(intervalF9))/mean(intervalF9)),((intervalF9[3]-mean(intervalF9))/mean(intervalF9)),
                           ((intervalF9[4]-mean(intervalF9))/mean(intervalF9)),((intervalF9[5]-mean(intervalF9))/mean(intervalF9)),
                           genotype5, taste5,paste(runinfo,"Well 9", sep=" "))}

if (length(feeds9)==0){
  aveW9 <- cbind(aveF9$`Mean Length of all Events`[1],aveL9$`Mean Length of all Events`[1],mean(aveF9$`Mean Intensity(per event)`),
                 mean(aveL9$`Mean Intensity(per event)`),aveF9$`Total time (All Events)`[1],aveL9$`Total time (All Events)`[1],
                 aveF9$`Mean interval between Events`[1],aveL9$`Mean interval between Events`[1], 
                 length(aveF9$`Event Number`),length(aveL9$`Event Number`),
                 aveF9$`Feed AUC (a.u^2)`, aveL9$`Lick AUC (a.u^2)`,0,0,0,0,NA,NA,genotype5, taste5, runinfo,9)
  deltaintw9 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype5,taste5,paste(runinfo,"Well 9", sep=" "))}
colnames(deltaintw9) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                          "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW9) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                     "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                     "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                     "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                     "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")

###WELL 10 CHUNK
#grouping feeding events- Well10
j <- 1
intervalF10 <- NULL
if (length(feeds10$feeds)>0){
feeds10$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds10$feeds)-1)){
  if (length(feeds10$feeds)==0){
    break
  }
  else if ((feeds10[i+1,2])-(feeds10[i,2])<5){
    feeds10$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds10$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds10$feeds[i+1]-feeds10$feeds[i]
    intervalF10 <- rbind(intervalF10,interval)
  }
}
aucfeed10 <- split(feeds10,feeds10$EventNumber)
aucfeed10 <- aucfeed10[sapply(aucfeed10, function(d) max(d[[3]]) >mode10+100)]
feeds10 <- do.call("rbind",aucfeed10)
feeds10 <- feeds10[order(feeds10$feeds),]
j <- 1
if (length(feeds10$feeds)>0){
  feeds10$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds10$feeds)-1)){
  if (length(feeds10$feeds)==0){
    break
  }
  else if ((feeds10[i+1,2])-(feeds10[i,2])<5){
    feeds10$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds10$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds10$feeds[i+1]-feeds10$feeds[i]
    intervalF10 <- rbind(intervalF10,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds10$feeds)>0){
feeds10 <- group_by(feeds10, EventNumber) %>% filter(n()>1)
#calculating mean
aveF10 <-aggregate(feeds10[, 3], list(feeds10$EventNumber), max)
aveF10 <-cbind(aveF10, aggregate(feeds10[, 3], list(feeds10$EventNumber), length)) 
aveF10 <- aveF10[,-3]
w10f1 <- feeds10[feeds10$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF10 <- cbind(aveF10, mean(aveF10[,3]), mean(feeds10$datum.feeds.),sum(aveF10[,3]), mean(intervalF10),
                auc((1:length(feeds10$feeds)),feeds10$datum.feeds.),genotype5,taste6)
colnames(aveF10) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
aveF10 <- aveF10[mixedorder(aveF10$`Event Number`),]}
if (length(feeds10$feeds)==0){
  aveF10 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF10) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
}

#grouping Licking events- Well 10
j <- 1
intervalL10 <- NULL
if (length(licks10$licks)>0){
licks10$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks10$licks)-1)){
  if (length(licks10$licks)==0){
    break
  }
  else if ((licks10[i+1,2])-(licks10[i,2])<5){
    licks10$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks10$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks10$licks[i+1]-licks10$licks[i]
    intervalL10 <- rbind(intervalL10,interval)
  }
}
auclick10 <- split(licks10,licks10$EventNumber)
auclick10 <- auclick10[sapply(auclick10, function(d) max(d[[3]]) < mode10+100)]
licks10 <- do.call("rbind",auclick10)
licks10 <- licks10[order(licks10$licks),]
j <- 1
intervalL10 <- NULL
if (length(licks10$licks)>0){
  licks10$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks10$licks)-1)){
  if (length(licks10$licks)==0){
    break
  }
  else if ((licks10[i+1,2])-(licks10[i,2])<5){
    licks10$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks10$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks10$licks[i+1]-licks10$licks[i]
    intervalL10 <- rbind(intervalL10,interval)
  }
}

#get rid of events only 1 second long
if (length(licks10$licks)>0){
licks10 <- group_by(licks10, EventNumber) %>% filter(n()>1)
#calculating mean
aveL10 <-aggregate(licks10[, 3], list(licks10$EventNumber), max)
aveL10 <-cbind(aveL10, aggregate(licks10[, 3], list(licks10$EventNumber), length)) 
aveL10 <- aveL10[,-3]
aveL10 <- cbind(aveL10, mean(aveL10[,3]), mean(licks10$datum.licks.),sum(aveL10[,3]), mean(intervalL10), 
                auc((1:length(licks10$licks)),licks10$datum.licks.),genotype5,taste6)
colnames(aveL10) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
aveL10 <- aveL10[mixedorder(aveL10$`Event Number`),]}
if (length(licks10$licks)==0){
  aveL10 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL10) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
}
if (is.null(intervalF10)==TRUE){
  intervalF10 <- NA
}
if (length(feeds10)>0){
  aveW10 <- cbind(aveF10$`Mean Length of all Events`[1],aveL10$`Mean Length of all Events`[1],mean(aveF10$`Mean Intensity(per event)`),
                  mean(aveL10$`Mean Intensity(per event)`),aveF10$`Total time (All Events)`[1],aveL10$`Total time (All Events)`[1],
                  aveF10$`Mean interval between Events`[1],aveL10$`Mean interval between Events`[1], 
                  length(aveF10$`Event Number`),length(aveL10$`Event Number`),
                  aveF10$`Feed AUC (a.u^2)`, aveL10$`Lick AUC (a.u^2)`,max(w10f1$datum.feeds.),length(w10f1$datum.feeds.), 
                  auc((1:length(w10f1$feeds)),w10f1$datum.feeds.),intervalF10[1],feeds10$feeds[1],
                  sum((aggregate(licks10[, 2], list(licks10$EventNumber), max) < feeds10$feeds[1]),na.rm=TRUE),
                  genotype5, taste6,runinfo,10)
  aucvec10 <- lapply(aucfeed10, function(d) (auc(d[[2]],d[[3]])))
  aucvec10 <- unlist(aucvec10)
  aucvec10 <- as.vector(aucvec10)
  deltaintw10 <- data.frame(aveF10$`Event Number`, ((aveF10$`Mean Intensity(per event)`
                                                   -mean(aveF10$`Mean Intensity(per event)`))/mean(aveF10$`Mean Intensity(per event)`)),((aveF10$`Length of Event`-mean(aveF10$`Length of Event`))/
                                                                                                                                         mean(aveF10$`Length of Event`)),(aucvec10-mean(aucvec10))/mean(aucvec10),
                            ((intervalF10[1]-mean(intervalF10))/mean(intervalF10)),((intervalF10[2]-mean(intervalF10))/mean(intervalF10)),((intervalF10[3]-mean(intervalF10))/mean(intervalF10)),
                            ((intervalF10[4]-mean(intervalF10))/mean(intervalF10)),((intervalF10[5]-mean(intervalF10))/mean(intervalF10)),
                            genotype5, taste6,paste(runinfo,"Well 10", sep=" "))}

if (length(feeds10)==0){
  aveW10 <- cbind(aveF10$`Mean Length of all Events`[1],aveL10$`Mean Length of all Events`[1],aveF10$`Mean Intensity (all Events)`[1],
                  aveL10$`Mean Intensity (all Events)`[1],aveF10$`Total time (All Events)`[1],aveL10$`Total time (All Events)`[1],
                  aveF10$`Mean interval between Events`[1],aveL10$`Mean interval between Events`[1], 
                  length(aveF10$`Event Number`),length(aveL10$`Event Number`),
                  aveF10$`Feed AUC (a.u^2)`, aveL10$`Lick AUC (a.u^2)`,0,0,0,0,NA,NA,genotype5, taste6,runinfo,10)
  deltaintw10 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype5,taste6,paste(runinfo,"Well 10", sep=" "))}  

colnames(deltaintw10) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                           "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW10) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                      "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                      "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                      "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                      "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")

###WELL 11 CHUNK
#grouping feeding events- Well11
j <- 1
intervalF11 <- NULL
if (length(feeds11$feeds)>0){
feeds11$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds11$feeds)-1)){
  if (length(feeds11$feeds)==0){
    break
  }
  else if ((feeds11[i+1,2])-(feeds11[i,2])<5){
    feeds11$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds11$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds11$feeds[i+1]-feeds11$feeds[i]
    intervalF11 <- rbind(intervalF11,interval)
  }
}
aucfeed11 <- split(feeds11,feeds11$EventNumber)
aucfeed11 <- aucfeed11[sapply(aucfeed11, function(d) max(d[[3]]) >mode11+ 100)]
feeds11 <- do.call("rbind",aucfeed11)
feeds11 <- feeds11[order(feeds11$feeds),]
j <- 1
intervalF11 <- NULL
if (length(feeds11$feeds)>0){
  feeds11$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds11$feeds)-1)){
  if (length(feeds11$feeds)==0){
    break
  }
  else if ((feeds11[i+1,2])-(feeds11[i,2])<5){
    feeds11$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds11$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds11$feeds[i+1]-feeds11$feeds[i]
    intervalF11 <- rbind(intervalF11,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds11$feeds)>0){
feeds11 <- group_by(feeds11, EventNumber) %>% filter(n()>1)
#calculating mean
aveF11 <-aggregate(feeds11[, 3], list(feeds11$EventNumber), max) 
aveF11 <-cbind(aveF11, aggregate(feeds11[, 3], list(feeds11$EventNumber), length)) 
aveF11 <- aveF11[,-3]
w11f1 <- feeds11[feeds11$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF11 <- cbind(aveF11, mean(aveF11[,3]), mean(feeds11$datum.feeds.),sum(aveF11[,3]), mean(intervalF11), 
                auc((1:length(feeds11$feeds)),feeds11$datum.feeds.),genotype6,taste5)
colnames(aveF11) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
aveF11 <- aveF11[mixedorder(aveF11$`Event Number`),]}
if (length(feeds11$feeds)==0){
  aveF11 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveF11) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events", "Feed AUC (a.u^2)","Genotype","Tastant")
}

#grouping Licking events- Well 11
j <- 1
intervalL11 <- NULL
if (length(licks11$licks)>0){
licks11$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks11$licks)-1)){
  if (length(licks11$licks)==0){
    break
  }
  else if ((licks11[i+1,2])-(licks11[i,2])<5){
    licks11$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks11$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks11$licks[i+1]-licks11$licks[i]
    intervalL11 <- rbind(intervalL11,interval) 
  }
}
auclick11 <- split(licks11,licks11$EventNumber)
auclick11 <- auclick11[sapply(auclick11, function(d) max(d[[3]]) < mode11+100)]
licks11 <- do.call("rbind",auclick11)
licks11 <- licks11[order(licks11$licks),]
j <- 1
intervalL11 <- NULL
if (length(licks11$licks)>0){
  licks11$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks11$licks)-1)){
  if (length(licks11$licks)==0){
    break
  }
  else if ((licks11[i+1,2])-(licks11[i,2])<5){
    licks11$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks11$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks11$licks[i+1]-licks11$licks[i]
    intervalL11 <- rbind(intervalL11,interval) 
  }
}

#get rid of events only 1 second long
if (length(licks11$licks)>0){
licks11 <- group_by(licks11, EventNumber) %>% filter(n()>1)
#calculating mean
aveL11 <-aggregate(licks11[, 3], list(licks11$EventNumber), max) 
aveL11 <-cbind(aveL11, aggregate(licks11[, 3], list(licks11$EventNumber), length)) 
aveL11 <- aveL11[,-3]
aveL11 <- cbind(aveL11, mean(aveL11[,3]), mean(licks11$datum.licks.),sum(aveL11[,3]), mean(intervalL11), 
                auc((1:length(licks11$licks)),licks11$datum.licks.),genotype6,taste5)
colnames(aveL11) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                      "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                      "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
aveL11 <- aveL11[mixedorder(aveL11$`Event Number`),]}
if (length(licks11$licks)==0){
  aveL11 <- data.frame(matrix(ncol=10, nrow=1))
  colnames(aveL11) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events", "Lick AUC (a.u^2)","Genotype","Tastant")
}
if (is.null(intervalF11)==TRUE){
  intervalF11 <- NA
}
if (length(feeds11)>0){
  aveW11 <- cbind(aveF11$`Mean Length of all Events`[1],aveL11$`Mean Length of all Events`[1],mean(aveF11$`Mean Intensity(per event)`),mean(aveL11$`Mean Intensity(per event)`),aveF11$`Total time (All Events)`[1],aveL11$`Total time (All Events)`[1],
                  aveF11$`Mean interval between Events`[1],aveL11$`Mean interval between Events`[1], 
                  length(aveF11$`Event Number`),length(aveL11$`Event Number`),
                  aveF11$`Feed AUC (a.u^2)`, aveL11$`Lick AUC (a.u^2)`,max(w11f1$datum.feeds.),length(w11f1$datum.feeds.), 
                  auc((1:length(w11f1$feeds)),w11f1$datum.feeds.),intervalF11[1],feeds11$feeds[1],
                  NA,
                  genotype6,taste5, runinfo,11)
  aucvec11 <- lapply(aucfeed11, function(d) (auc(d[[2]],d[[3]])))
  aucvec11 <- unlist(aucvec11)
  aucvec11 <- as.vector(aucvec11)
  deltaintw11 <- data.frame(aveF11$`Event Number`, ((aveF11$`Mean Intensity(per event)`
                                                     -mean(aveF11$`Mean Intensity(per event)`))/mean(aveF11$`Mean Intensity(per event)`)),((aveF11$`Length of Event`-mean(aveF11$`Length of Event`))/
                                                                                                                                             mean(aveF11$`Length of Event`)),(aucvec11-mean(aucvec11))/mean(aucvec11),
                            ((intervalF11[1]-mean(intervalF11))/mean(intervalF11)),((intervalF11[2]-mean(intervalF11))/mean(intervalF11)),((intervalF11[3]-mean(intervalF11))/mean(intervalF11)),
                            ((intervalF11[4]-mean(intervalF11))/mean(intervalF11)),((intervalF11[5]-mean(intervalF11))/mean(intervalF11)),
                            genotype6, taste5,paste(runinfo,"Well 11", sep=" "))}

if (length(feeds11)==0){
  aveW11 <- cbind(aveF11$`Mean Length of all Events`[1],aveL11$`Mean Length of all Events`[1],mean(aveF11$`Mean Intensity(per event)`),
                  mean(aveL11$`Mean Intensity(per event)`),aveF11$`Total time (All Events)`[1],aveL11$`Total time (All Events)`[1],
                  aveF11$`Mean interval between Events`[1],aveL11$`Mean interval between Events`[1], 
                  length(aveF11$`Event Number`),length(aveL11$`Event Number`),
                  aveF11$`Feed AUC (a.u^2)`, aveL11$`Lick AUC (a.u^2)`,0,0, 
                  0,0,NA,NA,genotype6,taste5, runinfo,11)  
  deltaintw11 <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype6,taste5,paste(runinfo,"Well 11", sep=" "))
}
colnames(deltaintw11) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                           "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
colnames(aveW11) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                      "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                      "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                      "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                      "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")

###WELL 12 CHUNK
#grouping feeding events- Well12
j <- 1
intervalF12 <- NULL
if (length(feeds12$feeds)>0){
feeds12$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds12$feeds)-1)){
  if (length(feeds12$feeds)==0){
    break
  }
  else if ((feeds12[i+1,2])-(feeds12[i,2])<5){
    feeds12$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds12$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds12$feeds[i+1]-feeds12$feeds[i]
    intervalF12 <- rbind(intervalF12,interval)
  }
}
aucfeed12 <- split(feeds12,feeds12$EventNumber)
aucfeed12 <- aucfeed12[sapply(aucfeed12, function(d) max(d[[3]]) >mode12+ 100)]
feeds12 <- do.call("rbind",aucfeed12)
feeds12 <- feeds12[order(feeds12$feeds),]
j <- 1
intervalF12 <- NULL
if (length(feeds12$feeds)>0){
  feeds12$EventNumber <- paste("Feeding Event Number: ",j,sep="")}
for (i in 1:(length(feeds12$feeds)-1)){
  if (length(feeds12$feeds)==0){
    break
  }
  else if ((feeds12[i+1,2])-(feeds12[i,2])<5){
    feeds12$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    feeds12$EventNumber[i+1] <- paste("Feeding Event Number: ",j,sep="")
    interval <-  feeds12$feeds[i+1]-feeds12$feeds[i]
    intervalF12 <- rbind(intervalF12,interval)
  }
}

#get rid of events only 1 second long
if (length(feeds12$feeds)>0){
feeds12 <- group_by(feeds12, EventNumber) %>% filter(n()>1)
#calculating mean
aveF12 <-aggregate(feeds12[, 3], list(feeds12$EventNumber), max) 
aveF12 <-cbind(aveF12, aggregate(feeds12[, 3], list(feeds12$EventNumber), length)) 
aveF12 <- aveF12[,-3]
w12f1 <- feeds12[feeds12$EventNumber %in% c("Feeding Event Number: 1"), ]
aveF12 <- cbind(aveF12, mean(aveF12[,3]), mean(feeds12$datum.feeds.),sum(aveF12[,3]), mean(intervalF12), 
                auc((1:length(feeds12$feeds)),feeds12$datum.feeds.),genotype6,taste6)
colnames(aveF12) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                     "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                     "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
aveF12 <- aveF12[mixedorder(aveF12$`Event Number`),]}
if (length(feeds12$feeds)==0){
  aveF12 <- data.frame(matrix(ncol=10, nrow=1))
  colnames(aveF12) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)",
                        "Mean interval between Events","Feed AUC (a.u^2)", "Genotype","Tastant")
}

#grouping Licking events- Well 12
j <- 1
intervalL12 <- NULL
if (length(licks12$licks)>0){
licks12$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks12$licks)-1)){
  if (length(licks12$licks)==0){
    break
  }
  if ((licks12[i+1,2])-(licks12[i,2])<5){
    licks12$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks12$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks12$licks[i+1]-licks12$licks[i]
    intervalL12 <- rbind(intervalL12,interval)
  }
}
auclick12 <- split(licks12,licks12$EventNumber)
auclick12 <- auclick12[sapply(auclick12, function(d) max(d[[3]]) < mode12+100)]
licks12 <- do.call("rbind",auclick12)
licks12 <- licks12[order(licks12$licks),]
j <- 1
intervalL12 <- NULL
if (length(licks12$licks)>0){
  licks12$EventNumber <- paste("Licking Event Number: ",j,sep="")}
for (i in 1:(length(licks12$licks)-1)){
  if (length(licks12$licks)==0){
    break
  }
  if ((licks12[i+1,2])-(licks12[i,2])<5){
    licks12$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
  }
  else{
    j <- j+1
    licks12$EventNumber[i+1] <- paste("Licking Event Number: ",j,sep="")
    interval <-  licks12$licks[i+1]-licks12$licks[i]
    intervalL12 <- rbind(intervalL12,interval)
  }
}

#get rid of events only 1 second long
if (length(licks12$licks)){
licks12 <- group_by(licks12, EventNumber) %>% filter(n()>1)
#calculating mean
aveL12 <-aggregate(licks12[, 3], list(licks12$EventNumber), max) 
aveL12 <-cbind(aveL12, aggregate(licks12[, 3], list(licks12$EventNumber), length)) 
aveL12 <- aveL12[,-3]
aveL12 <- cbind(aveL12, mean(aveL12[,3]), mean(licks12$datum.licks.),sum(aveL12[,3]), mean(intervalL12), 
                auc((1:length(licks12$licks)),licks12$datum.licks.),genotype6,taste6)
colnames(aveL12) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                     "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)", 
                     "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
aveL12 <- aveL12[mixedorder(aveL12$`Event Number`),]}
if (length(licks12$licks)==0){
  aveL12 <- data.frame(matrix(ncol=10,nrow=1))
  colnames(aveL12) <- c("Event Number", "Mean Intensity(per event)","Length of Event", 
                        "Mean Length of all Events", "Mean Intensity (all Events)", "Total time (All Events)", 
                        "Mean interval between Events","Lick AUC (a.u^2)", "Genotype","Tastant")
}
if (is.null(intervalF12)==TRUE){
  intervalF12 <- NA
}
if (length(feeds12)>0){
  aveW12 <- cbind(aveF12$`Mean Length of all Events`[1],aveL12$`Mean Length of all Events`[1],
                  mean(aveF12$`Mean Intensity(per event)`),mean(aveL12$`Mean Intensity(per event)`),
                  aveF12$`Total time (All Events)`[1],aveL12$`Total time (All Events)`[1],aveF12$`Mean interval between Events`[1],
                  aveL12$`Mean interval between Events`[1], length(aveF12$`Event Number`),length(aveL12$`Event Number`),
                  aveF12$`Feed AUC (a.u^2)`, aveL12$`Lick AUC (a.u^2)`,max(w12f1$datum.feeds.),length(w12f1$datum.feeds.), 
                  auc((1:length(w12f1$feeds)),w12f1$datum.feeds.),intervalF12[1], feeds12$feeds[1],
                  sum((aggregate(licks12[, 2], list(licks12$EventNumber), max) < feeds12$feeds[1]),na.rm=TRUE),
                  genotype6, taste6, runinfo,12)
  aucvec12 <- lapply(aucfeed12, function(d) (auc(d[[2]],d[[3]])))
  aucvec12 <- unlist(aucvec12)
  aucvec12 <- as.vector(aucvec12)
  deltaintw12 <- data.frame(aveF12$`Event Number`, ((aveF12$`Mean Intensity(per event)`
                                                     -mean(aveF12$`Mean Intensity(per event)`))/mean(aveF12$`Mean Intensity(per event)`)),((aveF12$`Length of Event`-mean(aveF12$`Length of Event`))/
                                                                                                                                             mean(aveF12$`Length of Event`)),(aucvec12-mean(aucvec12))/mean(aucvec12),
                            ((intervalF12[1]-mean(intervalF12))/mean(intervalF12)),((intervalF12[2]-mean(intervalF12))/mean(intervalF12)),((intervalF12[3]-mean(intervalF12))/mean(intervalF12)),
                            ((intervalF12[4]-mean(intervalF12))/mean(intervalF12)),((intervalF12[5]-mean(intervalF12))/mean(intervalF12)),
                            genotype6, taste6,paste(runinfo,"Well 12", sep=" "))}
if (length(feeds12)==0){
  aveW12 <- cbind(aveF12$`Mean Length of all Events`[1],aveL12$`Mean Length of all Events`[1],
                  mean(aveF12$`Mean Intensity (all Events)`),mean(aveL12$`Mean Intensity (all Events)`),
                  aveF12$`Total time (All Events)`[1],aveL12$`Total time (All Events)`[1],aveF12$`Mean interval between Events`[1],
                  aveL12$`Mean interval between Events`[1], length(aveF12$`Event Number`),length(aveL12$`Event Number`),
                  aveF12$`Feed AUC (a.u^2)`, aveL12$`Lick AUC (a.u^2)`,0,0, 
                  0,0, NA,NA, genotype6, taste6, runinfo,12)  
  deltaintw12 <-  data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,genotype6,taste6,paste(runinfo,"Well 12", sep=" "))
}
colnames(deltaintw12) <- c('Event.Number','Delta.Peak.Intensity', 'Delta.Duration', 'Delta.AUC', "Delta interval event 1", "Delta interval event 2","Delta interval event 3",
                           "Delta interval event 4","Delta interval event 5",'Genotype','Tastant','Run.Well.ID')
                       
colnames(aveW12) <- c("Mean Feeding Length","Mean Tasting Length","Mean Peak Feeding Intensity","Mean Peak Tasting Intensity",
                      "Total Feeding Time","Total Tasting Time","Interval b/w Feeding", "Interval b/w Tasting", 
                      "Feeding Events","Tasting Events",  "Feed AUC (a.u^2)", "Lick AUC (a.u^2)",
                      "First Feed Event Intensity", "First Feed Event Duration", "First Feeding Event AUC",
                      "Interval bw Feed#1/#2","First Feed timestamp","#Contacts before 1st Feeding Event", "Genotype", "Tastant", "Run", "Well")
#masterlistdelta <- smartbind(masterlistdelta,deltaintw1,deltaintw2,deltaintw3,deltaintw4,deltaintw5,deltaintw6,deltaintw7,
                            # deltaintw8,deltaintw9,deltaintw10,deltaintw11,deltaintw12)
#write.csv(masterlistdelta, "/Users/Vaibhav/Desktop/FLIC FINAL DATA REPO updated AUC/masterlistdelta.csv")
deltastats <- rbind(deltaintw1,deltaintw2,deltaintw3,deltaintw4,deltaintw5,deltaintw6,deltaintw7,deltaintw8,deltaintw9,deltaintw10,deltaintw11,deltaintw12)
write.csv(deltastats,"Delta Stats.csv")
platestats <- rbind(aveW1[1,],aveW2[1,],aveW3[1,],aveW4[1,],aveW5[1,],aveW6[1,],aveW7[1,],aveW8[1,],aveW9[1,],aveW10[1,],aveW11[1,],aveW12[1,])  
write.csv(platestats, "Plate Stats.csv")  
  