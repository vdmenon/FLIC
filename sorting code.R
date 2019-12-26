#Contact: Vaibhav Menon (vmeno002@ucr.edu)
#The purpose of this script is to generate a code that quantifies feeding events and licking events. The
#skeleton of the workhorse function here, feedslicksplots, is a heavily adapted version of Edric Tam's
#makePlotsandTables.R large function, used in the Elife paper Joseph et al, 2017. 
#The basis for separation of feeding from licking is through an intensity threshold (feeding > 100 A.U. 
#above baseline, tasting < 100 A.U. above baseline) that is consistent with other papers using FLIC assays.
#This code can definitely be implemented more efficiently, but it extracted the information we were interested in.
#I have included below the code used to extract information from Well 1 of each plate. One can use looping methods
#and adjust well-specific variables (tastant,genotype) to apply this code to all 12 wells.

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
if (length(feeds1$feeds)==0){ #if there are no feeding events in well
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
