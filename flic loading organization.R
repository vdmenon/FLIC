#This script reads through an accompanying text file that contains annotative information about flies being tested. Includes
#info such as genotype, tastant, date, run number, etc. Following the format of that text file critical for this script to work properly,
#as it relies on specific line numbers.
#description= "/Users/Vaibhav/Desktop/grad school/year1/winter 2018/dahanukar rotation/FLIC data sorting code/DFM1/20190129-1 description.txt"
library(stringr)
descriptionfile=file.choose()
description <- file.path(descriptionfile)

#Pulling genotype information for each well from the description text file
genotype1line <- readLines(description, n=2)
genotype1 <- grep("1:", genotype1line, value =TRUE)
genotype1 <- gsub("1: ","",genotype1)

genotype2line <- readLines(description, n=3)
genotype2 <- grep("2:", genotype2line, value =TRUE)
genotype2 <- gsub("2: ","", genotype2)

genotype3line <- readLines(description, n=4)
genotype3 <- grep("3:", genotype3line, value =TRUE)
genotype3 <- gsub("3: ","",genotype3)

genotype4line <- readLines(description, n=5)
genotype4 <- grep("4:", genotype4line, value =TRUE)
genotype4 <- gsub("4: ","",genotype4)

genotype5line <- readLines(description, n=6)
genotype5 <- grep("5:", genotype5line, value =TRUE)
genotype5 <- gsub("5: ","",genotype5)

genotype6line <- readLines(description, n=7)
genotype6 <- grep("6:", genotype6line, value =TRUE)
genotype6 <- gsub("6: ","",genotype6)
#Labeling A wells
wellnumbersA <- c(1,3,5,7,9,11)
genecharvectorA<- c(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6)

#Labeling B wells
genecharvectorB<- c(genotype1, genotype2, genotype3, genotype4, genotype5, genotype6)
wellnumbersB <- c(2,4,6,8,10,12)
genedf <- data.frame(genecharvectorA,genecharvectorB)

#copying tastant information
#DFM1
tastant1 <- read.table(description, header = FALSE, sep = "_", skip = 10, nrows = 1)
taste1 <- grep("A: ", tastant1$V1, value =TRUE)
taste1 <- gsub("A: ","", taste1)

tastant2 <- read.table(description, header = FALSE, sep = "_", skip = 11, nrows = 1)
taste2 <- grep("B: ", tastant2$V1, value =TRUE)
taste2 <- gsub("B: ","", taste2)

#DFM2
tastant3 <- read.table(description, header = FALSE, sep = "_", skip = 14, nrows = 1)
taste3 <- grep("A: ", tastant3$V1, value =TRUE)
taste3 <- gsub("A: ","", taste3)

tastant4 <- read.table(description, header = FALSE, sep = "_", skip = 15, nrows = 1)
taste4 <- grep("B: ", tastant4$V1, value =TRUE)
taste4 <- gsub("B: ","", taste4)

#DFM3
tastant5 <- read.table(description, header = FALSE, sep = "_", skip = 18, nrows = 1) 
taste5 <- grep("A: ", tastant5$V1, value =TRUE)
taste5 <- gsub("A: ","", taste5)

tastant6 <- read.table(description, header = FALSE, sep = "_", skip = 19, nrows = 1)
taste6 <- grep("B: ", tastant6$V1, value =TRUE)
taste6 <- gsub("B: ","", taste6)

tastantvectorDFM1A <-c(taste1)
tastantvectorDFM1B <-c(taste2)
tastantvectorDFM2A <-c(taste3)
tastantvectorDFM2B <-c(taste4)
tastantvectorDFM3A <-c(taste5)
tastantvectorDFM3B <-c(taste6)
