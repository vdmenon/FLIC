library(pROC)
library(ROCR)
library(ggplot2)
#load the training data
training.data.raw <- read.csv("/Users/Vaibhav/Desktop/FLIC FINAL DATA REPO updated AUC/CSHL test set.csv",stringsAsFactors=FALSE)
#Handling NA's and converting values in feature columns to numeric objects 
#Good column
training.data.raw$Good <- as.factor(training.data.raw$Good)
#mean feeding length
training.data.raw$Mean.Feeding.Length <- as.numeric(training.data.raw$Mean.Feeding.Length)
training.data.raw$Mean.Feeding.Length[is.na(training.data.raw$Mean.Feeding.Length)] <- 0
#mean tasting legnth
training.data.raw$Mean.Tasting.Length <- as.numeric(training.data.raw$Mean.Tasting.Length)
training.data.raw$Mean.Tasting.Length[is.na(training.data.raw$Mean.Tasting.Length)] <- 0
#Mean Feeding intensity
training.data.raw$Mean.Peak.Feeding.Intensity[is.na(training.data.raw$Mean.Peak.Feeding.Intensity)] <- 0
training.data.raw$Mean.Peak.Feeding.Intensity <- as.numeric(training.data.raw$Mean.Peak.Feeding.Intensity)
#Mean Tasting Intensity
training.data.raw$Mean.Peak.Tasting.Intensity[is.na(training.data.raw$Mean.Peak.Tasting.Intensity)] <- 0
training.data.raw$Mean.Peak.Tasting.Intensity <- as.numeric(training.data.raw$Mean.Peak.Tasting.Intensity)
#Total Feeding Time
training.data.raw$Total.Feeding.Time[is.na(training.data.raw$Total.Feeding.Time)] <- 0
training.data.raw$Total.Feeding.Time <- as.numeric(training.data.raw$Total.Feeding.Time)
#Total Tasting Time
training.data.raw$Total.Tasting.Time[is.na(training.data.raw$Total.Tasting.Time)] <- 0
training.data.raw$Total.Tasting.Time <- as.numeric(training.data.raw$Total.Tasting.Time)
#Interval between Feeding
training.data.raw$Interval.b.w.Feeding[is.na(training.data.raw$Interval.b.w.Feeding)] <- 0
training.data.raw$Interval.b.w.Feeding <- as.numeric(training.data.raw$Interval.b.w.Feeding)
#Interval between Tasting
training.data.raw$Interval.b.w.Tasting[is.na(training.data.raw$Interval.b.w.Tasting)] <- 0
training.data.raw$Interval.b.w.Tasting <- as.numeric(training.data.raw$Interval.b.w.Tasting)
#Number of feeding events
training.data.raw$Feeding.Events[is.na(training.data.raw$Feeding.Events)] <- 0
training.data.raw$Feeding.Events <- as.numeric(training.data.raw$Feeding.Events)
#Number of tasting events
training.data.raw$Tasting.Events[is.na(training.data.raw$Tasting.Events)] <- 0
training.data.raw$Tasting.Events <- as.numeric(training.data.raw$Tasting.Events)
#Feed Area under curve
training.data.raw$Feed.AUC..a.u.2. [is.na(training.data.raw$Feed.AUC..a.u.2.)] <- 0
training.data.raw$Feed.AUC..a.u.2. <- as.numeric(training.data.raw$Feed.AUC..a.u.2.)
#Feed Area under curve
training.data.raw$Lick.AUC..a.u.2. [is.na(training.data.raw$Lick.AUC..a.u.2.)] <- 0
training.data.raw$Lick.AUC..a.u.2. <- as.numeric(training.data.raw$Lick.AUC..a.u.2.)

#Generate split sets of equal size for train and test set
smp_size <- floor(0.5*nrow(training.data.raw))
set.seed(1234) #reproducibility
training_shuffle <- sample(seq_len(nrow(training.data.raw)),size=smp_size)
trainA <- training.data.raw[training_shuffle,]
testA <- training.data.raw[-training_shuffle,]

##creating model
model <- glm(Good ~ Mean.Feeding.Length+Mean.Tasting.Length+Mean.Peak.Feeding.Intensity+Mean.Peak.Tasting.Intensity+Total.Feeding.Time+Total.Tasting.Time+Interval.b.w.Feeding+Interval.b.w.Tasting+Feeding.Events+Tasting.Events,family=binomial(link='logit'),data=trainA,maxit = 150)
summary(model)
#analysis of deviance
anova(model)

#predict value between 0 and 1
fitted.results <- predict(model,newdata=subset(testA,select=c(1,2,3,4,5,6,7,8,9,10)),type='response')
testA$ModelResult=fitted.results
g <- roc(Good ~ModelResult, data=testA, print.auc=TRUE)
rocplot <- ggroc(g)
rocplot
rocplot + xlab("False Positive Rate (Specificity)")+ ylab("True Positive Rate (Sensitivity)") +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0), color="darkgrey", linetype="dashed")
print(auc(g))

#threshold of predictved value=0.5 separating good from bad
fitted.results <- ifelse(fitted.results > 0.5,"g","b")
testA$ModelResult=fitted.results
misClasificError <- mean(testA$ModelResult != testA$Good)
print(paste('Accuracy',1-misClasificError))
