# Methylation Array Classifier (MAC)
# Classifier code be Dr Reza Rafiee, 2015-2016
# Shiny web code Dr Matthew Bashton, 2016
# 450K Classifer, Software Version 1.3.2 (Successful version)
# Baed on NMF projection and SVM algorithm

# Input: NMB samples include: beta values of 10,000 probes from 450K methylation
# profiling

# Output: Classifier confidence and subgroup labels for input samples 4-group
# classifier (WNT,SHH, Grp3 and Grp4) - 4 metagenes

# For display later  
AppName <- "Methylation Array Classifier (MAC)"
AppVersion <- "1.3.2"

# Load librarys
library(e1071)
library(parallel)

## Load code for metagene projection
source("NMF_Functions_minfi.R")

## Get same 10000 probes that went into 434 classifier
## Load 10000 probes
load("Entire_10000_June2015.RData")

# Obtained by MME fit distribution from the final probability estimates of
# validation cohort (n=276 Hoves.), Jan. 2016

## CHANGE FOR RELEASE
#GoldCohort.Threshold1 <- 0.7045746  
GoldCohort.Threshold1 <- 0.72

## 220 Training set - 22 December 2015/Updated 14 January 2016
Trainingset450k_4Metagene_WithSubgroup <- as.matrix(read.csv("220TrainingCohort450KSubgrouping_4MetagenesANDSubgroupLabels_14Jan2016.csv",header=T,row.names=1))
# Remove subgroup col
Trainingset450k4M <- Trainingset450k_4Metagene_WithSubgroup[,1:4]
# Extract 220 subgroup labels
labels220 <- as.character(Trainingset450k_4Metagene_WithSubgroup[,5])
# Convert to factor
subgroup_labels <- factor(labels220)
y1 <- subgroup_labels
#-----------------------------------------------------------------

## Give names to groups
grp.4 <- ifelse(y1 == 1, "WNT", ifelse(y1 == 2,"SHH", ifelse(y1== 3 ,"Grp3","Grp4")))
## Load W matrix
bar <- get(load('2goldStandard_W_EntireCohort_10000.RData'))
avgW.4 <- bar[[4]]
Groups <- factor(grp.4, levels=c("WNT","SHH","Grp3","Grp4"))
trainH <- Trainingset450k4M 


## Preparing Input Dataset
# Input dataset: 15 samples of Volker Hovstetadt dataset
tmp <- read.csv("www/GSE54880_test_set.csv")
hov.final <- as.matrix(tmp[,-1])
rownames(hov.final) <- tmp[,1]

# check the probe names, choose only those 10,000 probes which are matched with the reference (n=434 Classifier) 
hov.10 <- hov.final[rownames(hov.final) %in% rownames(disc.10),]

# Match features
hov.match <- DW.MP.Match.and.Select(disc.10, hov.10)
## Match W matrix
hov.H <- DW.MP.Factors.Project.C(hov.10, avgW.4)  # applying the Moore-Penrose pseduoinverse of Wm (avgw.4) to the input data.


## Creating the classifier model using SVM
Optimised_cost <- 1.4   
Optimised_Gamma <- 0.02  

## Further analysis to assess confidence of calls
x <- 1000 ## Number of iterations

train.beta <- trainH 
amount <- round(0.80*nrow(train.beta)) 

sel2<- lapply(1:x, function(i) {
  set.seed(i)
  sample(1:nrow(train.beta), amount, replace=F) #nrow(t(train.beta))
})

## MB this bit causes a delay
Radial.svms <- mclapply(1:x, 
                        mc.cores=4,
                        function(i)  svm(x = trainH[sel2[[i]],],  #t(trainH)[sel2[[i]],],
                                         y = Groups[sel2[[i]]], scale = F,
                                         tolerance = 0.00001, type = "C-classification",
                                         kernel = "radial",cost = Optimised_cost,
                                         gamma=Optimised_Gamma, probability = T,
                                         seed=123456) 
)


## Test on input samples (hov)
Radial.tests <- mclapply(1:x,
                         mc.cores=4,
                         function(i) predict(Radial.svms[[i]],
                                             newdata=t(hov.H), # 4 Metagenes of input dataset 
                                             decision.values = T,
                                             probability = T)
)
prob.test <- (lapply(1:x,
                     function(i) attr(Radial.tests[[i]], "probabilities"))
)

####################################### Creating Pobes2 #############################################
k <- FALSE

for (j in 1:x) # the number of iterations
{
  predProbTemp <-prob.test[[j]] # j iteration
  predProbTemp <- t(predProbTemp)
  predProbTemp <- predProbTemp[c("WNT", "SHH", "Grp3","Grp4"),] # order the matrix 
  predProbTemp <- t(predProbTemp)

  if (k == FALSE) # Making defult tables
  {
    
    predProbabilities <- matrix(ncol = 4, nrow =nrow(predProbTemp)*x, 0.0)
    predProbabilities <- predProbTemp 
    k <- TRUE
  }
  else
  {
    predProbabilities <- rbind(predProbabilities,predProbTemp)
  }
}


probs2 <- matrix(ncol=nrow(predProbTemp),nrow=x,0.0)
colnames(probs2) <- rownames(predProbTemp)

for (ttt in 1:nrow(predProbTemp)) # number of samples
{
  mmm <- matrix(ncol = 4, nrow =x, 0.0)
  colnames(mmm) <- c("WNT","SHH", "Grp3", "Grp4")   # exact order of prob.test columns
  
  gg <- 0
  for (fftt in 1:x)
  {
    gg <- gg + 1
    mmm[gg,] <- predProbabilities[ttt+nrow(predProbTemp)*(fftt-1),]
  }
  
  ProbSubgroup <- apply(mmm[,1:4],1,max)
  probs2[,ttt] <- ProbSubgroup
  
}
####################################### creating pobes2 #############################################
# End of "Creating the final model for whole training set with selected parameters"
final.model <- svm(x = trainH, y = factor(Groups), scale=F,  
                   tolerance = 0.00001, type="C-classification", kernel = "radial",
                   probability = T, seed = 123456, cost = Optimised_cost, gamma = Optimised_Gamma)

#Testing with the Hoves... dataset
hov.calls <- predict(final.model, newdata= t(hov.H), probability = TRUE, decision.values = TRUE)
test.pred <- hov.calls

prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
maxProbs <- apply(prob.test,1,max)

Total.No.of.Samples <- ncol(hov.H)

maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("WNT", "SHH", "Grp3", "Grp4"))

## identify the threshold
threshold <- GoldCohort.Threshold1  # obtained by MME fit distribution of final probability estimates, 07 January 2016

# Tmp re-level.
levels(maxProbsWhich) <- c(1,2,3,4)

maxProbsCol <- ifelse(maxProbsWhich==1,"blue",ifelse(maxProbsWhich==2,"red",
                                                     ifelse(maxProbsWhich==3,"yellow2","darkgreen")))

maxProbsCol2 <- ifelse(maxProbsCol=="yellow2","#EEEE0066", ifelse(maxProbsCol=="blue","#0000FF66",
                                                                  ifelse(maxProbsCol=="darkgreen","#00640066","#FF000066")))

# MB Output for classification table 
levels(maxProbsWhich) <- c("WNT", "SHH", "Grp3", "Grp4")
results.df <- data.frame(names(maxProbsWhich), as.character(maxProbsWhich), maxProbs, row.names = NULL, stringsAsFactors = FALSE)
colnames(results.df) <- c("Sample", "Subgroup", "Confidence") 

# New plot code
cat(paste("Removing data points below threshold", threshold, "from graph:\n"))
index <- maxProbs > threshold
cat(names(maxProbs[!index]), "\n")
new.probs2 <- probs2[,index]
new.maxProbs <- maxProbs[index]
new.maxProbsWhich <- maxProbsWhich[index]
new.maxProbsCol <- maxProbsCol[index]
new.maxProbsCol2 <- maxProbsCol2[index]
new.Total.No.of.Samples <- length(maxProbs[index])

par(mfrow=c(1,1))
#par(mar=c(6,4,2,1) + 0.1)
par(mar=c(6,4,4,1) + 0.1)
par(cex=1.3)
par(cex.axis=1)

heading <- paste("Medulloblastoma subgroup call confidence intervals for", new.Total.No.of.Samples, "samples")

boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",new.probs2[,order(new.maxProbsWhich, new.maxProbs)],outpch=NA,ylim=c(0,1),las=2,
        col=new.maxProbsCol2[order(new.maxProbsWhich,new.maxProbs)] )

abline(col="grey",lty = 1, h = threshold)

# How many subgroups of each colour are we plotting
tmp <- table(new.maxProbsCol)
desired_col_order <-c("blue", "red", "yellow2", "darkgreen")
to_sort <- names(tmp)
# Re order by correct sub group col order using match on the desired_col_order vector
tmp <- tmp[to_sort[order(match(to_sort,desired_col_order))]]
# Index of where to draw the sub group deviders via cumsum
grp.sum <- cumsum(tmp)
# Add 0.5 to grp.sum for abline
grp.sum <- grp.sum + 0.5
# Index out final element of grp.sum to get rid of unwanted final abline
grp.sum <- grp.sum[1:length(grp.sum)-1]
# Check
grp.sum
abline(v=grp.sum)
#lines(col="black",lwd=2,new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
points(col=new.maxProbsCol[order(new.maxProbsWhich,new.maxProbs)],pch=19, new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4"), col=c("blue", "red", "yellow2", "darkgreen"), pch=19)   
axis(2, las=2)             

