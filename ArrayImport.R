# 450K test zip file based import and normalisation test code
# Matthew Bashton

#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("minfi")
#biocLite("IlluminaHumanMethylation450kmanifest")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#biocLite("IlluminaHumanMethylationEPICmanifest")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#install.packages("RColorBrewer")

# Below needs to start before reactive code
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(RColorBrewer)

# Test code
setwd("~/MAC")

# Unzip a zip
unzip("test_set.zip", exdir = ".", overwrite = TRUE, junkpaths = TRUE)

# Workout how many lines to skip on csv
lines <- readLines("sampleSheet.csv")
lines_to_skip <- grep("Sample_Name", lines) - 1

# Read in csv
targets <- read.csv("sampleSheet.csv", skip=lines_to_skip,header= TRUE, sep=",")

# Use info from the sampleSheet to compose file names which are made up of
# Sentrix ID and position
targets$Basename <- paste(targets$Sentrix_ID,"_",targets$Sentrix_Position,sep="")

# Read in 450k data from idat files to create a RGset both baseDir and targets 
# are set above, this can take some time for large datasets; RG set has data
# from both read and green channels of array.
RGset <- read.metharray.exp(targets=targets, verbose = TRUE)

# Store phenotype data by extracting it with pData accessor
pd <- pData(RGset)

# Density plots for QC in R session
densityPlot(RGset, main = "Beta value distribution", xlab = "Beta value")
densityBeanPlot(RGset, sampNames = pd$Sample_Name)

# Normalise all they arrays 
# Noob is not fast, but is faster than SWAN!
Mset <- preprocessNoob(rgSet = RGset, dyeCorr = TRUE, verbose = TRUE)

# Extract probes we need for the classifyer
load("Entire_10000_June2015.RData")

# Extract beta values into numeric matrix, cols are samples rows are probes
hov.final <- getBeta(Mset)

# check the probe names, choose only those 10,000 probes which are matched with the reference (n=434 Classifier) 
hov.10 <- hov.final[rownames(hov.final) %in% rownames(disc.10),]

# Change col names for those in pd cols should be same order as rows in pd
colnames(hov.10) <- pd$Sample_Name
