# 450K test zip file based import and normalisation test code
# Matthew Bashton

#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("minfi")
#biocLite("IlluminaHumanMethylation450kmanifest")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#biocLite("FDb.InfiniumMethylation.hg19")
#biocLite("IlluminaHumanMethylationEPICmanifest")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#install.packages("RColorBrewer")
#biocLite("CopyNumber450kData")
#biocLite("conumee")

# Below needs to start before reactive code
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(RColorBrewer)
library(CopyNumber450kData)
library(conumee)

# Test code
setwd("~/Downloads/test/")

# Unzip a zip
unzip("BreakingBad.zip", exdir = ".", overwrite = TRUE, junkpaths = TRUE)

# Workout how many lines to skip on csv
lines <- readLines("~/Downloads/test/sampleSheet.csv")
lines_to_skip <- grep("Sample_Name", lines) - 1

# Read in csv
targets <- read.csv("sampleSheet.csv", skip=lines_to_skip,header= TRUE, sep=",")

# Use info from the sampleSheet to compose file names which are made up of
# Sentrix ID and position
targets$Basename <- paste(targets$Sentrix_ID,"_",targets$Sentrix_Position,sep="")

# Read in 450k data from idat files to create a RGset both baseDir and targets 
# are set above, this can take some time for large datasets; RG set has data
# from both read and green channels of array.
RGset <- read.metharray.exp(base = "~/Downloads/test/", targets=targets, verbose = TRUE)

# Store phenotype data by extracting it with pData accessor
pd <- pData(RGset)

# Density plots for QC in R session
densityPlot(RGset, main = "Beta value distribution", xlab = "Beta value")
densityBeanPlot(RGset, sampNames = pd$Sample_Name)

# Normalise all they arrays 
# Noob is not fast, but is faster than SWAN!  illumina is very fast!
Mset <- preprocessIllumina(rgSet = RGset, bg.correct = TRUE)

# Extract probes we need for the classifyer
load("Entire_10000_June2015.RData")

# Extract beta values into numeric matrix, cols are samples rows are probes
hov.final <- getBeta(Mset)

# check the probe names, choose only those 10,000 probes which are matched with the reference (n=434 Classifier) 
hov.10 <- hov.final[rownames(hov.final) %in% rownames(disc.10),]

# Change col names for those in pd cols should be same order as rows in pd
colnames(hov.10) <- pd$Sample_Name

# 850k to 450k downgrade
# Sample names on Mset
colnames(Mset) <- pd$Sample_Name
load("Probes450k.Rdata")
Mset450klike <- Mset[rownames(Mset) %in% Probes450k,]
EPIC_overlap <- rownames(Mset450klike)

# Conumee data pre load this
data(exclude_regions)
data(detail_regions)
head(detail_regions, n = 2)
anno <- CNV.create_anno(exclude_regions = exclude_regions, detail_regions = detail_regions)
minfi.data <- CNV.load(Mset450klike)
# Get control data
load("~/MAC/Mset.CB.RData")

# Subset controls so they only have probes which are also on EPIC not 450k exclusive
Mset.CB <- Mset.CB[rownames(Mset.CB) %in% EPIC_overlap, ]
controls.data <-  CNV.load(Mset.CB)

i <- 1
x <- CNV.fit(minfi.data[i], controls.data, anno)
