# Shiny server for Methylation Array Classifier (MAC)
# Classifier code be Dr Reza Rafiee, 2015-2016
# Shiny web code Dr Matthew Bashton, 2016
# 450K Classifer, Software Version 1.3.2 (Successful version)
# Based on NMF projection and SVM algorithm

# Input: A zip file of idats and a sampleSheet.csv (as one would get from array provider)

# Output: Classifier confidence and subgroup labels for input samples 4-group
# classifier (WNT,SHH, Grp3 and Grp4) - 4 metagenes

# For display later  
AppName <- "MAC: Methylation Array Classifier"
AppVersion <- "1.3.2-2.0"

# Load all librarys now to speed up reactive part of code
library(shiny)
library(e1071) #for SVM classifier
library(parallel)  # For mclapply speeds up probability estimation 
library(gtools) # Needed for numerically rather than lexicographically sorted strings
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#library(RColorBrewer)

# Set Max upload size to 512MB!  Should allow ball park 48 samples max (not tested)
options(shiny.maxRequestSize=512*1024^2)

# Threshold for classification 
threshold <- 0.7045746

# Cut-off for percentage of failed probes to reject an array from returning a classification
percent_pval_cutoff <- 0.06

### Get input file name from UI
shinyServer(function(input, output) {
  
  ## Load code for metagene projection
  source("NMF_Functions_minfi.R")
  
  ## Get same 10000 probes that went into 434 classifier
  ## Load 10000 probes
  load("Entire_10000_June2015.RData")
  
  # Obtained by MME fit distribution from the final probability estimates of
  # validation cohort (n=276 Hoves.), Jan. 2016
  
  GoldCohort.Threshold1 <- threshold
  
  ## Preload datasets
  ## 220 Training set - 22 December 2015/Updated 14 January 2016
  Trainingset450k_4Metagene_WithSubgroup <- as.matrix(read.csv("220TrainingCohort450KSubgrouping_4MetagenesANDSubgroupLabels_14Jan2016.csv",header=T,row.names=1))
  # Remove subgroup col
  Trainingset450k4M <- Trainingset450k_4Metagene_WithSubgroup[,1:4]
  # Extract 220 subgroup labels
  labels220 <- as.character(Trainingset450k_4Metagene_WithSubgroup[,5])
  # Convert to factor
  subgroup_labels <- factor(labels220)
  y1 <- subgroup_labels
  
  
  ## Give names to groups
  grp.4 <- ifelse(y1 == 1, "WNT", ifelse(y1 == 2,"SHH", ifelse(y1== 3 ,"Grp3","Grp4")))
  ## Load W matrix
  bar <- get(load('2goldStandard_W_EntireCohort_10000.RData'))
  avgW.4 <- bar[[4]]
  Groups <- factor(grp.4, levels=c("WNT","SHH","Grp3","Grp4"))
  trainH <- Trainingset450k4M 
  
  
  #############################################################################
  ######################## Reactive classifier function #######################
  #############################################################################
  
  classifier <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    cat("input file is here:\n")
    cat(inFile$datapath, "\n")
    
    withProgress(message = 'Processing data', value = 0, {
      
      # Start the clock
      cat("Timing start\n")
      ptm <- proc.time()
      
      ## Get Input Dataset
      # Unzip a zip
      
      # Basedir
      basedir <- dirname(inFile$datapath)
      
      incProgress(0.10, detail = paste("Unziping archive"))
      unzip(inFile$datapath, exdir = basedir, overwrite = TRUE, junkpaths = TRUE)
      
      # Old CSV version
      #tmp <- read.csv(inFile$datapath)
      #hov.final <- as.matrix(tmp[,-1])
      #rownames(hov.final) <- tmp[,1]
      
      # Workout how many lines to skip on csv
      lines <- readLines(paste0(basedir,"/sampleSheet.csv"))
      lines_to_skip <- grep("Sample_Name", lines) - 1
      
      # Read in csv
      incProgress(0.10, detail = paste("Reading sampleSheet.csv"))
      targets <- read.csv(paste0(basedir,"/sampleSheet.csv"), skip=lines_to_skip,header= TRUE, sep=",")
      
      # Use info from the sampleSheet to compose file names which are made up of
      # Sentrix ID and position
      targets$Basename <- paste(targets$Sentrix_ID,"_",targets$Sentrix_Position,sep="")
      
      # Read in 450k data from idat files to create a RGset both baseDir and targets 
      # are set above, this can take some time for large datasets; RG set has data
      # from both read and green channels of array.
      incProgress(0.10, detail = paste("Loading .idat data"))
      RGset <- read.metharray.exp(base = basedir, targets=targets, verbose = TRUE)
      
      # Store phenotype data by extracting it with pData accessor
      pd <- pData(RGset)
      
      # Use this code later
      # Density plots for QC in R session
      #densityPlot(RGset, main = "Beta value distribution", xlab = "Beta value")
      #densityBeanPlot(RGset, sampNames = pd$Sample_Name)
      
      # Normalise all they arrays 
      # Noob is not fast, but is faster than SWAN!
      incProgress(0.10, detail = paste("Normalising arrays"))
      # Noob is slow use Illumina is much faster but need to state ref explicidly
      Mset <- preprocessNoob(rgSet = RGset, dyeCorr = TRUE, verbose = TRUE)
      # Some times produces NAs Do not use!
      # Mset <- preprocessIllumina(rgSet = RGset, bg.correct = TRUE)
      
      # Extract probes we need for the classifyer
      load("Entire_10000_June2015.RData")
      
      # Extract beta values into numeric matrix, cols are samples rows are probes
      hov.final <- getBeta(Mset)
      
      # check the probe names, choose only those 10,000 probes which are matched with the reference (n=434 Classifier) 
      incProgress(0.10, detail = paste("Selecting 10K probes"))
      hov.10 <- hov.final[rownames(hov.final) %in% rownames(disc.10),]
      
      # Change col names for those in pd cols should be same order as rows in pd
      colnames(hov.10) <- pd$Sample_Name
      
      # Work out %age of hov.10 probes with detection p-value above 0.05 
      incProgress(0.10, detail = paste("Acquiring detection p-values"))
      pvals <- detectionP(RGset)
      colnames(pvals) <- pd$Sample_Name
      # Get only our 10k of intrest
      pvals.10 <- pvals[rownames(pvals) %in% rownames(disc.10),]
      TF.pvals.10 <- pvals.10 > 0.05
      TF.pvals <- pvals > 0.05
      # For 10k only
      #percentfail <- round(apply(TF.pvals.10, 2, function(x) mean(x) * 100 ),digits = 2)
      # Fow whole array
      percentfail <- round(apply(TF.pvals, 2, function(x) mean(x) * 100 ),digits = 2)
      
      # Match features
      incProgress(0.10, detail = paste("Match features"))
      hov.match <- DW.MP.Match.and.Select(disc.10, hov.10)
      ## Match W matrix
      incProgress(0.10, detail = paste("Match W matrix"))
      hov.H <- DW.MP.Factors.Project.C(hov.10, avgW.4)  # applying the Moore-Penrose pseduoinverse of Wm (avgw.4) to the input data.
      
      
      ## Creating the classifier model using SVM
      Optimised_cost <- 1.4   
      Optimised_Gamma <- 0.02  
      
      ## Further analysis to assess confidence of calls
      incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 1"))
      x <- 1000 ## Number of iterations
      
      train.beta <- trainH 
      amount <- round(0.80*nrow(train.beta)) 
      
      sel2<- lapply(1:x, function(i) {
        set.seed(i)
        sample(1:nrow(train.beta), amount, replace=F) #nrow(t(train.beta))
      })
      
      ## MB this bit causes a delay
      incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 2"))
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
      incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 3"))
      Radial.tests <- mclapply(1:x,
                               mc.cores=4,
                               function(i) predict(Radial.svms[[i]],
                                                   newdata=t(hov.H), # 4 Metagenes of input dataset 
                                                   decision.values = T,
                                                   probability = T)
      )
      incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 4"))
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
      
      # Let user know we're done
      setProgress(value = 1, message = "Done!")
      cat("\nDone classification\n")
      
    }) # End with progress
    
    # Stop the clock
    time<- (proc.time() - ptm)
    
    # Return all the things:
    classified_data <- list(results.df,
                            Total.No.of.Samples,
                            threshold,
                            probs2,
                            maxProbsWhich,
                            maxProbs,
                            maxProbsCol,
                            maxProbsCol2,
                            time,
                            RGset,
                            pd,
                            percentfail)
    
    names(classified_data) <- c("results.df",
                                "Total.No.of.Samples",
                                "threshold",
                                "probs2",
                                "maxProbsWhich",
                                "maxProbs",
                                "maxProbsCol",
                                "maxProbsCol2",
                                "time",
                                "RGset",
                                "pd",
                                "percentfail")
    
    return(classified_data)
    
    #########################################################################
    ################## End of classifier reactive function ##################
    #########################################################################
    
  }) # End reactive classifier function
  
  
  # Output classification_table #####
  output$classification_table <- renderDataTable({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    results.df <- classified_data$results.df
    # Change name to Subgroup Call of 2nd col
    colnames(results.df)[2] <- "Subgroup Call"
    
    # Add percent failed column 
    percentfail <- classified_data$percentfail
    results.df[,4] <- percentfail 
    colnames(results.df)[4] <- "% failed probes"
    
    # Add Array QC col
    results.df[,5] <- "Pass"
    colnames(results.df)[5] <- "Array QC"
     
    # Apply threshold and label samples as unclassifiable 
    thresholded_results.df <- results.df
    i <- 1
    for (i in 1:nrow(results.df)) {
      if (!results.df[i,"Confidence"] > threshold) {
        thresholded_results.df[i,"Subgroup Call"] <- "Unclassifiable" 
        thresholded_results.df[i,"Confidence"] <- NA
        thresholded_results.df[i,"Array QC"] <- "Pass"
      }
    }
    
    # Apply percent_pval_cutoff to percentfail
    i <- 1
    for (i in 1:nrow(results.df)) {
      if (results.df[i,"% failed probes"] > percent_pval_cutoff) {
        thresholded_results.df[i,"Subgroup Call"] <- "-" 
        thresholded_results.df[i,"Confidence"] <- NA
        thresholded_results.df[i,"Array QC"] <- "Fail"
      }
    }
    
    # Convert to percentage (because medics)
    thresholded_results.df[,3] <- as.character(as.numeric(thresholded_results.df[,3])*100)
    colnames(thresholded_results.df)[3] <- "Probability %"
    thresholded_results.df[is.na(thresholded_results.df)] <- "-"
    # Sort via sample ID (correctly)
    thresholded_results.df <- thresholded_results.df[mixedorder(thresholded_results.df[,1]),]
    return(thresholded_results.df)
    
  })
  
  output$downloadClassification <- downloadHandler(
    filename = "MB_classification.csv",
    content =  function(file) {
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      results.df <- classified_data$results.df
      # Change name to Subgroup Call of 2nd col
      colnames(results.df)[2] <- "Subgroup Call"
    
      # Add percent failed column 
      percentfail <- classified_data$percentfail
      results.df[,4] <- percentfail 
      colnames(results.df)[4] <- "% failed probes"
      
      # Add Array QC col
      results.df[,5] <- "Pass"
      colnames(results.df)[5] <- "Array QC"
      
      # Apply threshold and label samples as unclassifiable 
      thresholded_results.df <- results.df
      for (i in 1:nrow(results.df)) {
        if (!results.df[i,"Confidence"] > threshold) {
          thresholded_results.df[i,"Subgroup Call"] <- "Unclassifiable" 
          thresholded_results.df[i,"Confidence"] <- NA
          thresholded_results.df[i,"Array QC"] <- "Pass"
        }
      }
      
      # Apply percent_pval_cutoff to percentfail
      for (i in 1:nrow(results.df)) {
        if (results.df[i,"% failed probes"] > percent_pval_cutoff) {
          thresholded_results.df[i,"Subgroup Call"] <- "-" 
          thresholded_results.df[i,"Confidence"] <- NA
          thresholded_results.df[i,"Array QC"] <- "Fail"
        }
      }
      
      # Convert to percentage (because medics)
      thresholded_results.df[,3] <- as.character(as.numeric(thresholded_results.df[,3])*100)
      thresholded_results.df[is.na(thresholded_results.df)] <- "-"
      # Sort via sample ID (correctly)
      thresholded_results.df <- thresholded_results.df[mixedorder(thresholded_results.df[,1]),]
      colnames(thresholded_results.df)[3] <- "Probability %"
      write.csv(thresholded_results.df, file, row.names = FALSE) 
      
    }
  )
  
  ###################################
  
  
  # Output graph ####################
  ## MB totally reworked to get sane graph of WNT, SHH, Grp3, Grp4
  output$classifierPlot <- renderPlot({
    
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    Total.No.of.Samples <- classified_data$Total.No.of.Samples 
    probs2 <- classified_data$probs2
    maxProbsWhich <- classified_data$maxProbsWhich
    maxProbs <- classified_data$maxProbs
    maxProbsCol <- classified_data$maxProbsCol
    maxProbsCol2 <- classified_data$maxProbsCol2
    percentfail <- classified_data$percentfail
    
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
    new.percentfail <- percentfail[index]
    
    # New remove Array QC fail samples
    cat(paste("Removing data points below percent_pval_cut_off", percent_pval_cutoff, "from graph:\n"))
    index <- !new.percentfail > percent_pval_cutoff
    cat(names(new.maxProbs[!index]), "\n")
    final.probs2 <- new.probs2[,index]
    final.maxProbs <- new.maxProbs[index]
    final.maxProbsWhich <- new.maxProbsWhich[index]
    final.maxProbsCol <- new.maxProbsCol[index]
    final.maxProbsCol2 <- new.maxProbsCol2[index]
    final.Total.No.of.Samples <- length(new.maxProbs[index])
    
    par(mfrow=c(1,1))
    #par(mar=c(6,4,2,1) + 0.1)
    par(mar=c(6,4,4,1) + 0.1)
    par(cex=1.3)
    par(cex.axis=1)
    
    heading <- paste("Medulloblastoma subgroup call confidence intervals for", final.Total.No.of.Samples, "samples")
    
    boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",final.probs2[,order(final.maxProbsWhich, final.maxProbs)],outpch=NA,ylim=c(0,1),las=2,
            col=final.maxProbsCol2[order(final.maxProbsWhich,final.maxProbs)] )
    
    abline(col="grey",lty = 1, h = threshold)
    
    # How many subgroups of each colour are we plotting
    tmp <- table(final.maxProbsCol)
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
    points(col=final.maxProbsCol[order(final.maxProbsWhich,final.maxProbs)],pch=19, final.maxProbs[order(final.maxProbsWhich,final.maxProbs)])
    legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4"), col=c("blue", "red", "yellow2", "darkgreen"), pch=19)   
    axis(2, las=2)     
    
  }) 
  # End output graph ################
  
  
  # Output graph download ####################
  ## MB totally reworked to get sane graph of WNT, SHH, Grp3, Grp4
  output$PlotDownload <- downloadHandler(
    
    filename = "MB_classification.png",
    content = function(file) {
      
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      Total.No.of.Samples <- classified_data$Total.No.of.Samples 
      probs2 <- classified_data$probs2
      maxProbsWhich <- classified_data$maxProbsWhich
      maxProbs <- classified_data$maxProbs
      maxProbsCol <- classified_data$maxProbsCol
      maxProbsCol2 <- classified_data$maxProbsCol2
      
      # Code to remove samples below threshold from plot
      cat(paste("Removing data points below threshold", threshold, "from graph:\n"))
      index <- maxProbs > threshold
      cat(names(maxProbs[!index]), "\n")
      new.probs2 <- probs2[,index]
      new.maxProbs <- maxProbs[index]
      new.maxProbsWhich <- maxProbsWhich[index]
      new.maxProbsCol <- maxProbsCol[index]
      new.maxProbsCol2 <- maxProbsCol2[index]
      new.Total.No.of.Samples <- length(maxProbs[index])
      
      # New remove Array QC fail samples
      cat(paste("Removing data points below percent_pval_cut_off", percent_pval_cutoff, "from graph:\n"))
      index <- !new.percentfail > percent_pval_cutoff
      cat(names(new.maxProbs[!index]), "\n")
      final.probs2 <- new.probs2[,index]
      final.maxProbs <- new.maxProbs[index]
      final.maxProbsWhich <- new.maxProbsWhich[index]
      final.maxProbsCol <- new.maxProbsCol[index]
      final.maxProbsCol2 <- new.maxProbsCol2[index]
      final.Total.No.of.Samples <- length(new.maxProbs[index])
      
      heading <- paste("Medulloblastoma subgroup call confidence intervals for", final.Total.No.of.Samples, "samples")
      
      png(file, height = 1280, width = 1440)
      par(mfrow=c(1,1))
      par(mar=c(7,4,4,1) + 0.1)
      par(cex=2)
      par(cex.axis=1)
      boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",final.probs2[,order(final.maxProbsWhich, final.maxProbs)],outpch=NA,ylim=c(0,1),las=2,
              col=final.maxProbsCol2[order(final.maxProbsWhich,final.maxProbs)] )
      
      abline(col="grey",lty = 1, h = threshold)
      
      # How many subgroups of each colour are we plotting
      tmp <- table(final.maxProbsCol)
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
      abline(v=grp.sum)
      #lines(col="black",lwd=2,new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
      points(col=final.maxProbsCol[order(final.maxProbsWhich,final.maxProbs)],pch=19, final.maxProbs[order(final.maxProbsWhich,final.maxProbs)])
      legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4"), col=c("blue", "red", "yellow2", "darkgreen"), pch=19)   
      axis(2, las=2)
      dev.off()
      
    })
  # End output graph download ################
  
  # Output density plot ####################
  output$DensityPlot <- renderPlot({
    
    # Get data
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    # Extract RGset
    RGset <- classified_data$RGset
    pd <- classified_data$pd
    
    # Plot
    densityPlot(RGset, main = "Beta value distribution")
    
  })
  # End output density plot

  # Download density plot ####################
  output$DenistyPlotDownload <- downloadHandler(
    
    filename = "Beta_value_distribution.png",
    content = function(file) {
      
      # Get data
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      # Extract RGset
      RGset <- classified_data$RGset
      pd <- classified_data$pd
      
      # Plot
      png(file, height = 1280, width = 1440)
      par(cex=2)
      par(cex.axis=1)
      densityPlot(RGset, main = "Beta value distribution")
      dev.off()
      
    })
  # End density plot download
  
  # Output density bean plot ####################
  output$DensityBeanPlot <- renderPlot({
    
    # Get data
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    # Extract RGset
    RGset <- classified_data$RGset
    pd <- classified_data$pd
    
    # Plot
    par(mar=c(5.1,9.1,4.1,2.1))
    densityBeanPlot(RGset, sampNames = pd$Sample_Name, main = "Beta value density bean plot")
    
  })
  # End output bean density plot
  
  # Download density bean plot ####################
  output$DenistyBeanPlotDownload <- downloadHandler(
    
    filename = "Beta_value_bean_plot.png",
    content = function(file) {
      
      # Get data
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      # Extract RGset
      RGset <- classified_data$RGset
      pd <- classified_data$pd
      
      # Plot
      par(mar=c(5.1,9.1,4.1,2.1))
      png(file, height = 1280, width = 1440)
      par(cex=2)
      par(cex.axis=1)
      densityBeanPlot(RGset, sampNames = pd$Sample_Name, main = "Beta value density bean plot")
      dev.off()
      
    })
  # End density bean plot download
  
  # Output time taken ###############
  
  output$time <- renderText({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    time <- classified_data$time
    
    c("Classification took", format(time[3]), "seconds")
  })
  
  ###################################
  

  
}) # End shinyServer
