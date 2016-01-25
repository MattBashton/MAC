### Functions

# Function to generate nice green black red heatmaps

odd <- function(x) x!=as.integer(x/2)*2
even <- function(x) x==as.integer(x/2)*2

# Function colorpanel

colorpanel <- function(n,low='green',mid='black',high='red')
   {
     if(even(n)) warning("n is even: colors panel will not be symmetric")

     # convert to rgb
     low <- col2rgb(low)
     mid <- col2rgb(mid)
     high <- col2rgb(high)

     # determine length of each component
     lower <- floor(n/2)
     upper <- n - lower

     red  <- c(
               seq(low[1,1], mid [1,1], length=lower),
               seq(mid[1,1], high[1,1], length=upper)
               )/255

     green <- c(
                seq(low[3,1], mid [3,1], length=lower),
                seq(mid[3,1], high[3,1], length=upper)
                )/255
     blue <- c(
              seq(low[2,1], mid [2,1], length=lower),
               seq(mid[2,1], high[2,1], length=upper)
               )/255


     rgb(red,blue,green)
   }

# Generate green-black-red colourscale

bluepink <- function(n) colorpanel(n, 'blue', 'white', 'red' )


DW.MP.Match.and.Select<-function(
   input1.matrix,
   input2.matrix)
{


# Match the genes of the first dataset on the second and select those rows from the second

# start of methodology

   #print(c("Running MP.Match.and.Select... on: ", input1.matrix, " ", input2.matrix))

# Read input datasets

   m1 <- (input1.matrix)
   gs.names1 <- rownames(m1)
   gs.descs1 <- rownames(m1)
   sample.names1 <- colnames(m1)

   m2 <- (input2.matrix)
   gs.names2 <- rownames(m2)
   gs.descs2 <- rownames(m2)
   sample.names2 <- colnames(m2)

# Match features to first dataset and create matching m2 dataset

   gs.names3 <- intersect(gs.names1, gs.names2)

   locations2 <- match(gs.names3, gs.names2, nomatch=0)
   gs.names2 <- gs.names2[locations2]
   gs.descs2 <- gs.descs2[locations2]
   m2 <- m2[locations2, ]

# Save dataset


   V<-input2.matrix[gs.names2,]
   #exprs(V) <- m2
   return(as.matrix(V))

}

 DW.MP.Factors.Project.C <-function(
   input.eset,
   factors.list,
   postprojnorm = TRUE,
   #output.file,
   method = "pseudo-inverse") {
# method: pseudo-inverse, nnls-solver

   library(MASS)

# start of methodology

  # print(c("Running MP.Factors.Project... on: "))

# Read input dataset


   m <- input.eset
   gs.names <- rownames(input.eset)
   gs.descs <- rownames(input.eset)
   sample.names <- colnames(input.eset)

# Read factors dataset

   #W <- as.data.frame(factors.list@fit@W)
   W <- factors.list
   W.row.names <- row.names(W)
   W.row.descs <- row.names(W)
   W.names <- names(W)

# Match features to first dataset and create matching m2 dataset

   overlap <- intersect(gs.names, W.row.names)

   #print(c("Size of Input dataset=", length(gs.names), " genes"))
   #print(c("Size of W matrix (rows)=", length(W.row.names), " genes"))
   #print(c("Size of overlap=", length(overlap), " genes"))

   locations.m <- match(overlap, gs.names, nomatch=0)
   m2 <- m[locations.m, ]

   locations.W <- match(overlap, W.row.names, nomatch=0)

   W2<-W[locations.W, ]

   W2<-as.matrix(W2)

   if (method == "pseudo-inverse") {

# Project input dataset using factors input

   H <- ginv(W2) %*% m2

# three potential ways to deal with negative values created in the approximated projection
#
# I:
#   max.H <- max(H)
#   min.H <- min(H)
#   H <- (H - min.H)/(max.H - min.H)
#
# II:
#   H <- ifelse(H < 0, 0, H)
#
# III:
#  n.col <- length(H[1,])
#  for (i in 1:n.col) {
#        max.H <- max(H[,i])
#        min.H <- min(H[,i])
#        H[,i] <- (H[,i] - min.H)/(max.H - min.H)
#  }

 # print(c("projecting using pseudo-inverse..."))

 } else if  (method == "nnls-solver") {  # using a non-negative least square solver

   H <- matrix(0, nrow=length(W2[1,]), ncol= length(m2[1,]))

   for (i in 1:length(m2[1,])) {
     H[, i] <- MP.nnls.fit(W2, m2[, i], wsqrt=1, eps=0, rank.tol=1e-07)
   }

 # print(c("projecting using NNLS solver..."))


  } else {
    stop("unknown method")
  }

# Normalize projected dataset to the unit hypersphere

  if (postprojnorm == TRUE) {
     n.col <- length(H[1,])
     for (i in 1:n.col) {
        S.2 <- sqrt(sum(H[,i]*H[,i]))
#        S.2 <- sum(H[,i])
        H[,i] <- H[,i]/S.2
     }
  }


# Save projected dataset

   V <- data.frame(H)
   names(V) <- sample.names
   row.names(V) <- W.names
   return(V)

}


# Build SVM and tree models

#   print("Building SVM model...")

   one.over <- function(x) { return(100/length(x)) }
#   class.number.list <- split(rep(1, length(class.list.train)) , class.list.train)
#   class.weights  <- sapply(class.number.list, one.over)
#   print(c("class.weights=", class.weights))

#   svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel, cost = cost, gamma = gamma, probability = T)

#   print("Computing train set predictions...")

#   train.pred <- predict(object = svm.model, newdata = t(m.train), decision.values = T, probability = T)

#   dec.vals.train <- attr(train.pred, "decision.values")
#   prob.train <- signif(attr(train.pred, "probabilities"), digits=2)
#  confidence.vector <- vector(length=n.train, mode="numeric")
#   bscore <- vector(length=n.train, mode = "numeric")
#   max.k <- length(prob.train[1,])
#   random.pred.conf <- ((max.k - 1)/max.k)^2 + (max.k - 1)*(1/max.k)^2
#   for (ii in 1:n.train) {
#      probs <- sort(prob.train[ii,], decreasing=T)
#      confidence.vector[ii] <-  1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2))/random.pred.conf
#      confidence.vector[ii] <- signif(confidence.vector[ii], digits=3)
#      if (class.list.train[ii] == as.character(train.pred[ii])) {
#         bscore[ii] <- signif((1 - probs[1])^2, digits=2)
#      } else {
#         bscore[ii] <- signif(probs[1]^2, digits=2)
#      }
#   }
#   confidence.call <- ifelse(confidence.vector >= high.conf.thres, " H ", " L ")
#   error.call <- ifelse(class.list.train == as.character(train.pred), "   ", " * ")
#   no.call <- ifelse(confidence.vector >= high.conf.thres, 0, 1)
#   real.error <- ifelse(((no.call == 0) & (error.call == " * ")), 1, 0)
#   correct.call <- ifelse(((no.call == 0) & (error.call == "   ")), 1, 0)

#   col.symbols.train <- paste(confidence.call, error.call)
#   class.names <- names(data.frame(prob.train))
#   Brier.train <- signif(mean(bscore), digits=2)

#   train.results <- data.frame(cbind(as.character(sample.names.train), class.list.train, as.character(train.pred), error.call, confidence.call, confidence.vector, no.call, real.error, correct.call, prob.train, bscore))
#   names(train.results)[1] <- "Train Sample Name"
#   names(train.results)[2] <- "Actual"
#   names(train.results)[3] <- "Predicted"
#   names(train.results)[4] <- "Error (*)"
#   names(train.results)[5] <- "Conf (H/L)"
#   names(train.results)[6] <- "Conf"
#   names(train.results)[7] <- "No Call"
#   names(train.results)[8] <- "Real Error"
#   names(train.results)[9] <- "Correct Call"

#   names(train.results)[10 + length(class.phen.train)] <- "Brier score"
#   print(train.results)
#  print(c("Brier score (Train) = ", Brier.train))

#   write("Training Results \n", file = prediction.results.file, append = F)
#   write.table(train.results, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")

#   write(c("\n\n Brier score (Train) = ", Brier.train), file = prediction.results.file, append = T)

#   no.call.list <- split(no.call, class.list.train)
#   real.error.list <- split(real.error, class.list.train)
#   correct.call.list <- split(correct.call, class.list.train)
#   count.class <- c(sapply(no.call.list, length), length(no.call))
#   no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
#   real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
#   correct.call.class <- c(sapply(correct.call.list, sum), sum(correct.call))
#   train.pred.high.conf <- ifelse(no.call == 0,  as.character(train.pred), "-- no call")
#   print(c("train.pred.high.conf =", train.pred.high.conf))

#      no.call.class.pct <- no.call.class/count.class
#      real.error.class.pct <- real.error.class/count.class
#      correct.call.class.pct <- correct.call.class/count.class
#   perf.table.train <- data.frame(cbind(c(names(no.call.list), "Total"), count.class, no.call.class, no.call.class.pct, real.error.class, real.error.class.pct, correct.call.class, correct.call.class.pct))
#  names(perf.table.train) <-  c("Class", "Count", "No Call", "No Call (%)", "Real Error", "Real Error (%)", "Correct Call", "Correct Call (%)")
#   write.table(perf.table.train, file = prediction.results.file, append = T, quote=F, row.names=F, sep = "\t")
#   print(perf.table.train)

#   conf.table.train <- table(class.list.train, train.pred.high.conf)
#   conf.table.train <- data.frame(cbind(row.names(conf.table.train), conf.table.train))
#   print(conf.table.train)
#   write("\n\n Confusion Matrix (Train) \n", file = prediction.results.file, append = T)
#   write.table(conf.table.train, file = prediction.results.file, append = T, quote=F, #row.names=F, sep = "\t")

#   print("Building SVM model completed. Predicting test data...")


