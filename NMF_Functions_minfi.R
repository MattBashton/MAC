### Functions

DW.MP.Match.and.Select<-function(
   input1.matrix,
   input2.matrix)
{


# Match the genes of the first dataset on the second and select those rows from the second

# start of methodology

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
   return(as.matrix(V))

}

 DW.MP.Factors.Project.C <-function(
   input.eset,
   factors.list,
   postprojnorm = TRUE,
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
 