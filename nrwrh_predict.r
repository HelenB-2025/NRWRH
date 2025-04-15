
#########build network-based network########
#by incorporating new similarity matrice in terms of shared common targets/drugs
#into drug and target network matrix respectively.
#require drug chemical similarity data as AD1, target sequence similarity data 
#as AT1, and drug-target association data as B.
combinedsimilarity <- function(B, AD1, AT1, gamad1,gamat1,m, n) {
  # Initialize AD2 (drug-drug similarity based on common targets
  AD2 <- diag(rowSums(t(B)))  # t(B) because drugs are columns
  
  for (i in 2:m) {
    for (j in 1:(i - 1)) {
      k <- sum(B[, i] == 1 & B[, j] == 1 & B[, i] == B[, j])
      AD2[i, j] <- k
    }
  }
  
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      AD2[i, j] <- AD2[j, i]
    }
  }
  
  # Initialize AT2 (target-target similarity based on common drugs)
  AT2 <- diag(rowSums(B))
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      k <- sum(B[i, ] == 1 & B[j, ] == 1 & B[i, ] == B[j, ])
      AT2[i, j] <- k
    }
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      AT2[i, j] <- AT2[j, i]
    }
  }
  
  diag(AT2)[diag(AT2) == 0] <- 1
  diag(AD2)[diag(AD2) == 0] <- 1
  
  # Normalize AT2 and AD2
  NAT2 <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      NAT2[i, j] <- AT2[i, j] / sqrt(AT2[i, i]) / sqrt(AT2[j, j])
    }
  }
  
  NAD2 <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      NAD2[i, j] <- AD2[i, j] / sqrt(AD2[i, i]) / sqrt(AD2[j, j])
    }
  }

  # Combine similarities
  AD <- (gamad1 *AD1 + (1-gamad1) * NAD2) 
  AT <- (gamat1 * AT1 + (1-gamat1) * NAT2)
  
  return(list(AD = AD, AT = AT))
}



#########For prediction#########
#it calls the previously built function to generate network-based network
nrwrhdrugtarget <- function(backprobability,lambda,eta,gamad1,gamat1,drugID) {
  # Load data
  B<-as.matrix(read.table("code/example/interaction.txt"))
  AT1 <- as.matrix(read.table("code/example/sequencesimilarity.txt"))
  AD1 <- as.matrix(read.table("code/example/chemicalsimilarity.txt"))
drug=as.data.frame(AD1)
target=as.data.frame(AT1)
adj=as.data.frame(B)
  #HB: modify B by removing first column which is drug ID
  n <- nrow(B)  # number of targets
  m <- ncol(B)  # number of drugs
  
  #call previously built function to standardize drug and target own similarity matrix
  sim_list <- combinedsimilarity(B, AD1, AT1, gamad1,gamat1, m, n)
  AD <- sim_list$AD
  AT <- sim_list$AT
  
  # Construct transition matrices
  MTD <- matrix(0, n, m)
  MDT <- matrix(0, m, n)
  MT  <- matrix(0, n, n)
  MD  <- matrix(0, m, m)
  
  for (i in 1:n) {
    row_sum <- sum(B[i, ])
    if (row_sum != 0) {
      MTD[i, ] <- lambda * B[i, ] / row_sum
    }
  }
  
  for (i in 1:m) {
    col_sum <- sum(B[, i])
    if (col_sum != 0) {
      MDT[i, ] <- lambda * B[, i] / col_sum
    }
  }
  
  for (i in 1:n) {
    row_sum <- sum(AT[i, ])
    if (sum(B[i, ]) == 0 && row_sum != 0) {
      MT[i, ] <- AT[i, ] / row_sum
    } else if (row_sum != 0) {
      MT[i, ] <- (1 - lambda) * AT[i, ] / row_sum
    }
  }
  
  for (i in 1:m) {
    row_sum <- sum(AD[i, ])
    if (sum(B[, i]) == 0 && row_sum != 0) {
      MD[i, ] <- AD[i, ] / row_sum
    } else if (row_sum != 0) {
      MD[i, ] <- (1 - lambda) * AD[i, ] / row_sum
    }
  }
  
  M1 <- cbind(MT, MTD)
  M2 <- cbind(MDT, MD)
  M  <- rbind(M1, M2)

  # Initialize probabilities
  u0 <- rep(1/length(u0),n)
  
  v0 <- rep(0, m)
  v0[which(colnames(AD1)==drugID)] <- 1

  p0 <- c((1 - eta) * u0, eta * v0)
  p1 <- p0
  p <- (1 - backprobability) * t(M) %*% p1 + backprobability * p0
  
  while (max(abs(p - p1)) > 1e-10) {
    p1 <- p
    p <- (1 - backprobability) * t(M) %*% p1 + backprobability * p0
  }
  
  u <- p[1:n]

  # Rank targets
  ss <- as.data.frame(u)
  dline<-data.frame(known=adj[,which(colnames(adj)==drugID)])
  concate<-cbind(ss,dline,ind=seq(1:n))
  comp <- concate[order(-concate$u), ]
  
  #first order among unknown associated target
  first<-comp[which(comp$known==0)[1],]
  top=row.names(adj[first$ind,])  
  
  return(top)
}
