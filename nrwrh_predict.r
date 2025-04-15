
#########For leave-one-out cross validation#########
#compare fold enrichment with RWRH
valid <- function(backprobability,lambda,eta,gamad1,gamat1) {
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
  
  #extract all drug ID
  drug<-colnames(drug)
  dn<-length(drug)
  
  ###Leave-one-out cross-validation by masking one known association
  #for each drug
  for(l in 1:dn){
    drugID<-drug[l]
    
    #for each drug's known associated target
    if(sum(adj[,which(colnames(adj)==drugID)])>0){
      dv<-adj[drugID]
      ass_tar <- c(rownames(dv)[dv > 0])
      kt<-length(ass_tar)
      
      for(k in 1:kt){
        adj_mask<-adj
        adj_mask[ass_tar[k],drugID]<-0
        B_mask<-as.matrix(adj_mask)
        
        #call previously built function to standardize drug and target own similarity matrix
        sim_list <- combinedsimilarity(B_mask,AD1,AT1,gamad1,gamat1,m,n)
        AD <- sim_list$AD
        AT <- sim_list$AT
        
        # Construct transition matrices
        MTD <- matrix(0, n, m)
        MDT <- matrix(0, m, n)
        MT  <- matrix(0, n, n)
        MD  <- matrix(0, m, m)
        
        for (i in 1:n) {
          row_sum <- sum(B_mask[i, ])
          if (row_sum != 0) {
            MTD[i, ] <- lambda * B_mask[i, ] / row_sum
          }
        }
        
        for (i in 1:m) {
          col_sum <- sum(B_mask[, i])
          if (col_sum != 0) {
            MDT[i, ] <- lambda * B_mask[, i] / col_sum
          }
        }
        
        for (i in 1:n) {
          row_sum <- sum(AT[i, ])
          if (sum(B_mask[i, ]) == 0 && row_sum != 0) {
            MT[i, ] <- AT[i, ] / row_sum
          } else if (row_sum != 0) {
            MT[i, ] <- (1 - lambda) * AT[i, ] / row_sum
          }
        }
        
        for (i in 1:m) {
          row_sum <- sum(AD[i, ])
          if (sum(B_mask[, i]) == 0 && row_sum != 0) {
            MD[i, ] <- AD[i, ] / row_sum
          } else if (row_sum != 0) {
            MD[i, ] <- (1 - lambda) * AD[i, ] / row_sum
          }
        }
        
        M1 <- cbind(MT, MTD)
        M2 <- cbind(MDT, MD)
        M  <- rbind(M1, M2)
        
        
        
        #matrix used by RWRH without integrating shared drug/targets 
        MTD_rwrh <- matrix(0, n, m)
        MDT_rwrh <- matrix(0, m, n)
        MT_rwrh  <- matrix(0, n, n)
        MD_rwrh  <- matrix(0, m, m)
        for (i in 1:n) {
          row_sum <- sum(AT1[i, ])
          if (sum(B_mask[i, ]) == 0 && row_sum != 0) {
            MT_rwrh[i, ] <- AT1[i, ] / row_sum
          } else if (row_sum != 0) {
            MT_rwrh[i, ] <- (1 - lambda) * AT1[i, ] / row_sum
          }
        }
        for (i in 1:m) {
          row_sum <- sum(AD1[i, ])
          if (sum(B_mask[, i]) == 0 && row_sum != 0) {
            MD_rwrh[i, ] <- AD1[i, ] / row_sum
          } else if (row_sum != 0) {
            MD_rwrh[i, ] <- (1 - lambda) * AD1[i, ] / row_sum
          }
        }
        M1_rwrh <- cbind(MT_rwrh, MTD)
        M2_rwrh <- cbind(MDT, MD_rwrh)
        M_rwrh  <- rbind(M1_rwrh, M2_rwrh)
        
        
        # Initialize probabilities
        #for NRWRH
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
        
        #for RWRH
        p1_rwrh <- p0
        p_rwrh <- (1 - backprobability) * t(M_rwrh) %*% p1_rwrh + backprobability * p0
        while (max(abs(p_rwrh - p1_rwrh)) > 1e-10) {
          p1_rwrh <- p_rwrh
          p_rwrh <- (1 - backprobability) * t(M_rwrh) %*% p1_rwrh + backprobability * p0
        }
        u_rwrh <- p_rwrh[1:n]
        
      
        # Rank targets
        ss <- as.data.frame(u)
        tname<-data.frame(tname=rownames(adj_mask))
        concate<-cbind(ss,tname)
        comp <- concate[order(-concate$u), ]
        #the order of masked target among unknown associated target
        ord<-which(comp["tname"] == ass_tar[k])
        enrich<-n/2/ord
        
        #rwrh
        ss_rwrh <- as.data.frame(u_rwrh)
        concate_rwrh<-cbind(ss_rwrh,tname)
        comp_rwrh <- concate_rwrh[order(-concate_rwrh$u_rwrh), ]
        #the order of masked target among unknown associated target
        ord_rwrh<-which(comp_rwrh["tname"] == ass_tar[k])
        enrich_rwrh<-n/2/ord_rwrh
        
        if(l==1&k==1){
          #create matrix for saving each drug's each target fold enrichment
          out_rwrh<-enrich_rwrh;out<-enrich
        }else{
          out_rwrh<-c(out_rwrh,enrich_rwrh);out<-c(out,enrich)
        }
      }#target loop end
    }#if have associated target condition ened
  }#drug loop end
  all<-list(rwrh=out_rwrh,nrwrh=out)
  return(all)
}
