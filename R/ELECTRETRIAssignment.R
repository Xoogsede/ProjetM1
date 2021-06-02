ELECTRETRIAssignment <- function(performanceTable, 
                                 indefferenceThresholdVector, 
                                 preferenceThresholdVector, 
                                 vetoThresholdVector=NULL,
                                 cutThreshold =NULL,
                                 affectationMethod = NULL,
                                 categoriesLowerProfiles,  
                                 categoriesRanks, 
                                 criteriaWeights, 
                                 criteriaMinMax,
                                 alternativesIDs = NULL, 
                                 criteriaIDs = NULL, 
                                 categoriesIDs = NULL)
{
  
    ## check the input data
   { 
      if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
        stop("wrong performanceTable, should be a matrix or a data frame")
      
      
      if (!(is.vector(indefferenceThresholdVector)))
          stop("indefferenceThresholdVector should be a vector")
          
      if (!(is.vector(preferenceThresholdVector)))
          stop("preferenceThresholdVector should be a vector")
          
      if (!((is.null(vetoThresholdVector) || (is.vector(vetoThresholdVector)))))
        stop("vetoThresholdVector should be a vector")
  
      if (!((is.matrix(categoriesLowerProfiles)||(is.data.frame(categoriesLowerProfiles)))))
        stop("categoriesLowerProfiles should be a matrix or a data frame")
      
      if (!(is.vector(categoriesRanks)))
        stop("categoriesRanks should be a vector")
      
      if(is.null(names(categoriesRanks)))
        stop("categoriesRanks should be named")
      
      if(!all(sort(categoriesRanks) == 1:length(categoriesRanks)))
        stop("categoriesRanks should contain a permutation of the category indices (from 1 to the number of categories)")
      
      if (!(is.vector(criteriaMinMax)))
        stop("criteriaMinMax should be a vector")
      
      if (!(is.vector(criteriaWeights)))
        stop("criteriaWeights should be a vector")
      
      if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
        stop("alternativesIDs should be a vector")
      
      if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
        stop("criteriaIDs should be a vector")
      
      if (!(is.null(categoriesIDs)) && !(is.vector(categoriesIDs)))
        stop("categoriesIDs should be a vector")
  }
  

  
    # declaration of the different constants to simplify the calculations.
   { 
      G <- performanceTable  
      
      B <- categoriesLowerProfiles 
      catR <- categoriesRanks
      
      K <- criteriaWeights
      Crt <- criteriaIDs 
      F <- 1:ncol(performanceTable)                # criteria indices  
      
      
      crtmm <- criteriaMinMax
      A <- alternativesIDs        # Alternatives
      catIDs <- categoriesIDs
      
      n <- length(F)    # number of criteria 
      m <- length(A)    # number of alternatives  
      h <- nrow(categoriesLowerProfiles)      # number of categori lower profil
      numCat <- length(categoriesRanks)
      
      q=indefferenceThresholdVector
      p=preferenceThresholdVector
      
      v=vetoThresholdVector
      if (is.null(v)){              # Creation of the null vector for to avoid calculation problem
        v <- rep(0, n)
        vetoThresholdVector <- v
      }
      
      lambda <- cutThreshold
      if (is.null(lambda)){
        lambda <-0.76
      }
      
      if (is.null(affectationMethod)){
        affectationMethod <- "pess"  # The default method will be the pessimistic one           
      }
     
    }
    
    
    
    # Construction of concordance index table.
   { 
      c_j <- matrix(0, nrow = nrow(categoriesLowerProfiles) *2, ncol = ncol(performanceTable))          # Matrix of the partial concordance indices 
      cj <- NULL
      concornames <- NULL                             # For  partial Concordance names
      
      while(is.null(cj)) {
      # Table of partial concordance indices 
      for (i in seq_along(alternativesIDs)) {      # Loop for each alternative
        for (j in seq_along(F)) {    # Loop for each criteria
          if (criteriaMinMax[j]=="min"){
            # For  maximisation
          for (l in seq_along(categoriesLowerProfiles[,1])) {  # Loop for each category profile
            
            # Calculation of the c_j(a_i,b_h)  
            if (!is.null(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                length(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) == 1 && 
                !is.na(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])){
              c_j[l,j] <- 0
            }else if(!is.null((categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j]) < performanceTable[i,j]) && 
                     length((categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j]) < performanceTable[i,j]) == 1 && 
                     !is.na((categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j]) < performanceTable[i,j]) && 
                     (categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j]) < performanceTable[i,j]){
              c_j[l,j] <- 1
            }else {
              c_j[l,j] <- (preferenceThresholdVector[j]+performanceTable[i,j]-categoriesLowerProfiles[l,j])/(preferenceThresholdVector[j]-indefferenceThresholdVector[j])
            }
            
            # Calculation of the c_j(b_h, a_i) 
            if (!is.null(performanceTable[i,j] >= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                length(performanceTable[i,j] >= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) == 1 && 
                !is.na(performanceTable[i,j] >= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                performanceTable[i,j] >= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])){
              c_j[l+2,j] <- 0
            }else if(!is.null(performanceTable[i,j] < (indefferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                     length(performanceTable[i,j] < (indefferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) == 1 && 
                     !is.na(performanceTable[i,j] < (indefferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                     performanceTable[i,j] < (indefferenceThresholdVector[j]+categoriesLowerProfiles[l,j])){
              c_j[l+2,j] <- 1
            }else {
              c_j[l+2,j] <- (preferenceThresholdVector[j]-performanceTable[i,j]+categoriesLowerProfiles[l,j])/(preferenceThresholdVector[j]-indefferenceThresholdVector[j])
            }
          } 
        }
        else if (criteriaMinMax[j]=="max"){
        # For minimisation
          for (l in 1:h) {  # Loop for each category profile
              
              # Calculation of the c_j(a_i,b_h)  
              if (!is.null(performanceTable[i,j] >= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) && 
                  length(performanceTable[i,j] >= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) == 1 && 
                  !is.na(performanceTable[i,j] >= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) && 
                  performanceTable[i,j] >= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])){
                c_j[l,j] <- 0
              }else if(!is.null((categoriesLowerProfiles[l,j]+indefferenceThresholdVector[j]) > performanceTable[i,j]) && 
                       length((categoriesLowerProfiles[l,j]+indefferenceThresholdVector[j]) > performanceTable[i,j]) == 1 && 
                       !is.na((categoriesLowerProfiles[l,j]+indefferenceThresholdVector[j]) > performanceTable[i,j]) && 
                       (categoriesLowerProfiles[l,j]+indefferenceThresholdVector[j]) > performanceTable[i,j]){
                c_j[l,j] <- 1
              }else {
                c_j[l,j] <- (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j]-performanceTable[i,j])/(preferenceThresholdVector[j]-indefferenceThresholdVector[j])
              }
              
              # Calculation of the c_j(b_h, a_i) 
              if (!is.null(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  length(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) == 1 && 
                  !is.na(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])){
                c_j[l+2,j] <- 0
              }else if(!is.null(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j])) && 
                       length(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j])) == 1 && 
                       !is.na(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j])) && 
                       performanceTable[i,j] > (categoriesLowerProfiles[l,j]-indefferenceThresholdVector[j])){
                c_j[l+2,j] <- 1
              }else {
                c_j[l+2,j] <- (preferenceThresholdVector[j]+performanceTable[i,j]-categoriesLowerProfiles[l,j])/(preferenceThresholdVector[j]-indefferenceThresholdVector[j])
              }
            } 
        }
          
        }
        # Track of the partial concordance indices c_j(a,b_h) and c_j(b_h, a)
        cj<- rbind(cj,c_j)        
        colnames(cj) <- colnames(performanceTable)
      }
      
        # partial Concordance names
        for (z in 2:length(categoriesRanks)-1) {
          for (i in seq_along(alternativesIDs)) {
            concornames <- cbind(concornames,paste("cj(a",i,",","b",z,")", sep=""), paste("cj(b",z,",","a",i,")", sep=""))
          } 
        }
        ab <- seq(1,nrow(cj), 2)
        ba <- seq(2,nrow(cj), 2)
        cj_ordered <- rbind(cj[ab,], cj[ba,]) #cj(a,b) ordered
        rownames(cj_ordered) <- t(concornames)
        
        }
    }
    
    
    
    # Calculation of the global concordance indices C(a,b_h) et C(b_h, a)
   { 
      cj_ordered<- data.frame(cj_ordered)
      cj_ordered$Cglob <- t((criteriaWeights%*%t(cj_ordered))/sum(criteriaWeights))       # Matrix of the global concordance indices
      cj_ordered
      }
      
    
    
    # Construction of discordance index table.
   { 
      
      d_j <- matrix(0, nrow = nrow(categoriesLowerProfiles)*2, ncol = ncol(performanceTable))          # Matrix of the partial discordance indices 
      dj <- NULL
      discornames <- NULL                             # For  partial discordance names
      
      # Table of partial discordance indices 
      for (i in seq_along(alternativesIDs)) {      # Loop for each alternative
        for (j in seq_along(F)) {    # Loop for each criteria
          if (criteriaMinMax[j]=="min"){
            # For  maximisation
            for (l in seq_along(categoriesLowerProfiles[,1])) {  # Loop for each category profile
              
              # Calculation of the d_j(a_i,b_h)  
              if (!is.null(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  length(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) == 1 && 
                  !is.na(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])){
                d_j[l,j] <- 0
              }else if(!is.null((categoriesLowerProfiles[l,j]-vetoThresholdVector[j]) >= performanceTable[i,j]) && 
                       length((categoriesLowerProfiles[l,j]-vetoThresholdVector[j]) >= performanceTable[i,j]) == 1 && 
                       !is.na((categoriesLowerProfiles[l,j]-vetoThresholdVector[j]) >= performanceTable[i,j]) && 
                       (categoriesLowerProfiles[l,j]-vetoThresholdVector[j]) >= performanceTable[i,j]){
                d_j[l,j] <- 1
              }else {
                d_j[l,j] <- (performanceTable[i,j]-categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])/(vetoThresholdVector[j]-preferenceThresholdVector[j])
              }
              
              # Calculation of the d_j(b_h, a_i) 
              if (!is.null(performanceTable[i,j] <= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                  length(performanceTable[i,j] <= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) == 1 && 
                  !is.na(performanceTable[i,j] <= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                  performanceTable[i,j] <= (preferenceThresholdVector[j]+categoriesLowerProfiles[l,j])){
                d_j[l+2,j] <- 0
              }else if(!is.null(performanceTable[i,j] > (vetoThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                       length(performanceTable[i,j] > (vetoThresholdVector[j]+categoriesLowerProfiles[l,j])) == 1 && 
                       !is.na(performanceTable[i,j] > (vetoThresholdVector[j]+categoriesLowerProfiles[l,j])) && 
                       performanceTable[i,j] > (vetoThresholdVector[j]+categoriesLowerProfiles[l,j])){
                d_j[l+2,j] <- 1
              }else {
                d_j[l+2,j] <- (performanceTable[i,j]-categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])/(vetoThresholdVector[j]-preferenceThresholdVector[j])
              }
            } 
          }
          else if (criteriaMinMax[j]=="max"){
            # For minimisation
            for (l in 1:h) {  # Loop for each category profile
              
              # Calculation of the d_j(a_i,b_h)  
              if (!is.null(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) && 
                  length(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) == 1 && 
                  !is.na(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])) && 
                  performanceTable[i,j] <= (categoriesLowerProfiles[l,j]+preferenceThresholdVector[j])){
                d_j[l,j] <- 0
              }else if(!is.null((categoriesLowerProfiles[l,j]+vetoThresholdVector[j]) < performanceTable[i,j]) && 
                       length((categoriesLowerProfiles[l,j]+vetoThresholdVector[j]) < performanceTable[i,j]) == 1 && 
                       !is.na((categoriesLowerProfiles[l,j]+vetoThresholdVector[j]) < performanceTable[i,j]) && 
                       (categoriesLowerProfiles[l,j]+vetoThresholdVector[j]) < performanceTable[i,j]){
                d_j[l,j] <- 1
              }else {
                d_j[l,j] <- (performanceTable[i,j]-categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])/(vetoThresholdVector[j]-preferenceThresholdVector[j])
              }
              
              # Calculation of the d_j(b_h, a_i) 
              if (!is.null(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  length(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) == 1 && 
                  !is.na(performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])) && 
                  performanceTable[i,j] > (categoriesLowerProfiles[l,j]-preferenceThresholdVector[j])){
                d_j[l+2,j] <- 0
              }else if(!is.null(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-vetoThresholdVector[j])) && 
                       length(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-vetoThresholdVector[j])) == 1 && 
                       !is.na(performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-vetoThresholdVector[j])) && 
                       performanceTable[i,j] <= (categoriesLowerProfiles[l,j]-vetoThresholdVector[j])){
                d_j[l+2,j] <- 1
              }else {
                d_j[l+2,j] <- (categoriesLowerProfiles[l,j]-performanceTable[i,j]-preferenceThresholdVector[j])/(vetoThresholdVector[j]-preferenceThresholdVector[j])
              }
            } 
          }
          
        }
        # Track of the partial discordance indices d_j(a,b_h) and d_j(b_h, a)
        dj<- rbind(dj,d_j)                         
        colnames(dj) <- colnames(performanceTable)
      }
      
      # partial discordance names
      for (z in 2:length(categoriesRanks)-1) {
        for (i in seq_along(alternativesIDs)) {
          discornames <- cbind(discornames,paste("dj(a",i,",","b",z,")", sep=""), paste("dj(b",z,",","a",i,")", sep=""))
        } 
      }
      dj_ordered <- rbind(dj[seq(1,nrow(dj), 2),],dj[seq(2,nrow(dj), 2),]) #dj(a,b) ordered
      rownames(dj_ordered) <- t(discornames)
    }
    
    
    
    # Calculation of the degree of credibility
   { 
      sgma <- NULL                                         # Declaration of the credibility index vector           
      credibnames <- NULL                                  # Declaration of the credibility indices name           
      for (i in seq_along(dj_ordered[,1])) {
      if (max(dj_ordered[i,])==0) {
        sgma <- rbind(sgma,  cj_ordered$Cglob[i])
      }else if (max(dj_ordered[i,])==1) {
        sgma <- rbind(sgma, 0)
      }else if (max(dj_ordered) < 1) {
        for (j in 1:ncol(dj_ordered)) {
          if(cj_ordered$Cglob[i] < dj_ordered[i,j]){
            sgma <- rbind(sgma,cj_ordered$Cglob[i]*prod((1-dj_ordered[i,j])/(1-cj_ordered$Cglob[i])))
        } 
        }
      }
      }
      # Credibility row names
      for (z in 2:length(categoriesRanks)-1) {
        for (i in seq_along(alternativesIDs)) {
          credibnames <- cbind(credibnames,paste("sgma(a",i,",","b",z,")", sep=""), paste("sgma(b",z,",","a",i,")", sep=""))
        } 
      }
      rownames(sgma) <- t(credibnames)
      colnames(sgma) <- "credibility index"
    }
    
    
    
    # dertermination of the preference relation between a_i and b_h
   { 
      Outranks <- NULL
      # I     : Indefference (aSb_h and b_hSa ???aIb_h)
      # > (<) : Preference (a>b => aSb_h ; a<b => b_hSa)
      # R     : incomparable (not aSb_h and not b_hSa)
      for (i in seq(1,nrow(sgma), 2)) {
        if (sgma[i] >= lambda && sgma[i+1] >= lambda) {            # If sigma(a, b_h) < lambda => 0 means aSb_h
          Outranks<-rbind("I", Outranks)
        } else if (sgma[i] >= lambda && sgma[i+1] < lambda) {
          Outranks<-rbind("<", Outranks)                # else sigma(a, b_h) b_hSa
        } else if (sgma[i] < lambda && sgma[i+1] >= lambda) {
          Outranks<-rbind(">", Outranks)                # else sigma(a, b_h) b_hSa
        }else if (sgma[i] < lambda && sgma[i+1] < lambda) {
          Outranks<-rbind("R", Outranks)                # else sigma(a, b_h) b_hSa
        }
      } 
      outrankTable <- NULL
      Outranks <- data.frame(Outranks)
      for (i in seq(1, nrow(Outranks), ncol(categoriesLowerProfiles))) {
        outrankTable <- rbind(t(Outranks$Outranks[i:(i+ncol(categoriesLowerProfiles)-1)]), outrankTable)
      }
      colnames(outrankTable) <- rownames(performanceTable)
      rownames(outrankTable) <- rownames(categoriesLowerProfiles)
      outrankTable
       } 
     
    
    
    # Assignment  
   { 
      # pessimistic assignment
      # a) compare alternatives ("a") successively to "b(i)" , for i=p,p-1, ..., 0,
      # b) let "b(h)" = the first profile such that "a" outranks "b(h).",
      # affect "a" to the category C(h+1).
        
      
      if (!(affectationMethod=="opt")) {
        Pessimistic <-  matrix (0,nrow(categoriesLowerProfiles)+1 ,nrow(performanceTable) )
        rownames(Pessimistic) <- names(categoriesRanks)
        colnames(Pessimistic) <- rownames(performanceTable)
      
        # Pessimistic assignment procedure:
        k<-1
        
        for (i in seq_along(alternativesIDs)) {
          while (k<nrow(categoriesLowerProfiles)) {
            if (outrankTable[k,i]=="<" &&  k==1)  {
              Pessimistic[k,i]=rownames(Pessimistic)[k]
              Pessimistic[-k,i]=""
              k=k+1
            }else if(k>=1 && k<nrow(categoriesLowerProfiles)){
              if ((outrankTable[k,i]==">" && outrankTable[k+1,i]==">") || 
                  (outrankTable[k,i]=="I" && outrankTable[k+1,i]==">") || 
                  (outrankTable[k,i]=="R" && outrankTable[k+1,i]==">") ||
                  (outrankTable[k,i]=="R" && outrankTable[k+1,i]=="R")){
                Pessimistic[k+2,i]=rownames(Pessimistic)[k+2]
                Pessimistic[-(k+2),i]=""
                k=k+1 
              }else if ((outrankTable[k-1,i]=="R" && outrankTable[k,i]=="R") ||
                        (outrankTable[k-1,i]=="R" && outrankTable[k,i]==">") ||
                        (outrankTable[k-1,i]==">" && outrankTable[k,i]==">") ||
                        (outrankTable[k-1,i]=="I" && outrankTable[k,i]==">")){
                Pessimistic[k+1,i]=rownames(Pessimistic)[k+1]
                Pessimistic[-(k+1),i]=""
                k=k+1
              }
            }
            k=k+1
            
          }
          k<- 1 
        }
        Pessimistic <- data.frame(Pessimistic)
        rownames(Pessimistic) <- NULL
    
      
        Pessimistic <- data.frame(Pessimistic)
        rownames(Pessimistic) <- NULL
        assingTable <- NULL
        for (i in 1:ncol(Pessimistic)) {
          for (l in 1:nrow(Pessimistic)) {
            if (!(Pessimistic[l,i]=="")) {
              assingTable <- cbind(assingTable, Pessimistic[l,i])
            }
          }
        }
        assingTable <- data.frame(assingTable)
        colnames(assingTable) <- colnames(Pessimistic)
        rownames(assingTable) <- "Pessimistic"
        Pessimistic <- assingTable
    
        
      }else if(affectationMethod=="opt"){
        Optimistic <-  matrix (0,nrow(categoriesLowerProfiles)+1 ,nrow(performanceTable) )
        rownames(Optimistic) <- names(categoriesRanks)
        colnames(Optimistic) <- rownames(performanceTable)
        
        # optimistic assignment,
        # a) compare alternatives ("a") successively to "b(i)" ,for  i=1, 2, ..., p+1,
        # b) let "b(h)" = the first profile such that "b(h)" outranks ("a"),
        # affect "a" to the category C(h). So we start in the lower row of the outrankTable
        
        # Optimistic assignment procedure: 
        k<-nrow(categoriesLowerProfiles)
        
          for (i in 1:nrow(performanceTable)) {
            while (k>1) {
              if (outrankTable[k,i]==">" &&  k==nrow(categoriesLowerProfiles))  {
                  Optimistic[k+1,i]=rownames(Optimistic)[k+1]
                  Optimistic[-(k+1),i]=""
                  k=k-1
                }else if(k>1 && k<=nrow(categoriesLowerProfiles)){
                  if ((outrankTable[k,i]=="<" && outrankTable[k-1,i]==">") || 
                      (outrankTable[k,i]=="I" && outrankTable[k-1,i]==">") || 
                      (outrankTable[k,i]=="R" && outrankTable[k-1,i]==">") ||
                      (outrankTable[k,i]=="I" && outrankTable[k-1,i]==">" )){
                      Optimistic[-k,i]=rownames(Optimistic)[k]
                      Optimistic[k,i]=""
                      k=k-1 
                  }else if ((outrankTable[k,i]=="R" && outrankTable[k-1,i]=="R") ||
                            (outrankTable[k,i]=="R" && outrankTable[k-1,i]=="<") ||
                            (outrankTable[k,i]=="<" && outrankTable[k-1,i]=="<") ||
                            (outrankTable[k,i]=="I" && outrankTable[k-1,i]=="<")){
                      Optimistic[k-1,i]=rownames(Optimistic)[k-1]
                      Optimistic[-(k-1),i]=""
                      k=k-1
                    }
                }
                  k=k-1
                
            }
         k<- nrow(categoriesLowerProfiles) 
          }
        
        Optimistic <- data.frame(Optimistic)
        rownames(Optimistic) <- NULL
        assingTable <- NULL
        for (i in seq_along(Optimistic[1,])) {
          for (l in seq_along(Optimistic[,1])) {
            if (!(Optimistic[l,i]=="")) {
              assingTable <- cbind(assingTable, Optimistic[l,i])
            }
          }
        }
        assingTable <- data.frame(assingTable)
        colnames(assingTable) <- colnames(Optimistic)
        rownames(assingTable) <- "Optimistic"
        Optimistic <- assingTable
        
      }   
        
        
        assignments <- NULL
        if (!(affectationMethod=="opt")) {
          assignments <- Pessimistic
        }else{
          assignments <-Optimistic
        }
    }
    
    
    return(assignments)             
} 
