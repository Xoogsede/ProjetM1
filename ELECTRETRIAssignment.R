ELECTRETRIAssignment <- function(performanceTable, 
                                 indefferenceThresholdVector, preferenceThresholdVector, vetoThresholdVector=NULL,
                                 cutThreshold =NULL,
                                 affectationMethod = NULL,
                                 categoriesLowerProfiles,  
                                 categoriesRanks, criteriaWeights, criteriaMinMax,
                                 alternativesIDs = NULL, criteriaIDs = NULL, categoriesIDs = NULL)
  
  
  
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
    
    if (!(is.null(categoriesIDs) || is.vector(categoriesIDs)))
      stop("categoriesIDs should be a vector")
    
    if (!((is.null(affectationMethod) && (affectationMethod == 'opt')) || (is.null(affectationMethod) && (affectationMethod == 'pess'))))
      stop("affectationMethod should be opt (resp. pess) for optimistic (resp. pessimistic) method")
            
  }
  

  # declaration of the different constants to simplify the calculations.
  {
    G=performanceTable  
    
    B <- categoriesLowerProfiles 
    catR <- categoriesRanks
    
    K <- criteriaWeights
    Crt <- criteriaIDs 
    F <- 1:length(Crt)                # criteria indices  
    
    # Declare the number of criteria, the number of alternatives and number of categories
    numCrit <- dim(performanceTable)[2]
    numAlt <- dim(performanceTable)[1]
    numCat <- length(categoriesRanks)
    
    crtmm <- criteriaMinMax
    A <- alternativesIDs        # Alternatives
    catIDs <- categoriesIDs
    
    n <- length(F)    # number of criteria 
    m <- length(A)    # number of alternatives  
    h <- nrow(B)      # number of categori lower profil
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
    
   finalAffectation <- NULL                 # This for the final affectation table
    
  }
  

  # Construction of concordance index table.
  {
    c_j <- matrix(0, nrow = h*2, ncol = n)          # Matrix of the partial concordance indices 
    cj <- NULL
    concornames <- NULL                             # For  partial Concordance names
    
    # Table of partial concordance indices 
    
    for (i in 1:m) {      # Loop for each alternative
      for (j in F) {    # Loop for each criteria
        for (l in 1:h) {  # Loop for each category profile
          
          # Calculation of the c_j(a_i,b_h)  
          if (B[l,j]-G[i,j]>=p[j]){
            c_j[l,j] <- 0
          }else if(B[l,j]-G[i,j]<=q[j]){
            c_j[l,j] <- 1
          }else {
            c_j[l,j] <- (p[j]+G[i,j]-B[h,j])/(p[j]-q[j])
          }
          
          # Calculation of the c_j(b_h, a_i) 
          if (G[i,j]-B[l,j]>=p[j]){
            c_j[l+2,j] <- 0
          }else if(G[i,j]-B[l,j]<=q[j]){
            c_j[l+2,j] <- 1
          }else {
            c_j[l+2,j] <- (p[j]-G[i,j]+B[h,j])/(p[j]-q[j])
          }
        } 
      }
      
      # Track of the partial concordance indices c_j(a,b_h) and c_j(b_h, a)
      cj<- rbind(cj,c_j)        
      colnames(cj) <- colnames(G)
    }
    
      # partial Concordance names
      for (z in 2:numCat-1) {
        for (i in 1:m) {
        concornames <- cbind(concornames,paste("cj(a",i,",","b",z,")", sep=""), paste("cj(b",z,",","a",i,")", sep=""))
        } 
      }
      cj_ordered <- rbind(cj[seq(1,nrow(cj), 2),],cj[seq(2,nrow(cj), 2),]) #cj(a,b) ordered
      rownames(cj_ordered) <- t(concornames)
    
  }
  

  # Calculation of the global concordance indices C(a,b_h) et C(b_h, a)
  {
    
    Cglob <- t((K%*%t(cj))/sum(K))       # Matrix of the global concordance indices
    Cglob_ordered <- t(t(c(Cglob[seq(1,nrow(Cglob), 2)],Cglob[seq(2,nrow(Cglob), 2)])))
    colnames(Cglob_ordered) <- "GlobalConcordance"
    rownames(Cglob_ordered) <-(t(concornames))
    
    # Creation of the global concordance table
    Cglob1 <- NULL
    Cglob2 <- NULL
    CabNames <- NULL
    
    for (i in seq(1,nrow(Cglob), nrow(Cglob_ordered)/length(A)/2)) {
      Cglob1 <- rbind(Cglob1,c(Cglob_ordered[i], Cglob_ordered[i+1]))
    }
    for (i in 1:length(A)) {
      Cglob2 <- rbind(Cglob2,c(Cglob1[i,], Cglob1[i+length(A),]))
    }
   
    for (i in 1:h) {
      CabNames <- cbind(CabNames,paste("cj(ai",",","b",i,")", sep=""), paste("cj(b",i,",","ai",")", sep=""))
    }  
    rownames(Cglob2) <- A
    colnames(Cglob2) <- CabNames
  }
  

  # Construction of discordance index table.
  {
    # v_j <- B-row(G)
    d_j <- matrix(0, nrow = h*2, ncol = n)          # Matrix of the partial discordance indices 
    dj <- NULL                                      # Table of partial discordance indices 
    
    for (i in 1:m) {      # Loop for each alternative
      for (j in F) {    # Loop for each criteria
        for (l in 1:h) {  # Loop for each category profile
          
          # Calculation of the d_j(a_i,b_h)  
          if (G[i,j] <= (B[l,j] + p[j])){
            d_j[l,j] <- 0
          }else if(G[i,j] > B[l,j] + v[j]){
            d_j[l,j] <- 1
          }else {
            d_j[l,j] <- (B[h,j]-G[i,j]-p[j])/(v[j]-p[j])
          }
          
          # Calculation of the d_j(b_h, a_i) 
          if (B[l,j] <= G[i,j] +p[j]){
            d_j[l+2,j] <- 0
          }else if( B[l,j] > G[i,j]+v[j]){
            d_j[l+2,j] <- 1
          }else {
            d_j[l+2,j] <- (G[i,j]-B[h,j]-p[j])/(v[j]-p[j])
          }
        } 
      }
      dj<- rbind(dj,d_j)                      # Track of the partial concordance indices d_j(a,b_h) and d_j(b_h, a)
    }
    dj_ordered <- rbind(dj[seq(1,nrow(dj), 2),],dj[seq(2,nrow(dj), 2),])
    
    
  }
  
  # Calculation of the degree of credibility
  {
    sgma <- NULL                                         # Declaration of the credibility index vector           
    for (i in 1:nrow(dj)) {
      tstVect <- (dj[i,] > Cglob[i] & dj[i,]==1)         # Booleen vector for testing if the particular cases of sigma
      if (TRUE %in% tstVect){                            # We test if the tstVect vector contain a TRUE value which means    
        sgma <- rbind(sgma, 0)                           # there is j such that dj(a, b) > C(a, b) and dj(a, b)=1 <=> sigma(a,b)=0
      }else if (Cglob[i]==1){                            # If C(a,b)=1 <=> sigma(a,b)=1 and dj(a,b)=0 for all j such that 
        sgma <- rbind(sgma, 1)                           # dj(a,b) > C(a,b)
      }else if (!(FALSE %in% (dj[i,] < Cglob[i]))){      # If for all j, dj(a, b) < C(a,b) then sigma(a,b) = C(a,b)
        sgma <- rbind(sgma, Cglob[i])
      }else{
        sgma <- rbind(sgma,Cglob[1]*prod((1-dj[i,])/(1-Cglob[i])))
      }
    }
    sgma_ordred <- t(t(c(sgma[seq(1,nrow(sgma), 2)], sgma[seq(2,nrow(sgma), 2)])))
    Tglob <- cbind(dj, Cglob, sgma)
    Tglob_ordered <- cbind(dj_ordered, Cglob_ordered, sgma_ordred)
  }
  
  # derternation of the preference relation between a_i and b_h
  {
  
    Surclass <- NULL
    for (i in seq(1,nrow(sgma_ordred), 2)) {
      if (sgma_ordred[i] >= lambda & sgma_ordred[i+1] >= lambda) {            # If sigma(a, b_h) < lambda => 0 means aSb_h
        Surclass<-rbind(0, Surclass)
      } else if (sgma_ordred[i] >= lambda & sgma_ordred[i+1] < lambda) {
        Surclass<-rbind(1, Surclass)                # else sigma(a, b_h) b_hSa
      } else if (sgma_ordred[i] < lambda & sgma_ordred[i+1] >= lambda) {
        Surclass<-rbind(-1, Surclass)                # else sigma(a, b_h) b_hSa
      }else if (sgma_ordred[i] < lambda & sgma_ordred[i+1] < lambda) {
        Surclass<-rbind("R", Surclass)                # else sigma(a, b_h) b_hSa
      }
    } 
  }    
    #################################### Le tableau des sigma à corriger###########################
    ###############################################################################################
    
    
  
  # fuctions 
  
  getCategory <- function(i)
    {
      if (affectationMethod == "pess") {
        
        for (k in (numCat-1):1)
        {
          cat <- names(categoriesRanks)[categoriesRanks == k]
          
          weightedSum <- 0
          
          for (crit in names(criteriaMinMax))
          {
            if (criteriaMinMax[crit] == "min")
            {
              if (performanceTable[i,crit] %<=% categoriesLowerProfiles[cat,crit])
                weightedSum <- weightedSum + criteriaWeights[crit]
            }
            else
            {
              if (performanceTable[i,crit] %>=% categoriesLowerProfiles[cat,crit])
                weightedSum <- weightedSum + criteriaWeights[crit]
            }
          }
          
          vetoActive <- FALSE
          
          if(!is.null(criteriaVetos))
          {
            for (crit in names(criteriaMinMax))
            {
              if(!is.na(criteriaVetos[cat,crit]) & !is.null(criteriaVetos[cat,crit]))
              {
                if (criteriaMinMax[crit] == "min")
                {
                  if (performanceTable[i,crit] %>=% criteriaVetos[cat,crit])
                  {
                    vetoActive <- TRUE
                    break
                  }
                }
                else
                {
                  if (performanceTable[i,crit] %<=% criteriaVetos[cat,crit])
                  {
                    vetoActive <- TRUE
                    break
                  }
                }
              }
            }
          }
          
          # stopping condition
          if(weightedSum < majorityThreshold || vetoActive)
            return(names(categoriesRanks)[categoriesRanks == (k + 1)])
        }
        # better than all profiles -> top categ
        return(names(categoriesRanks)[categoriesRanks == 1])
      }
    }
    
    
    
    
    
    
    
    if (affectationMethod == "opt") {
      for (i in seq(1,nrow(Surclass)/2,2)) {
        if(Surclass[i+1]=="-1"){
          finalAffectation <- rbind(Affect, cbind(alternativesIDs[i], categoriesIDs[numCat]))
          }else if (Surclass[i+1]=="0" & Surclass[i]=="-1") {
              finalAffectation <- rbind(Affect, cbind(alternativesIDs[i], categoriesIDs[numCat]))
          }
        }
      }
  
  
  # Example
  {
    # the performance table
    
    
    performanceTable <- rbind(
      c(2.63, 5.26, 52.63, 84.21, 26.32),
      c(0, 0, 492.31, 492.31, 0),
      c(0, 10.13, 20.25, 20.25, 405.06),
      c(0, 0, 0, 0, 0),
      c(0, 210.53, 280.7, 280.7, 1171.93))
    
    alternativesIDs <- c("T1","T2","T3","T4","T5")
    rownames(performanceTable) <- alternativesIDs
    # Cr1= "Bouchage du tronçon", Cr2= "Effondrement de la paroi du tronçon", Cr3="Ensablement du tronçon", Cr4="Détérioration de la capacité hydraulique", Cr5="Présence d'infiltrations"
    criteriaIDs <- c("Cr1", "Cr2", "Cr3", "Cr4", "Cr5")
    colnames(performanceTable) <- criteriaIDs
    
    
    # thresholds 
    
    indefferenceThresholdVector <- c(0.5, 0.5, 0.5, 0.5, 0.5)
    
    preferenceThresholdVector <- c(1, 1, 1, 1, 1)
    
    # vetos
    # vetoThresholdVector <- rbind(c(NA, NA, NA, NA, NA),c(NA, NA, NA, NA, NA),c(NA,NA,NA, NA, NA))
    
    colnames(vetoThresholdVector) <- colnames(performanceTable)
    rownames(vetoThresholdVector) <- c("Good","Medium","Bad")
    
    
    # lower profiles of the categories 
    # (best category in the first position of the list)
    
    categoriesLowerProfiles <- rbind(
      c(0.5, 2.23, 37.53, 37.96, 40.67),
      c(7, 29, 126.7, 128, 300))
    
    colnames(categoriesLowerProfiles) <- colnames(performanceTable)
    
    rownames(categoriesLowerProfiles)<-c("b1","b2")
    
    # the order of the categories, 1 being the best
    
    categoriesRanks <-c(1,2,3)
    categoriesIDs <- c("Good","Medium","Bad")
    names(categoriesRanks) <- categoriesIDs
    
    # criteria to minimize or maximize
    
    criteriaMinMax <- c("min","min","min", "min", "min")
    
    names(criteriaMinMax) <- colnames(performanceTable)
    
    
    # weights
    
    criteriaWeights <- c(0.1, 0.35, 0.1, 0.1, 0.35)
    
    names(criteriaWeights) <- colnames(performanceTable)
    
  }
  
