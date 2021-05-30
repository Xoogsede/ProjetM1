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
    F <- 1:ncol(G)                # criteria indices  
    
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
  
  # check if we have a lower profile for the worst category
  { 
    worstCat <- names(categoriesRanks)[categoriesRanks == length(categoriesRanks)]
    
    if(!(worstCat %in% rownames(categoriesLowerProfiles)))
    {
      categoriesLowerProfiles <- rbind(categoriesLowerProfiles, rep(NA,length(criteriaMinMax)))
      
      rownames(categoriesLowerProfiles)[length(categoriesRanks)] <- worstCat
    }
    
    if (!is.null(criteriaVetos))
    {
      if(!(worstCat %in% rownames(criteriaVetos)))
      {
        criteriaVetos <- rbind(criteriaVetos, rep(NA,length(criteriaMinMax)))
        
        rownames(criteriaVetos)[length(categoriesRanks)] <- worstCat
      }
    }
  }
  
  ## filter the data according to the given alternatives and criteria
  {
    if (!is.null(alternativesIDs)){
      performanceTable <- performanceTable[alternativesIDs,,drop=FALSE]
    } 
    
    if (!is.null(criteriaIDs)){
      performanceTable <- performanceTable[,criteriaIDs,drop=FALSE]
      criteriaWeights <- criteriaWeights[criteriaIDs,drop=FALSE]
      criteriaMinMax <- criteriaMinMax[criteriaIDs,drop=FALSE]
      categoriesLowerProfiles <- categoriesLowerProfiles[,criteriaIDs,drop=FALSE]
    }
    
    if ((!is.null(criteriaIDs)) && (!is.null(criteriaVetos))){
      criteriaVetos <- criteriaVetos[,criteriaIDs,drop=FALSE]  
    }
    
    if ((!is.null(categoriesIDs)) && (!is.null(criteriaVetos))){
      criteriaVetos <- criteriaVetos[categoriesIDs,,drop=FALSE]
    }
    
    if (!is.null(categoriesIDs)){
      categoriesLowerProfiles <- categoriesLowerProfiles[categoriesIDs,,drop=FALSE]
    }
    
    if (!is.null(categoriesIDs)){
      # filter out categories
      categoriesRanks <- categoriesRanks[names(categoriesRanks) %in% categoriesIDs]
      # check if we took out all categories
      if(length(categoriesRanks) == 0)
        stop('categoriesIDs have filtered out all categories')
      # order the remaining ones
      categoriesRanks <- sort(categoriesRanks)
      # store their order
      catOrder <- names(categoriesRanks)
      # adjust their indices to a range from 1 to the number of remaining categories
      categoriesRanks <- 1:length(categoriesRanks)
      # rename them
      names(categoriesRanks) <- catOrder
    }
  }

    # Construction of concordance index table.
  {
    c_j <- matrix(0, nrow = h*2, ncol = n)          # Matrix of the partial concordance indices 
    cj <- NULL
    concornames <- NULL                             # For  partial Concordance names
    
    # Table of partial concordance indices 
    for (i in 1:m) {      # Loop for each alternative
      for (j in F) {    # Loop for each criteria
        if (criteriaMinMax[j]=="max"){
          # For  maximisation
        for (l in 1:h) {  # Loop for each category profile
          
          # Calculation of the c_j(a_i,b_h)  
          if (G[i,j] <= (B[l,j]-p[j])){
            c_j[l,j] <- 0
          }else if((B[l,j]-q[j]) < G[i,j]){
            c_j[l,j] <- 1
          }else {
            c_j[l,j] <- (p[j]+G[i,j]-B[l,j])/(p[j]-q[j])
          }
          
          # Calculation of the c_j(b_h, a_i) 
          if (G[i,j] >= (p[j]+B[l,j])){
            c_j[l+2,j] <- 0
          }else if(G[i,j] < (q[j]+B[l,j])){
            c_j[l+2,j] <- 1
          }else {
            c_j[l+2,j] <- (p[j]-G[i,j]+B[l,j])/(p[j]-q[j])
          }
        } 
      }
      else if (criteriaMinMax[j]=="min"){
      # For minimisation
        for (l in 1:h) {  # Loop for each category profile
            
            # Calculation of the c_j(a_i,b_h)  
            if (G[i,j] >= (B[l,j]+p[j])){
              c_j[l,j] <- 0
            }else if((B[l,j]+q[j]) > G[i,j]){
              c_j[l,j] <- 1
            }else {
              c_j[l,j] <- (p[j]+B[l,j]-G[i,j])/(p[j]-q[j])
            }
            
            # Calculation of the c_j(b_h, a_i) 
            if (G[i,j] <= (B[l,j]-p[j])){
              c_j[l+2,j] <- 0
            }else if(G[i,j] > (B[l,j]-q[j])){
              c_j[l+2,j] <- 1
            }else {
              c_j[l+2,j] <- (p[j]+G[i,j]-B[l,j])/(p[j]-q[j])
            }
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
  cj_ordered<- data.frame(cj_ordered)
  cj_ordered$Cglob <- t((K%*%t(cj_ordered))/sum(K))       # Matrix of the global concordance indices
  cj_ordered
  }
  

  # Construction of discordance index table.
  {
  
  d_j <- matrix(0, nrow = h*2, ncol = n)          # Matrix of the partial discordance indices 
  dj <- NULL
  discornames <- NULL                             # For  partial discordance names
  
  # Table of partial discordance indices 
  for (i in 1:m) {      # Loop for each alternative
    for (j in F) {    # Loop for each criteria
      if (criteriaMinMax[j]=="max"){
        # For  maximisation
        for (l in 1:h) {  # Loop for each category profile
          
          # Calculation of the d_j(a_i,b_h)  
          if (G[i,j] > (B[l,j]-p[j])){
            d_j[l,j] <- 0
          }else if((B[l,j]-v[j]) >= G[i,j]){
            d_j[l,j] <- 1
          }else {
            d_j[l,j] <- (G[i,j]-B[l,j]-p[j])/(v[j]-p[j])
          }
          
          # Calculation of the d_j(b_h, a_i) 
          if (G[i,j] <= (p[j]+B[l,j])){
            d_j[l+2,j] <- 0
          }else if(G[i,j] > (v[j]+B[l,j])){
            d_j[l+2,j] <- 1
          }else {
            d_j[l+2,j] <- (G[i,j]-B[l,j]-p[j])/(v[j]-p[j])
          }
        } 
      }
      else if (criteriaMinMax[j]=="min"){
        # For minimisation
        for (l in 1:h) {  # Loop for each category profile
          
          # Calculation of the d_j(a_i,b_h)  
          if (G[i,j] <= (B[l,j]+p[j])){
            d_j[l,j] <- 0
          }else if((B[l,j]+v[j]) < G[i,j]){
            d_j[l,j] <- 1
          }else {
            d_j[l,j] <- (G[i,j]-B[l,j]-p[j])/(v[j]-p[j])
          }
          
          # Calculation of the d_j(b_h, a_i) 
          if (G[i,j] > (B[l,j]-p[j])){
            d_j[l+2,j] <- 0
          }else if(G[i,j] <= (B[l,j]-v[j])){
            d_j[l+2,j] <- 1
          }else {
            d_j[l+2,j] <- (B[l,j]-G[i,j]-p[j])/(v[j]-p[j])
          }
        } 
      }
      
    }
    # Track of the partial discordance indices d_j(a,b_h) and d_j(b_h, a)
    dj<- rbind(dj,d_j)                         
    colnames(dj) <- colnames(G)
  }
  
  # partial discordance names
  for (z in 2:numCat-1) {
    for (i in 1:m) {
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
  for (i in 1:nrow(dj_ordered)) {
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
  for (z in 2:numCat-1) {
    for (i in 1:m) {
      credibnames <- cbind(credibnames,paste("sgma(a",i,",","b",z,")", sep=""), paste("sgma(b",z,",","a",i,")", sep=""))
    } 
  }
  rownames(sgma) <- t(credibnames)
  colnames(sgma) <- "credibility index"
}

  
  # derternation of the preference relation between a_i and b_h
  {
  
    Outranks <- NULL
    # I     : Indefference (aSb_h and b_hSa ???aIb_h)
    # > (<) : Preference (a>b => aSb_h ; a<b => b_hSa)
    # R     : incomparable (not aSb_h and not b_hSa)
    for (i in seq(1,nrow(sgma), 2)) {
      if (sgma[i] >= lambda & sgma[i+1] >= lambda) {            # If sigma(a, b_h) < lambda => 0 means aSb_h
        Outranks<-rbind("I", Outranks)
      } else if (sgma[i] >= lambda & sgma[i+1] < lambda) {
        Outranks<-rbind("<", Outranks)                # else sigma(a, b_h) b_hSa
      } else if (sgma[i] < lambda & sgma[i+1] >= lambda) {
        Outranks<-rbind(">", Outranks)                # else sigma(a, b_h) b_hSa
      }else if (sgma[i] < lambda & sgma[i+1] < lambda) {
        Outranks<-rbind("R", Outranks)                # else sigma(a, b_h) b_hSa
      }
    } 
    outrankTable <- NULL
    Outranks <- data.frame(Outranks)
    for (i in seq(1, nrow(Outranks), ncol(B))) {
      outrankTable <- rbind(outrankTable, t(Outranks$Outranks[i:(i+ncol(B)-1)]))
    }
    colnames(outrankTable) <- rownames(G)
    rownames(outrankTable) <- rownames(B)
    outrankTable
   }    
    #################################### Le tableau des sigma à corriger###########################
    ###############################################################################################
    
    
  
  # fuctions 
{
  # pessimistic assignment
  # a) compare alternatives ("a") successively to "b(i)" , for i=p,p-1, ..., 0,
  # b) let "b(h)" = the first profile such that "a" outranks "b(h).",
  # affect "a" to the category C(h+1).
  
  
  if (!(affectationMethod=="opt")) {
    Pessimistic <-  matrix (0,nrow(B)+1 ,nrow(G) )
    rownames(Pessimistic) <- names(categoriesRanks)
    colnames(Pessimistic) <- rownames(G)
  
    # Pessimistic assignment procedure:
    for (i in 1:nrow(G)){
      for (j in 1:nrow(B)){
        if (outrankTable[j,i]=="<" && j==1 && !(max(Pessimistic[1:nrow(B)+1,i])==1)) {
          Pessimistic[j,i]=1
        }else{
          for (j in 1:nrow(B)) {
            if (!(outrankTable[j,i]=="<") && !(j==1) && !(max(Pessimistic[1:nrow(B)+1,i])==1)){
              Pessimistic[j,i]=1
            }else if ((outrankTable[j,i]==">" || outrankTable[j,i]=="R" || outrankTable[j,i]=="I") && j==nrow(B) && !(max(Pessimistic[1:nrow(B)+1,i])==1)) {
              Pessimistic[j+1,i]=1
            }else{
              Pessimistic[j+1,i]=0
            }
          }
        }
      }
    }
  
  View(Pessimistic)
  
  }else if(affectationMethod=="opt"){
    Optimistic <-  matrix (0,nrow(B)+1 ,nrow(G) )
    rownames(Optimistic) <- names(categoriesRanks)
    colnames(Optimistic) <- rownames(G)
    
    # optimistic assignment,
    # a) compare alternatives ("a") successively to "b(i)" ,for  i=1, 2, ..., p+1,
    # b) let "b(h)" = the first profile such that "b(h)" outranks ("a"),
    # affect "a" to the category C(h).
    
    # Optimistic assignment procedure: 
    for (i in 1:nrow(G)){
      for (j in nrow(B):1){
        if (outrankTable[j,i]==">" & j==nrow(B)) {
          Optimistic[j+1,i]=1
        }else{
          for (j in nrow(B):1) {
            if (!(outrankTable[j,i]==">") & !(j==1)){
              next
            }else if ((outrankTable[j,i]=="<" | outrankTable[j,i]=="R" | outrankTable[j,i]=="I") & j==1 & !(max(Optimistic[1:nrow(B)+1,i])==1)) {
              Optimistic[j,i]=1
            }else{
              Optimistic[j,i]=0
            }
          }
        }
      }
    }
    View(Optimistic)
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
  
