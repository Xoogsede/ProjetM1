
# Example

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
  # vetosThresholdVector <- rbind(c(NA, NA, NA, NA, NA),c(NA, NA, NA, NA, NA),c(NA,NA,NA, NA, NA))
  #vetoThresholdVector <- NULL
  
 # "if (is.null(vetoThresholdVector)){              # Creation of the null vector for to avoid calculation problem
  #  vetoThresholdVector <- rep(0, ncol(performanceTable))
 # }"
  
  
  
  # lower profiles of the categories 
  # (best category in the first position of the list)
  
  categoriesLowerProfiles <- rbind(
    c(0.5, 2.23, 37.53, 37.96, 40.67),
    c(7, 29, 126.7, 128, 300))
  
  colnames(categoriesLowerProfiles) <- colnames(performanceTable)
  
  rownames(categoriesLowerProfiles)<-c("b1","b2")
  
  # the order of the categories, 1 being the best
  
  categoriesRanks <-c(1,2,3)
  categoriesIDs <- as.vector(c("Good","Medium","Bad"))
  names(categoriesRanks) <- categoriesIDs
  
  # criteria to minimize or maximize
  
  criteriaMinMax <- c("min","min","min", "min", "min")
  
  names(criteriaMinMax) <- colnames(performanceTable)
  
  
  # weights
  
  criteriaWeights <- c(0.1, 0.35, 0.1, 0.1, 0.35)
  
  names(criteriaWeights) <- colnames(performanceTable)
  

  
  #cutThreshold =NULL
  #affectationMethod = NULL

 

  

  assignmentsOPT <- ELECTRETRIAssignment(performanceTable = performanceTable, 
                                      indefferenceThresholdVector = indefferenceThresholdVector, 
                                      preferenceThresholdVector = preferenceThresholdVector,
                                      categoriesLowerProfiles = categoriesLowerProfiles, affectationMethod = "opt",  
                                      categoriesRanks = categoriesRanks, 
                                      criteriaWeights = criteriaWeights, 
                                      criteriaMinMax = criteriaMinMax)
  
  assignmentsPESS <- ELECTRETRIAssignment(performanceTable = performanceTable, 
                                         indefferenceThresholdVector = indefferenceThresholdVector, 
                                         preferenceThresholdVector = preferenceThresholdVector,
                                         categoriesLowerProfiles = categoriesLowerProfiles, affectationMethod = "pess",  
                                         categoriesRanks = categoriesRanks, 
                                         criteriaWeights = criteriaWeights, 
                                         criteriaMinMax = criteriaMinMax)
  


