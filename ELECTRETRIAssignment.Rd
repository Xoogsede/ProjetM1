\name{ELECTRETRIAssignment}
\alias{ELECTRETRIAssignment}

\title{
Electri-Tri sorting method by Amir NAFI and Caty WEREY
}
\description{
The electre TRI method uses the optimistic or pessimistic assignment rule, with indifference or preference or veto thresholds attached to criteria. It is the MRSort function modified.
}
\usage{
ELECTRETRIAssignment(performanceTable, 
      indefferenceThresholdVector,
      preferenceThresholdVector,
      vetosThresholdVector,
      credibilityThreshold,
      affectationMethod,
      categoriesLowerProfiles, 
      categoriesRanks, criteriaWeights, criteriaMinMax,
      alternativesIDs, criteriaIDs,
      categoriesIDs)
}

\arguments{
  \item{performanceTable}{Matrix or data frame containing the performance table. 
  Each row corresponds to an alternative, and each column to a criterion. 
  Rows (resp. columns) must be named according to the IDs of the alternatives (resp. criteria).}
  
  \item{indefferenceThresholdVector}{Vector of the indefference  thresholds for each criteria.}
  
  \item{preferenceThresholdVector}{Vector of the preference  thresholds for each criteria.}
  
  \item{vetoThresholdVector}{Vector of the veto  thresholds for each criteria.
  NA values mean that no threshold is defined}
  
  \item{credibilityThreshold}{The value of the cut threshold (sigma)}
  
  \item{affectationMethod}{Choise of the assignment method.
  "opt" (resp. "pess") indicate that the optimistic (resp. pessimistic) methode has to be used}
  
  \item{categoriesLowerProfiles}{Matrix containing, in each row, the lower profiles of the categories. 
  The columns are named according to the criteria, and the rows are named according to the categories. 
  The index of the row in the matrix corresponds to the rank of the category.}
  
  \item{categoriesRanks}{A vector containing the ranks of the categories (1 for the best,
  with higher values for increasingly less preferred categories). The vector needs to be 
  named with the categories names, whereas the ranks need to be a range of values from 1
  to the number of categories.}
  
  \item{criteriaWeights}{Vector containing the weights of the criteria. 
  The elements are named according to the IDs of the criteria.}
  
  \item{criteriaMinMax}{Vector containing the preference direction on each of the criteria. 
  "min" (resp. "max") indicates that the criterion has to be minimized (maximized). 
  The elements are named according to the IDs of the criteria.}
 
  \item{alternativesIDs}{Vector containing IDs of alternatives, according to which the datashould be filtered.}
  \item{criteriaIDs}{Vector containing IDs of criteria, according to which the data should be filtered.}
  \item{categoriesIDs}{Vector containing IDs of categories, according to which the data should be filtered.}
}



\value{
The function returns a vector containing the assignments of the alternatives to the categories.
}

\references{
MRSort function and "Aide à la décision multicritère : introduction aux méthodes
d'analyse multicritère de type ELECTRE" by Amir NAFI and Caty WEREY
}


\examples{
# the performance table

performanceTable <- rbind(
  c(2.63, 5.26, 52.63, 84.21, 26.32),
  c(0, 0, 492.31, 492.31, 0),
  c(0, 10.13, 20.25, 20.25, 405.06),
  c(0, 0, 0, 0, 0),
  c(0, 210.53, 280.7, 280.7, 1171.93))


rownames(performanceTable) <- c("T1","T2","T3","T4","T5")
# Cr1= "Bouchage du tronçon", Cr2= "Effondrement de la paroi du tronçon", Cr3="Ensablement du tronçon", Cr4="Détérioration de la capacité hydraulique", Cr5="Présence d'infiltrations"
colnames(performanceTable) <- c("Cr1", "Cr2", "Cr3", "Cr4", "Cr5")


# thresholds 

indefferenceThresholdVector <- c(0.5, 0.5, 0.5, 0.5, 0.5)

preferenceThresholdVector <- c(1, 1, 1, 1, 1)

# vetos
# vetosThresholdVector <- rbind(c(NA, NA, NA, NA, NA),c(NA, NA, NA, NA, NA),c(NA,NA,NA, NA, NA))

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

names(categoriesRanks) <- c("Good","Medium","Bad")

# criteria to minimize or maximize

criteriaMinMax <- c("min","min","min", "min", "min")

names(criteriaMinMax) <- colnames(performanceTable)


# weights

criteriaWeights <- c(0.1, 0.35, 0.1, 0.1, 0.35)

names(criteriaWeights) <- colnames(performanceTable)


# MRSort

assignments<-MRSort(performanceTable, categoriesLowerProfiles,
                    categoriesRanks,criteriaWeights,
                    criteriaMinMax, 3, 
                    vetoThresholdVector = vetoThresholdVector)
 
print(assignments)

# un peu de filtrage

assignments<-MRSort(performanceTable, categoriesLowerProfiles, 
                     categoriesRanks, criteriaWeights,
                    criteriaMinMax, 2,
                    categoriesIDs = c("Medium","Bad"), 
                    criteriaIDs = c("Price","Time"), 
                    alternativesIDs = c("RER", "BUS"))

print(assignments)

}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}

