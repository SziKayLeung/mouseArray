##---------------------------------------------------------------------#
##
## Title: Run EWAS
##
## Purpose of script: Run linear regression within each cell type/bulk
##                    
##                    This scripted is adapted from one co-authored
##                    by EJH and EMW for the MRC schizophrenia project
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# to run: sbatch..... 

# within cell type models

# pathology: lm(DNAm ~ QCmetrics$Group + QCmetrics$Pathology + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch


# group*age: DNAm ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch


# include null models of each with anova to test that pathology/interaction term is improving the model fit



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(lme4)
library(lmerTest)
library(dplyr)


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

# if cell type = neun+ source functions without random terms
# if cell type = neun- source functions w/ random terms

runEWAS<-function(row,QCmetrics){
  
  nullLM<-lm(row ~ QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch)
  
  interactionLM<-lm(row ~ QCmetrics$Group + QCmetrics$Age*QCmetrics$Group + QCmetrics$Age + QCmetrics$Sex + QCmetrics$Batch)
  
  
  # extract case control main effect and cell proportion effect
  return(c(summary(interactionLM)$coefficients["QCmetrics$GroupWT",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$Age",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
           summary(interactionLM)$coefficients["QCmetrics$GroupWT:QCmetrics$Age",c(1,2,4)],
           
           
           # extract cell specific case control effect
           summary(nullLM)$coefficients["QCmetrics$GroupWT",c(1,2,4)],
           summary(nullLM)$coefficients["QCmetrics$Age",c(1,2,4)],
           summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
           
           # run anova
           anova(interactionLM, nullLM)[2,"Pr(>F)"]))
  
}



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
projDir <- args[1]
#projDir <- "/lustre/projects/Research_Project-191406/cellSortedEWAS"
cellType <- args[2]
#cellType <- "NEUNpos"


normData<-file.path(projDir, "2_normalised/normalisedData.rdat")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(projDir)
load(normData)
# output = QCmetrics & normBeta if bulk
# OR QCmetrics & celltypeNormbeta


# Check if 'Group' column exists
if ("Group" %in% colnames(QCmetrics)) {
  # If 'Group' exists, convert it to a factor
  QCmetrics$Group <- as.factor(QCmetrics$Group)
} else if ("Phenotype" %in% colnames(QCmetrics)) {
  # If 'Group' doesn't exist but 'Phenotype' exists, create 'Group' from 'Phenotype'
  QCmetrics$Group <- as.factor(QCmetrics$Phenotype)
  print("Using 'Phenotype' as 'Group'")
} else if ("Genotype" %in% colnames(QCmetrics)) {
  # If 'Group' and 'Phenotype' don't exist, but 'Genotype' exists, create 'Group' from 'Genotype'
  QCmetrics$Group <- as.factor(QCmetrics$Genotype)
  print("Using 'Genotype' as 'Group'")
} else {
  # If neither 'Group', 'Phenotype', nor 'Genotype' exist
  print("No valid group column found")
}

# subset on cell type if 2nd argument provided
if(!is.na(cellType) && !is.null(cellType)){
  print(paste0("running EWAS on ", cellType, " cell type..."))
  ## subset to cell type samples
  QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]
  
  # subset beta matrix to cell type specific samples
  celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
}else{
  print("running EWAS on bulk")
  celltypeNormbeta<-normBeta
  cellType <- "bulk"
} 

# use covariate chip position if bulk tissue 
if(cellType == "bulk"){
  QCmetrics$Batch <- as.factor(QCmetrics$Chip_Position)
}else{
  QCmetrics$Batch <- as.factor(QCmetrics$Batch)
}

QCmetrics$Sex <- as.factor(QCmetrics$Sex)


# take top 100 rows for debugging
#betasSub <- celltypeNormbeta[1:100,]
#row <- colMeans(celltypeNormbeta)

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

# 22 = number of output for each row from runEWAS() (12 (interactionLM) + 9 (nullLM) + 1 (anova))
outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 22, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("GroupWT_coeff", "GroupWT_SE", "GroupWT_P", 
                    "Age_coeff", "Age_SE", "Age_P",
                    "SexM_coeff", "SexM_SE", "SexM_P",
                    "GroupWT:QCmetrics$Age_coeff", "GroupWT:QCmetrics$SE", "GroupWT:QCmetrics$Age_P",
                    
                    "nullGroupWT_coeff", "nullGroupWT_SE", "nullGroupWT_P", 
                    "nullAge_coeff", "nullAge_SE", "nullAge_P",
                    "nullSexM_coeff", "nullSexM_SE", "nullSexM_P",
                    
                    "anovoP") 

filePath <- paste0(projDir, "/3_analysis/results/", cellType, "EWASout.rdat")
save(outtab, file = filePath)



#----------------------------------------------------------------------#
# Explore results
#----------------------------------------------------------------------#


bonfP <- 0.05/nrow(outtab)

# for each col, how many resultsa are lower than bonfP

countSig <- function(x){
  sig <- sum(x < bonfP)
  return(sig)
}

apply(outtab[,1:22], 2, countSig) 
