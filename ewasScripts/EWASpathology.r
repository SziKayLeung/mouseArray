##---------------------------------------------------------------------#
##
## Title: Run EWAS
##
## Purpose of script: Run linear regression within each cell type
##                    
##                    This scripted is adapted from one co-authored
##                    by EJH and EMW for the MRC schizophrenia project
##
## Adapted by SKL to also run bulk methylation                
##
##---------------------------------------------------------------------#

message("Script processed: ", format(Sys.time(),'%A, %B %d, %Y %H:%M:%S'),".")


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#

# to run: sbatch..... 

# within cell type models

# null model to check gorup term improves the model

# excluded age (Jon's advice)



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(dplyr))


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

# if cell type = neun+ source functions without random terms
# if cell type = neun- source functions w/ random terms

runEWAS<-function(row, QCmetrics, deconvolution, species){

  if(isTRUE(deconvolution)){
  
    if(species == "mouse"){
    
       nullLM<-lm(row ~ QCmetrics$Sex + QCmetrics$DN + QCmetrics$NEUN)

       fullLM<-lm(row ~ QCmetrics$Pathology + QCmetrics$Sex + QCmetrics$DN + QCmetrics$NEUN + QCmetrics$Pathology*QCmetrics$DN)
       
      # extract case control main effect and sex effect
      return(c(summary(fullLM)$coefficients["QCmetrics$Pathology",c(1,2,4)],
               summary(fullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
               
               # extract sex effect
               summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
               
               # extract interaction terms
               summary(fullLM)$coefficients["QCmetrics$Pathology:QCmetrics$DN",c(1,2,4)],
               
               # run anova
               anova(fullLM, nullLM)[2,"Pr(>F)"])) 
    }else{
    
      nullLM<-lm(row ~ QCmetrics$Sex + QCmetrics$IRF8 + QCmetrics$NEUN + QCmetrics$SOX10 + QCmetrics$TN)

      fullLM<-lm(row ~ QCmetrics$Pathology + QCmetrics$Sex + QCmetrics$IRF8 + QCmetrics$NEUN + QCmetrics$SOX10 + QCmetrics$TN + QCmetrics$Pathology*QCmetrics$IRF8 + QCmetrics$Pathology*QCmetrics$SOX10)
       
      # extract case control main effect and sex effect
      return(c(summary(fullLM)$coefficients["QCmetrics$Pathology",c(1,2,4)],
               summary(fullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
               
               # extract sex effect
               summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
               
               # extract interaction terms
               summary(fullLM)$coefficients["QCmetrics$Pathology:QCmetrics$IRF8",c(1,2,4)],
               summary(fullLM)$coefficients["QCmetrics$Pathology:QCmetrics$SOX10",c(1,2,4)],
               
               # run anova
               anova(fullLM, nullLM)[2,"Pr(>F)"]))
  
    }
  
  
  }else{
    
    nullLM<-lm(row ~ QCmetrics$Sex)
    
    fullLM<-lm(row ~ QCmetrics$Pathology + QCmetrics$Sex)
    
    # extract case control main effect and sex effect
    return(c(summary(fullLM)$coefficients["QCmetrics$Pathology",c(1,2,4)],
             summary(fullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
             
             # extract sex effect
             summary(nullLM)$coefficients["QCmetrics$SexM",c(1,2,4)],
             
             # run anova
             anova(fullLM, nullLM)[2,"Pr(>F)"]))
  
  }             
  
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

# Check if third argument exists and is non-empty (i.e. whether to run pathology with CETYGO score as covariate, in which case input = CEYTGO input)
if (length(args) < 3 || args[3] == "" || is.na(args[3])) {
  deconvolution_input <- NULL
  deconvolution = FALSE
  species=NULL
  message("No deconvolution input provided.")
} else {
  deconvolution_input <- args[3]
  deconvolution = TRUE
  message("Deconvolution input: ", deconvolution_input)
  
  if (length(args) < 4 || args[4] == "" || is.na(args[4])) {
    message("Need 4th argument about species")
  }else{
    species <- args[4]
  }
  
}



# Added matching pathology column (not previously in the normalisedData)
normData<-file.path(projDir, "2_normalised/normalisedDataPathology.rdat")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(projDir)
message("Loading normalised data and QC metrics: normData")
load(normData)


# Check if 'Group' column exists
if ("Group" %in% colnames(QCmetrics)) {
  # If 'Group' exists, convert it to a factor
  QCmetrics$Group <- as.factor(QCmetrics$Group)
} else if ("Phenotype" %in% colnames(QCmetrics)) {
  # If 'Group' doesn't exist but 'Phenotype' exists, create 'Group' from 'Phenotype'
  QCmetrics$Group <- as.factor(QCmetrics$Phenotype)
  message("Using 'Phenotype' as 'Group'")
} else if ("Genotype" %in% colnames(QCmetrics)) {
  # If 'Group' and 'Phenotype' don't exist, but 'Genotype' exists, create 'Group' from 'Genotype'
  QCmetrics$Group <- as.factor(QCmetrics$Genotype)
  message("Using 'Genotype' as 'Group'")
} else {
  # If neither 'Group', 'Phenotype', nor 'Genotype' exist
  message("No valid group column found")
}

# subset on cell type if 2nd argument provided
if(cellType != "bulk"){
  print(paste0("running EWAS on ", cellType, " cell type..."))
  ## subset to cell type samples
  QCmetrics<-QCmetrics[which(QCmetrics$Cell_Type == cellType),]
  
  # subset beta matrix to cell type specific samples
  celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
}else{
  message("Running EWAS on bulk")
  celltypeNormbeta<-normBeta
} 

# use covariate chip position if bulk tissue 
if(cellType == "bulk"){
  QCmetrics$Batch <- as.factor(QCmetrics$Chip_Position)
}else{
  QCmetrics$Batch <- as.factor(QCmetrics$Batch)
}

QCmetrics$Sex <- as.factor(QCmetrics$Sex)

# check consistent number of samples in QC metrics and normalised data
if(length(setdiff(colnames(celltypeNormbeta),QCmetrics$Basename)) != 0){
  print("Normalised beta values have samples that are not in the QC metrics")
}
if(length(setdiff(QCmetrics$Basename, colnames(celltypeNormbeta))) != 0){
  print("Normalised beta values have samples that are not in the QC metrics")
}

# check the order of sample names in QC metrics and the normalised data column names are the same
if(isFALSE(identical(colnames(celltypeNormbeta), QCmetrics$Basename))){
  print("Column names of normalised beta and QC metrics Basename are not in the same order; Reordered QC metrics")
  QCmetrics <- QCmetrics[match(colnames(celltypeNormbeta), QCmetrics$Basename),]
}

# input CETYGO input if using it as covariate
if(!is.null(deconvolution_input)){

  message("Reading in deconvolution input file")
  deconvolution_input <- get(load(deconvolution_input))
  
  # check that there is a CETYGO column
  if(isFALSE("CETYGO" %in% colnames(deconvolution_input))){
    print("CETGO score needs to be in the input deconvolution file")
  }
  
  # check same samples are in both the QCmetric and in the CETYGO output
  if(length(setdiff(row.names(deconvolution_input), row.names(QCmetrics))) != 0){
    print(setdiff(row.names(deconvolution_input), row.names(QCmetrics)))
    print(row.names(QCmetrics))
    print("Samples in the deconvolution file that are not in the QC metric, please check and re-run")
  }
  
  if(length(setdiff(row.names(QCmetrics), row.names(deconvolution_input))) != 0){
    print("Samples in the QC metric that are not in the deconvolution input, please check and re-run")
  }

  # merge by row.names
  QCmetrics <- merge(QCmetrics, deconvolution_input, all = T, by = 0)
}


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

message("Deconvolution = ", deconvolution)
message("Deconvolution on the species dataset = ", species) 

if(isTRUE(deconvolution) & species == "mouse"){
  outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics, deconvolution, species), ncol = 13, byrow = TRUE)
  outputcolumnNames=c("Pathology_coeff", "Pathology_SE", "Pathology_P","SexM_coeff", "SexM_SE", "SexM_P", "nullSexM_coeff", "nullSexM_SE", "nullSexM_P", 
"Pathology:DN_coeff", "Pathology:DN_SE", "Pathology:DN_P", "anovoP")
}else if(isTRUE(deconvolution) & species == "human"){
  outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics, deconvolution, species), ncol = 16, byrow = TRUE)
  outputcolumnNames=c("Pathology_coeff", "Pathology_SE", "Pathology_P","SexM_coeff", "SexM_SE", "SexM_P", "nullSexM_coeff", "nullSexM_SE", "nullSexM_P", 
"Pathology:IRF8_coeff", "Pathology:IRF8_SE", "Pathology:IRF8_P", "Pathology:SOX10_coeff", "Pathology:SOX10_SE", "Pathology:SOX10_P", "anovoP")
}else{
 outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics, deconvolution, species), ncol = 10, byrow = TRUE)
 outputcolumnNames=c("Pathology_coeff", "Pathology_SE", "Pathology_P","SexM_coeff", "SexM_SE", "SexM_P", "nullSexM_coeff", "nullSexM_SE", "nullSexM_P", "anovoP")
}


rownames(outtab)<-rownames(celltypeNormbeta)
colnames(outtab)<-outputcolumnNames

if(isTRUE(deconvolution)){
  filePath <- paste0(projDir, "/3_analysis/results/", cellType, "EWASpathout_", species, "cellTypeDeconvolution.rdat")
  save(outtab, file = filePath)
}else{
  filePath <- paste0(projDir, "/3_analysis/results/", cellType, "EWASpathout.rdat")
  save(outtab, file = filePath)
}

message("Output:", filePath)
message("EWAS finished!")