##---------------------------------------------------------------------#
##
## Title: Cluster cell types
##
## Purpose of script: confirm cell type labels are accurate using a variety
##                    of methods
##
##                    Note the calMouseMethQCmetrics.r script needs to be 
##                    have been run prior to this one
##
##                    The CellTypeQC.rmd script uses these metrics to output
##                    a html report
##
##                    This script is an adaptation of the one from brainFans
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#


# parameters and relative paths etc are loaded from the config.r 
# file in the project folder

#----------------------------------------------------------------------#
# # STOP IF NOT REQUESTED TO RUN THIS SCRIPT
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]

setwd(dataDir)

source("config.r")

if(!ctCheck){
	quit(save = "no", status = 0)
}


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

print("loading packages...")
library(MatrixGenerics)
library(ENmix)
library(data.table)


#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#

# N.B. Clustering is only performed on samples that have passed the following:
# intensity, pfilt and bscon checks

refDir <- args[2]

manifest <- file.path(refDir, "MouseMethylation-12v1-0_A2.csv")
normDir <- file.path(dataDir, "2_normalised")
QCDir <- file.path(normDir, "/QC")  

# these files are not currently filtered...
load(file = file.path(QCDir, "mraw.rdat"))
betas <- getB(mraw)

load(file = file.path(QCDir, "QCmetrics.rdat"))

man <- fread(manifest, skip=7, fill=TRUE, data.table=F)

QCSum <- read.csv(file.path(QCDir, "passQCStatusStage1AllSamples.csv"), stringsAsFactors = F)

passQC<-QCSum$Basename[QCSum[,"PassQC1"]]

QCmetrics.all<-QCmetrics
QCmetrics<-QCmetrics[match(passQC, QCmetrics$Basename),]

rawbetas<-betas[,match(passQC, colnames(betas))]

auto.probes<-man$IlmnID[man$CHR != "X" & man$CHR != "Y" & man$CHR != "MT"]

rawbetas<-rawbetas[row.names(rawbetas) %in% auto.probes,]

cellTypes<-unique(QCmetrics$Cell_Type)
cellTypes<-cellTypes[!is.na(cellTypes)]
cellTypes<-sort(cellTypes)

# filter out NAs
rawbetas<-na.omit(rawbetas)

# filter out SNPs
betas<-rawbetas[-grep("rs", rownames(rawbetas)),]




#----------------------------------------------------------------------#
# CALCULATE PCA
#----------------------------------------------------------------------#

print("Calculating PCs from beta matrix")

pca <- prcomp(t(betas))
betas.scores = pca$x
colnames(betas.scores) = paste(colnames(betas.scores), '_betas', sep='')
betas.pca<-pca$sdev^2/sum(pca$sdev^2)
betas.scores<-betas.scores[,which(betas.pca > 0.01)]

#----------------------------------------------------------------------#
# IDENTIFY OUTLIERS
#----------------------------------------------------------------------#

# Following methods rely on the ability to define an average profile of a cell type so need to check if any outliers first
print("Identifying outliers")

# calculate Studentized scores for each sample compared to it's cell type mean
studentPCA<-list()
studentPCA[[1]]<-matrix(NA, nrow = nrow(QCmetrics), ncol = ncol(betas.scores))
for(i in 1:nrow(QCmetrics)){
  keep<-QCmetrics$Basename[which(QCmetrics$Cell_Type == QCmetrics$Cell_Type[i])]
  ## exclude itself
  keep<-keep[keep != QCmetrics$Basename[i]]
  if(length(keep) > 1){
    cell.means<-colMeans(betas.scores[keep,])
    cell.sd<-apply(betas.scores[keep,], 2, sd)
    studentPCA[[1]][i,]<-c((betas.scores[i,]-cell.means)/cell.sd)
  }
}


# identify outliers from first 2 PCs only
studentCTOutlier<-rowMaxs(abs(studentPCA[[1]][,1:2])) > studentThres
QCmetrics$maxStudentPCA<-rowMaxs(abs(studentPCA[[1]][,1:2]))
QCmetrics$studentCTOutlier<-studentCTOutlier


#----------------------------------------------------------------------#
# CALCULATE CELL TYPE AVERAGE PROFILES
#----------------------------------------------------------------------#

# for each sample quantify how similar to cell type average it is
cellMeanPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier,], by = list(QCmetrics$Cell_Type[!QCmetrics$studentCTOutlier]), mean, na.rm = TRUE)
rownames(cellMeanPCA)<-cellMeanPCA[,1]
cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

cellSDPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier,], by = list(QCmetrics$Cell_Type[!QCmetrics$studentCTOutlier]), sd, na.rm = TRUE)
rownames(cellSDPCA)<-cellSDPCA[,1]
cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

maxSD<-rep(NA, nrow(QCmetrics))
for(i in 1:nrow(QCmetrics)){
  if(QCmetrics$Cell_Type[i] %in% colnames(cellMeanPCA)){
    maxSD[i]<-max(abs(betas.scores[i,1:2]-cellMeanPCA[1:2,QCmetrics$Cell_Type[i]])/cellSDPCA[1:2,QCmetrics$Cell_Type[i]])
  } 
}

#----------------------------------------------------------------------#
# CALCULATE INDIVIDUAL FACS SCORE
#----------------------------------------------------------------------#

keepCols<-c("Individual_ID", projVar)
keepCols<-keepCols[keepCols %in% colnames(QCmetrics)]
uniqueIDs<-unique(QCmetrics[,keepCols])
indFACSEff<-indFACSEff<-aggregate(maxSD[which(QCmetrics$Cell_Type != "Total")], by = list(QCmetrics$Individual_ID[which(QCmetrics$Cell_Type != "Total")]), FUN = median, na.rm = TRUE)
nFACs<-table(QCmetrics$Individual_ID[QCmetrics$Cell_Type != "Total"])

uniqueIDs<-cbind(uniqueIDs, indFACSEff$x[match(uniqueIDs$Individual_ID, as.character(indFACSEff$Group.1))], as.numeric(nFACs[uniqueIDs$Individual_ID]))
colnames(uniqueIDs)<-c(keepCols, "FACsEffiency", "nFACS")

write.csv(uniqueIDs, paste0(QCDir, "/IndividualFACsEffciencyScores.csv"))

# exclude individuals who hd very poor FACS sorts
QCmetrics$passFACS<-QCmetrics$Individual_ID %in% uniqueIDs$Individual_ID[which(uniqueIDs$FACsEffiency < 5)]


#----------------------------------------------------------------------#
# REPEAT OUTLIER DETECTION
#----------------------------------------------------------------------#

# calculate Studentized scores for each sample compared to it's cell type mean
studentPCA[[2]]<-matrix(NA, nrow = nrow(QCmetrics), ncol = ncol(betas.scores))
for(i in 1:nrow(QCmetrics)){
  if(QCmetrics$passFACS[i]){
    keep<-QCmetrics$Basename[which(QCmetrics$Cell_Type == QCmetrics$Cell_Type[i] & QCmetrics$passFACS)]
    # exclude itself
    keep<-keep[keep != QCmetrics$Basename[i]]
    if(length(keep) > 1){
      cell.means<-colMeans(betas.scores[keep,])
      cell.sd<-apply(betas.scores[keep,], 2, sd)
      studentPCA[[2]][i,]<-c((betas.scores[i,]-cell.means)/cell.sd)
    }
  }
}

# identify outliers from first 2 PCs only
studentPCA[[2]][is.na(studentPCA[[2]])]<-100
studentCTOutlier<-rowMaxs(abs(studentPCA[[2]][,1:2])) > studentThres
QCmetrics$maxStudentPCA2<-rowMaxs(abs(studentPCA[[2]][,1:2]))
QCmetrics$maxStudentPCA2[which(QCmetrics$maxStudentPCA2 == 100)]<-NA
QCmetrics$studentCTOutlier2<-studentCTOutlier


#----------------------------------------------------------------------#
# COMPARE TO CELL TYPE AVERAGE
#----------------------------------------------------------------------#

cellMeanPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell_Type[!QCmetrics$studentCTOutlier2]), mean, na.rm = TRUE)
rownames(cellMeanPCA)<-cellMeanPCA[,1]
cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))

cellSDPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell_Type[!QCmetrics$studentCTOutlier2]), sd, na.rm = TRUE)
rownames(cellSDPCA)<-cellSDPCA[,1]
cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))

maxSD<-rep(NA, nrow(QCmetrics))
for(i in 1:nrow(QCmetrics)){
  if(QCmetrics$Cell_Type[i] %in% colnames(cellMeanPCA)){
    maxSD[i]<-max(abs(betas.scores[i,1:2]-cellMeanPCA[1:2,QCmetrics$Cell_Type[i]])/cellSDPCA[1:2,QCmetrics$Cell_Type[i]])
  }
}

QCmetrics<-cbind(QCmetrics, maxSD)

# use Mahalanobis to calculate distance with cell type medians from PCAs 
# exclude outliers identified above from calculation of median profile


cellMedPCA<-aggregate(betas.scores[!QCmetrics$studentCTOutlier2,], by = list(QCmetrics$Cell_Type[!QCmetrics$studentCTOutlier2]), median, na.rm = TRUE)
rownames(cellMedPCA)<-cellMedPCA[,1]
cellMedPCA<-as.matrix(t(cellMedPCA[,-1]))

## calc cell type covariance matrices
cov_sigmaPCA<-list()
for(each in colnames(cellMedPCA)){
  ## only possible if > 1 sample
  if(sum(QCmetrics$Cell_Type == each & !QCmetrics$studentCTOutlier2) > 1){
    cov_sigmaPCA[[each]]<-cov(betas.scores[which(QCmetrics$Cell_Type == each & !QCmetrics$studentCTOutlier2),1:2], use = "p")
  }
}

mahDistPCA<-matrix(data = 10^9, ncol = ncol(cellMedPCA), nrow = nrow(QCmetrics))
colnames(mahDistPCA)<-colnames(cellMedPCA)
for(each in colnames(cellMedPCA)){
  if(!is.null(cov_sigmaPCA[[each]]) && det(cov_sigmaPCA[[each]]) > 1){
    mahDistPCA[,each]<-mahalanobis(betas.scores[,1:2], cellMedPCA[,each], cov_sigmaPCA[[each]], na.rm = TRUE,tol=1e-20)
  }
}

closestCellTypePCA<-colnames(mahDistPCA)[unlist(apply(mahDistPCA, 1, which.min))]
closestLabelledCellType<-QCmetrics$Cell_Type == closestCellTypePCA

# As two different antibodies used for neuronal cells used, accept if either predicted as the other
neunIndex<-which(QCmetrics$Cell_Type %in% neunCT)
closestLabelledCellType[neunIndex[closestCellTypePCA[neunIndex] %in% neunCT]]<-TRUE

# SATB2- are a composition of non neuronal cells; accept so long as predicted as non-neuronal
satb2NegIndex<-which(QCmetrics$Cell_Type == "SATB2-")
closestLabelledCellType[satb2NegIndex[!closestCellTypePCA[satb2NegIndex] %in% neunCT]]<-TRUE

QCmetrics<-cbind(QCmetrics, closestCellTypePCA, closestLabelledCellType)
write.csv(QCmetrics[which(closestLabelledCellType == "FALSE"),], paste0(QCDir, "/SamplesPredictedDiffCellTypePCAMahDist.csv"))

#----------------------------------------------------------------------#
# COMPARE TO CELL TYPE POLYTOPE
#----------------------------------------------------------------------#

# additionally check sample overlaps with mean profile of that cell type
# again exclude outliers identified above from calculating the average profile
print("Confirming sample overlap with group polytope")
outlierThres<-c(1.5,2,2.5,3)
pcaClassify<-list("predictCellType" = matrix(data = NA, nrow = nrow(QCmetrics), ncol = length(outlierThres)), "withinSDMean" = matrix(data = NA, nrow = nrow(QCmetrics), ncol = length(outlierThres)))
for(thres in outlierThres){
  studentCTOutlier<-QCmetrics$maxStudentPCA2 > thres
  
  cellMeanPCA<-aggregate(betas.scores[!studentCTOutlier,], by = list(QCmetrics$Cell_Type[!studentCTOutlier]), mean, na.rm = TRUE)
  rownames(cellMeanPCA)<-cellMeanPCA[,1]
  cellMeanPCA<-as.matrix(t(cellMeanPCA[,-1]))
  
  cellSDPCA<-aggregate(betas.scores[!studentCTOutlier,], by = list(QCmetrics$Cell_Type[!studentCTOutlier]), sd, na.rm = TRUE)
  rownames(cellSDPCA)<-cellSDPCA[,1]
  cellSDPCA<-as.matrix(t(cellSDPCA[,-1]))
  
  lowerBound<-cellMeanPCA-nSDThres*cellSDPCA
  upperBound<-cellMeanPCA+nSDThres*cellSDPCA
  
  
  # assumes cell types are distinct so just consider set of non-overlappping cell types
  pcaClassifyDistinctCT<-rep("", nrow(QCmetrics))
  for(each in predDistinctCT){
    if(each %in% colnames(upperBound)){
      indexClassify<-which(betas.scores[,1] < upperBound[1,each] & betas.scores[,1] > lowerBound[1,each] & betas.scores[,2] < upperBound[2,each] & betas.scores[,2] > lowerBound[2,each])
      pcaClassifyDistinctCT[indexClassify]<-paste(pcaClassifyDistinctCT[indexClassify], each, sep = ";")
    }
  }
  
  pcaClassify[["predictCellType"]][,match(thres, outlierThres)]<-pcaClassifyDistinctCT
  
  # can overlap with multiple CT so check if it if within 2SD of mean of it's labelled CT
  withinSDMean<-rep(FALSE, nrow(QCmetrics))
  for(i in 1:nrow(QCmetrics)){
    if(QCmetrics$Cell_Type[i] %in% colnames(upperBound)){
      withinSDMean[i]<-betas.scores[i,1] < upperBound[1,QCmetrics$Cell_Type[i]] & betas.scores[i,1] > lowerBound[1,QCmetrics$Cell_Type[i]] & betas.scores[i,2] < upperBound[2,QCmetrics$Cell_Type[i]] & betas.scores[i,2] > lowerBound[2,QCmetrics$Cell_Type[i]]
    }
  }
  
  pcaClassify[["withinSDMean"]][,match(thres, outlierThres)]<-withinSDMean
}

# pull out results for specified threshold
withinSDMean<-pcaClassify[["withinSDMean"]][,match(studentThres, outlierThres)]
pcaClassifyDistinctCT<-pcaClassify[["predictCellType"]][,match(studentThres, outlierThres)]

write.csv(QCmetrics[which(withinSDMean == "FALSE"),],paste0(QCDir, "/SamplesPCAOutlierFromCellType.csv"))

#----------------------------------------------------------------------#
# CLASSIFY CORRECT CELL TYPE
#----------------------------------------------------------------------#

# exclude individuals with sub optimal FACs
pcaClassifyDistinctCT[!QCmetrics$passFACS]<-FALSE
withinSDMean[!QCmetrics$passFACS]<-FALSE
QCmetrics<-cbind(QCmetrics, pcaClassifyDistinctCT,withinSDMean)

# keep all TOTAL samples
closestLabelledCellType[which(QCmetrics$Cell_Type == "Total")]<-TRUE
QCmetrics$predLabelledCellType<-closestLabelledCellType
withinSDMean[which(QCmetrics$Cell_Type == "Total")]<-TRUE
QCmetrics$withinSDMean<-withinSDMean

passCTCheck<-withinSDMean
passCTCheck[which(QCmetrics$Cell_Type == "Total")]<-TRUE
QCmetrics$passCTCheck<-passCTCheck


#----------------------------------------------------------------------#
# SAVE AND CLOSE
#----------------------------------------------------------------------#

save(betas.scores, mahDistPCA, pcaClassify, studentPCA, file = paste0(QCDir,"/PCAAcrossAllCellTypes.rdata"))

# add outcome to qc summary
QCmetrics<-QCmetrics[match(QCSum$Basename, QCmetrics$Basename),]

# add in data for samples excluded thus far
QCmetrics[which(is.na(QCmetrics$withinSDMean)),colnames(QCmetrics.all)]<-QCmetrics.all[which(is.na(QCmetrics$withinSDMean)),]

#QCSum<-QCSum[,-ncol(QCSum)]
QCSum<-cbind(QCSum, QCmetrics$passFACS, QCmetrics$passCTCheck,  QCmetrics$passFACS & QCmetrics$passCTCheck)
# retain TOTAL samples
QCSum[which(QCSum$Cell_Type == "Total" & QCSum$passQCS2 ==TRUE),ncol(QCSum)]<-TRUE
colnames(QCSum)[(ncol(QCSum)-2):ncol(QCSum)]<-c("passCTCheck", "passFACS", "passQCS3")

write.csv(QCSum, paste0(QCDir,"/passQCStatusStage3AllSamples.csv"))
write.csv(QCmetrics, paste0(QCDir,"/QCMetricsPostCellTypeClustering.csv"))

save(QCmetrics, file=file.path(QCDir, "QCmetricsPostCellTypeChecks.rdat"))

