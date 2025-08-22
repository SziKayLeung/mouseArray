#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion

# run on command line

## load modules
module load Pandoc
module load R/4.2.1-foss-2022a

## input paths
projDir=/lustre/projects/Research_Project-MRC190311/DNAm/dementiaMouse/APPKI_bulk
scriptPath=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/mouseArray/ewasScripts/EWASpathology.r
mouseDeconvolution=/lustre/projects/Research_Project-MRC190311/DNAm/dementiaMouse/APPKI_bulk/3_analysis/CETYGO/bulk_CETYGO_predProp.RData
humanDeconvolution=/lustre/projects/Research_Project-MRC190311/DNAm/dementiaMouse/APPKI_bulk/3_analysis/CETYGO/mouse_human_CETYGO_predProp.RData

# run EWAS without CETYGO covariate 
Rscript ${scriptPath} ${projDir} bulk &> EWASpath.log

# run EWAS with CETYGO score as covariate (mouse FANS as training data, mouse bulk as testing data)
Rscript ${scriptPath} ${projDir} bulk ${mouseDeconvolution} mouse &> ${projDir}/3_analysis/results/EWASpathMouseCetygo.log

# run EWAS with CETYGO score as covariate (human FANS as training data, mouse bulk as testing data)
Rscript ${scriptPath} ${projDir} bulk ${humanDeconvolution} human &> ${projDir}/3_analysis/results/EWASpathHumanCetygo.log