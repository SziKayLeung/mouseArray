#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=EWAS.err # error file
#SBATCH --output=EWAS.log # output file
#SBATCH --job-name=EWAS


#------------------------------------------------------

# 1. command line argument input is:

#   1. <filepath/to/projectFolder>
#   2. cellType
# e.g. to run sbatch EWASjobSubmission.sh <filepath/to/projectFolder> "NEUNpos"

#-----------------------------------------------------

## print start date and time
echo Job started on:
date -u
JOBNAME="EWAS"

echo Job sumbitted from:
echo $(pwd)

# Move the user to the project directory
cd $1


## load modules
module load Pandoc
module load R/4.2.1-foss-2022a

Rscript /lustre/projects/Research_Project-MRC148213/lsl693/scripts/mouseArray/ewasScripts/EWAS.r $1 $2

## print finish date and time
echo Job finished on:
date -u