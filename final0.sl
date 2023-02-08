#!/bin/tcsh
#SBATCH --job-name=fishbodysize0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=48:00:00
#SBATCH --error=efin0.err
#SBATCH --output=ofin0.out
#SBATCH --mail-type=END
#SBATCH --mail-user=oconnoyq@bc.edu

cd /scratch/oconnoyq
module load R/4.2.1gnu9.2.0
R -e 'setwd("/mmfs1/data/oconnoyq/R/x86_64-pc-linux-gnu-library/4.2")'
Rscript final0.R
