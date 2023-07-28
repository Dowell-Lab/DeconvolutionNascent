#!/bin/Rscript
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=02:00:00
#SBATCH --mail-user=zama8258@colorado.edu # nolint

install.packages(c('tidyverse', 'ggthemes', 'e1071',
                   'LiblineaR', 'MASS', 'nnls', 'progress',
                   'corrplot', 'rstan'),
                 repos = 'https://cloud.r-project.org')
