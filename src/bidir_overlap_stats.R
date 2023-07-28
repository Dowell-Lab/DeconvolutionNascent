### bidir_overlap_stats.R --- Analyze some bidir overlap stats
##
## Filename: bidir_overlap_stats.R
## Author: Zach Maas
## Created: Wed Feb  8 14:19:48 2023 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains some exploratory analysis about overlaps that
## GIGGLE finds.
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:

library('tidyverse')
library('data.table')
library('sparsepca')
library('distances')

## WARNING - This is a big file, needs ~4Gb of RAM probably
data_path <- '~/Dropbox/phd/research/dna_lab/deconascent/data/table.txt.gz'
mask <- fread(data_path)

## Get some basic analysis stats and make sparse
n_samples <- length(colnames(mask)) - 3
regions <- mask[,1:3]
regions[,roi:=paste0(regions[,1], regions[,2], regions[,3])]
mask <- mask[,4:n_samples]

per_sample_bidirs<- colSums(mask)
ggplot() +
    geom_density(aes(x=per_sample_bidirs))

per_bidir_samples <- rowSums(mask)
ggplot() +
    geom_density(aes(x=per_bidir_samples))

## distance_matrix <- distances::distances(mask)
## pca <- cmdscale(distance_matrix)

srr_samples <- c("SRR7010982",
                 "SRR4090102",
                 "SRR5364303",
                 "SRR10669536",
                 "SRR3713700",
                 "SRR6780907",
                 "SRR7616132",
                 "SRR10601203")
mask_merge <- data.table(overlaps=rep(0,nrow(sample_mask)))
for (srr in srr_samples) {
    sample_mask <- mask[,..srr]
    mask_merge <- mask_merge + sample_mask
}
cbind(regions, mask_merge)[overlaps > 0]

######################################################################
### bidir_overlap_stats.R ends here
