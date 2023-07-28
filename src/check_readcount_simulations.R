### check_readcount_simulations.R --- check readcount shuffles
##
## Filename: check_readcount_simulations.R
## Author: Zach Maas
## Created: Thu Nov  3 14:46:52 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to analyze iterations of our deconvolution
## model made by randomly subsampling reads from the cram stage of the
## file instead of mixing reads from the final counts tables.
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
library('ggthemes')
library('e1071') ## LibSVM Wrapper
library('LiblineaR') ## LibSVM Wrapper
library('MASS')
library('nnls')
library('progress') ## Progress bars for slow stuff

read_sample <- function(iter) {
    data_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/simulated_counts/merged_"
    sample <- read_delim(paste0(data_dir, iter, "_counts.txt"),
                         skip=1, show_col_types=FALSE)
    colnames(sample) <- c("Geneid", "Chr", "Start",
                          "End", "Strand", "Length",
                          paste0("run_", iter))
    sample <- sample %>% dplyr::select(c("Geneid", paste0("run_", iter)))
    return(sample)
}

samples <- read_sample(1)
for (sample_iter in seq(2,99)) {
    sample <- read_sample(sample_iter)
    samples <- inner_join(samples, sample, by="Geneid")
}

######################################################################
### check_readcount_simulations.R ends here
