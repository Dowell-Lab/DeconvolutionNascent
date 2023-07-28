### specs_from_scratch.R --- Implementing SPECS from scratch
##
## Filename: specs_from_scratch.R
## Author: Zach Maas
## Created: Mon Jun 26 09:57:40 2023 (-0600)
##
######################################################################
##
### Commentary:
##
## The math in the SPECS paper is not very complicated, but it is not
## super well written and the implementation provided by the authors
## lacks comments. This file is an attempt to implement the SPECS
## score from scratch in R. It is not intended to be complete.
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

## This implementation assumes 1 sample per tissue, so the impl is
## simplified compared to the paper

## Optional benchmarking
## library('microbenchmark')

## SPECS for a single sample
calc_specs_singlesample <- function(data, weight) {
    ## Generate a probability matrix
    probs<- matrix(0, nrow=nrow(data), ncol=ncol(data))
    ## Iterate over every feature
    for (i in seq_len(nrow(data))) {
        ## For every sample
        for (j in seq_len(ncol(data))) {
            ## Calculate the probability that the feature is greater than
            ## every other feature efficiently
            probs[i,j] <- sum(data[i,j] > data[i,][-j]) /
                (ncol(data)-1)
        }
    }
    specs <- probs / weight
    rownames(specs) <- rownames(data)
    colnames(specs) <- colnames(data)
    return(specs)
}

tf_idf <- function(data) {
    data_smoothed <- data + 1
    ## Calculate the term frequency
    tf <- data_smoothed / colSums(data_smoothed)
    ## Calculate the inverse document frequency
    idf <- log(nrow(data_smoothed) / rowSums(data_smoothed))
    ## Calculate the tf-idf
    tf_idf <- tf * idf
    rownames(tf_idf) <- rownames(data)
    colnames(tf_idf) <- colnames(data)
    return(tf_idf)
}

## Generate a fake matrix
test_data <- matrix(rpois(25*10000, 100), ncol=25)
colSums(calc_specs_singlesample(test_data, 1) > 0.95)

## lengths <- roi_info$length[filter_regions]
## reads_per_length <- counts_homo
## for (i in seq_len(ncol(reads_per_length))) {
##     reads_per_length[,i] <- reads_per_length[,i] / lengths
## }
## tpm_homo <- 10^6 * reads_per_length / colSums(reads_per_length)


homo_tf_idf <- tf_idf(counts_homo[subset_idxs,])
homo_tf_idf_df <- homo_tf_idf %>%
    as_tibble() %>%
    pivot_longer(colnames(.))

ggplot(homo_tf_idf_df) +
    geom_density(aes(x=value,color=name)) +
    scale_x_log10()

ggplot(homo_tf_idf_df) +
    stat_ecdf(aes(x=value,color=name)) +
    xlim(0,5e-5)

######################################################################
### specs_from_scratch.R ends here
