#!/bin/env Rscript
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=16:00:00
#SBATCH --mail-user=zama8258@colorado.edu # nolint
### run_svm_iter.R --- Run an SVM iteration
##
## Filename: run_svm_iter.R
## Author: Zach Maas
## Created: Wed Oct 19 13:22:03 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## This script runs a single SVR iteration for our titration
## experiments, the code is split up this way to easily parallelize
## things since running them iteratively is too slow to iterate how we
## want to.
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

## Load our shared utilities
suppressMessages(source("./decon_utils.R"))

suppressPackageStartupMessages(library("optparse")) ## Argument Parser
## suppressPackageStartupMessages(library("tidyverse")) ## For data munging
## suppressPackageStartupMessages(library("e1071")) ## LibSVM Wrapper
## suppressPackageStartupMessages(library("LiblineaR")) ## LibLinear Wrapper

## Argument parsing
option_list <- list(
    make_option(c("-i", "--iter"),
                default = 0.5,
                help = "Model iteration to run [default %default]"
                ),
    make_option(c("-m", "--model"),
                default = "nnls",
                help = "Model to run [default %default]"
                ),
    make_option(c("-p", "--param"),
                default = 0.5,
                help = "Parameter for SVR / models with a choosable param [default %default]"
                ),
    make_option(c("-o", "--output"),
                default = ".",
                help = "Output file path [default \"%default\"]"
                )
)
opt <- parse_args(OptionParser(option_list = option_list))
iter <- opt$iter
model_choice <- opt$model
param <- opt$param
outfile <- opt$output

## Pick which model to use
if (model_choice == "nu_svm") {
    print("Using nu_svm")
    model <- function(X, Y) {
        nu_svm(X, Y, param)
    }
} else if (model_choice == "eps_svm") {
    print("Using eps_svm")
    model <- function(X, Y) {
        eps_svm(X, Y, param)
    }
} else if (model_choice == "invert") {
    print("Using invert")
    model <- function(X, Y) {
        invert_matrix(X, Y)
    }
} else if (model_choice == "nnls") {
    print("Using nnls")
    model <- function(X, Y) {
        decon_nnls(X, Y)
    }
} else if (model_choice == "mcmc") {
    print("Using mcmc")
    model <- function(X, Y) {
        decon_mcmc(X, Y)
    }
}

## Load our data in
## load("../test_data.Rdata")
load("/scratch/Users/zama8258/deconascent/test_data.Rdata")
## Provides:
## counts_homo
## samples_sub
## mixing_fractions

## Generate our simulated matrices
mixing_fraction <- mixing_fractions[iter]
coefs_true <- seq(mixing_fraction, 1, length.out = 9)
## For holdout experiments...
                                        # coefs_true[5] <- 0
                                        # coefs_true[9] <- 0
coefs_true <- coefs_true / sum(coefs_true)
counts_hetero <- as.matrix(samples_sub[iter])
if (model_choice == "mcmc") {
    ## STAN wants an array as input
    counts_hetero <- array(counts_hetero)
}

testing <- FALSE
if (testing) {
    n_row <- 1000
    counts_homo <- counts_homo[1:n_row,1:ncol(counts_homo)]
    counts_hetero <- counts_hetero[1:n_row]
}

## Run Model
coefs <- model(counts_homo, counts_hetero)


rms <- sqrt(sum((coefs - coefs_true)^2))
write(
    x = paste0(coefs, collapse="\t"),
    file = outfile
)

######################################################################
### run_svm_iter.R ends here
