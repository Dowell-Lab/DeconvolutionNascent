### decon_utils.R --- Deconvolution Utilities
##
## Filename: decon_utils.R
## Author: Zach Maas
## Created: Tue Nov 29 16:00:52 2022 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains code to make deconvolution experiments easier
## and more replicable across different files.
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
library('corrplot')
library('rstan')
library('scales')
library('glmnet')
library('Matrix')
library('expm')
library('memoise')
library('patchwork')

## install.packages(c('tidyverse', 'ggthemes', 'e1071',
##                    'LiblineaR', 'MASS', 'nnls',
##                    'progress', 'corrplot', 'rstan',
##                    'glmnet', 'expm', 'Matrix', 'expm',
##                    'memoise', 'gtools', 'scales', 'patchwork'))

## install.packages(c('tidyverse', 'ggthemes', 'e1071', 'LiblineaR', 'MASS', 'nnls', 'progress', 'corrplot', 'rstan', 'glmnet', 'expm', 'Matrix', 'expm', 'memoise', 'gtools', 'scales', 'patchwork'))

###################################
## Deconvolution Implementations ##
###################################

scale_coefs <- function(coefs) {
    coefs <- ifelse(coefs > 0, coefs, 0) ## Force terms to be positive
    coefs <- coefs / sum(coefs) ## Force terms to sum to one
    return(coefs)
}

## Wrapper to run a linear nu-svm model
nu_svm <- function(X, Y, nu) {
    model <- svm(X, y=Y, kernel='linear',
                 type='nu-regression',
                 C=1, nu=nu)
    coefs <- coef(model)[-1] ## Remove the intercept term
    coefs <- scale_coefs(coefs)
    return(coefs)
}

eps_svm <- function(X, Y, eps) {
    ## Alternative impl using liblinear for speed
    model <- LiblineaR(X, Y, type=11, svr_eps=0.1)
    coefs <- model$W[1:length(model$W)-1]
    coefs <- scale_coefs(coefs)
    return(coefs)
}

invert_matrix <- function(X, Y) {
    coefs <- ginv(X) %*% Y
    coefs <- scale_coefs(coefs)
    return(coefs)
}

decon_nnls <- function(X, Y) {
    res <- nnls(X, Y)
    coefs <- scale_coefs(res$x)
    return(coefs)
}

decon_lasso <- function(X, Y) {
    cv_model <- cv.glmnet(X, Y, alpha=1)
    best_lambda <- cv_model$lambda.min
    res <- glmnet(X, Y, alpha=1, lambda=best_lambda)
    coefs <- coef(res)[-1]
    coefs[is.nan(coefs)] <- 0
    coefs <- scale_coefs(coefs)
    return(coefs)
}

decon_relasso <- function(X, Y) {
    cv_model <- cv.glmnet(X, Y, alpha=1, relax=TRUE)
    best_lambda <- cv_model$lambda.min
    res <- glmnet(X, Y, alpha=1, lambda=best_lambda, relax=TRUE)
    coefs <- coef(res)[-1]
    coefs <- scale_coefs(coefs)
    return(coefs)
}

decon_ridge <- function(X, Y) {
    cv_model <- cv.glmnet(X, Y, alpha=0)
    best_lambda <- cv_model$lambda.min
    res <- glmnet(X, Y, alpha=0, lambda=best_lambda)
    coefs <- coef(res)[-1]
    coefs <- scale_coefs(coefs)
    return(coefs)
}

decon_elastic <- function(X, Y, alpha) {
    cv_model <- cv.glmnet(X, Y, alpha=alpha)
    best_lambda <- cv_model$lambda.min
    res <- glmnet(X, Y, alpha=alpha, lambda=best_lambda)
    coefs <- coef(res)[-1]
    coefs <- scale_coefs(coefs)
    return(coefs)
}

num_cores <- 4
options(mc.cores=num_cores)
rstan_options(auto_write = TRUE)
model <- 'src/linreg.stan'
model <- '/scratch/Users/zama8258/deconascent/src/linreg.stan'
decon_mcmc <- function(X, Y) {
    count_data <- list(x=X,y=Y,N=nrow(X),K=ncol(X))
    fit <- stan('src/linreg.stan', # nolint
                data = count_data, # nolint
                iter=2000, algorithm="NUTS", cores=num_cores)
    posterior <- get_posterior_mean(fit)
    size_fit <- ncol(posterior)
    coefs <- posterior[1:ncol(X)+1,size_fit]
    coefs <- scale_coefs(coefs)
    return(coefs)
}

## Function to infer optimal nu from a number of nu values
get_best_nu <- function(X, Y, nus) {
    errors <- c()
    for (curr_nu in nus) {
        model <- nu_svm(X, Y, curr_nu)
        error <- sum(abs(model$residuals))
        errors <- c(errors, error)
    }
    return(nus[which.min(errors)])
}

## Pruning algorithm
make_sorted_ratio_genelist <- function(counts_homo, homo_ratio) {
    ## Calculate ratio of sample / max(other samples) per ROI
    for (celltype in colnames(homo_ratio)) {
        ## print(paste0("Calculating ", celltype))
        homo_ratio[,colnames(counts_homo)==celltype] <-
            counts_homo[,colnames(counts_homo)==celltype] /
            apply(counts_homo[,colnames(counts_homo)!=celltype], 1, max)
    }
    ## Sort descending by the largest total fold-ratios
    homo_ratio_sorted <- homo_ratio[order(-apply(homo_ratio, 1, max)),]
    return(homo_ratio_sorted)
}
make_sorted_ratio_genelist_memo <- memoise(make_sorted_ratio_genelist)
find_top_genes_per_celltype <- function(counts_homo, n) {
    print(paste0("Calculating top genes for n=",n," genes"))
    ## Copy our data and apply 90 percentile thresholding
    homo_ratio <- data.frame(counts_homo) %>%
        dplyr::filter(if_all(everything(),
                             ~.x <= quantile(.x,.99))) %>%
        as.matrix()
    homo_copy <- homo_ratio
    ## Get the sorted genelist (memoized)
    homo_ratio_sorted <- make_sorted_ratio_genelist_memo(homo_copy, homo_ratio)
    ## Find the indices corresponding to the max celltype in each
    hrs_max_idx <- apply(homo_ratio_sorted, 1, which.max)
    ## For each celltype, pick the top n rois that are celltype specific
    i <- 0
    rois <- c()
    for (celltype in colnames(homo_ratio)) {
        i <- i + 1
        rois <- c(rois, names(head(hrs_max_idx[hrs_max_idx == i], n)))
        if (any(is.na(rois))) {
            print(paste0("Exhausted distinct ROIs for ", celltype))
        }
    }
    ## Turn these rois into indices of the original matrix
    subset_idxs <- match(rois, rownames(counts_homo))
    subset_idxs <- na.omit(subset_idxs)
    return(subset_idxs)
}

#######################
## Utility Functions ##
#######################

## RMS calculations
get_rms <- function(expt, true) {
    rms <- sqrt(sum((expt - true)^2))
}

## Read in a single sample
read_sample <- function(iter, data_dir) {
    sample <- read_delim(paste0(data_dir, iter, "_counts.txt"),
                         skip=1, show_col_types=FALSE)
    colnames(sample) <- c("Geneid", "Chr", "Start",
                          "End", "Strand", "Length",
                          paste0("run_", iter))
    sample <- sample %>% dplyr::select(c("Geneid", paste0("run_", iter)))
    return(sample)
}

## Unified function to read in homogenous samples from a directory
read_merged_samples <- function(data_dir) {
    sample_hct <- read_delim(paste0(data_dir, "SRR6789175.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_hela <- read_delim(paste0(data_dir, "SRR7010982.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_mcf7 <- read_delim(paste0(data_dir, "SRR4090102.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_k562 <- read_delim(paste0(data_dir, "SRR5364303.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_esc <- read_delim(paste0(data_dir, "SRR10669536.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_kasumi <- read_delim(paste0(data_dir, "SRR3713700.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_cd4 <- read_delim(paste0(data_dir, "SRR6780907.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_jurkat <- read_delim(paste0(data_dir, "SRR7616132.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_bj5ta <- read_delim(paste0(data_dir, "SRR10601203.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    sample_lcl <- read_delim(paste0(data_dir, "SRR1552485.sorted.unstranded.bidir_counts.txt"), show_col_types=FALSE)
    samples_merged <- inner_join(sample_hct, sample_hela,
                                 by="GeneID") %>%
        inner_join(sample_mcf7, by="GeneID") %>%
        inner_join(sample_k562, by="GeneID") %>%
        inner_join(sample_esc, by="GeneID") %>%
        inner_join(sample_kasumi, by="GeneID") %>%
        inner_join(sample_cd4, by="GeneID") %>%
        inner_join(sample_jurkat, by="GeneID") %>%
        inner_join(sample_bj5ta, by="GeneID") %>%
        inner_join(sample_lcl, by="GeneID") %>%
        transmute(GeneID = GeneID,
                  hct = `SRR6789175.sorted.sorted.bam`,
                  hela = `SRR7010982.sorted.sorted.bam`,
                  mcf7 = `SRR4090102.sorted.sorted.bam`,
                  k562 = `SRR5364303.sorted.sorted.bam`,
                  esc = `SRR10669536.sorted.sorted.bam`,
                  kasumi = `SRR3713700.sorted.sorted.bam`,
                  cd4 = `SRR6780907.sorted.sorted.bam`,
                  jurkat = `SRR7616132.sorted.sorted.bam`,
                  bj5ta = `SRR10601203.sorted.sorted.bam`,
                  lcl = `SRR1552485.sorted.sorted.bam`) #%>%
    return(samples_merged)
}

## Generate a result tibble for use in deconvolution loops
gen_res_tbl <- function() {
    res_tbl <- tibble(
        mixing_fraction  = numeric(),
        rms              = numeric(),
        hct_actual       = numeric(),
        hela_actual      = numeric(),
        mcf7_actual      = numeric(),
        k562_actual      = numeric(),
        esc_actual       = numeric(),
        kasumi_actual    = numeric(),
        cd4_actual       = numeric(),
        jurkat_actual    = numeric(),
        bj5ta_actual     = numeric(),
        lcl_actual       = numeric(),
        hct_estimated    = numeric(),
        hela_estimated   = numeric(),
        mcf7_estimated   = numeric(),
        k562_estimated   = numeric(),
        esc_estimated    = numeric(),
        kasumi_estimated = numeric(),
        cd4_estimated    = numeric(),
        jurkat_estimated = numeric(),
        bj5ta_estimated  = numeric(),
        lcl_estimated    = numeric()
    )
    return(res_tbl)
}

## Alternate implementation with added values
gen_multires_tbl <- function() {
    ## Returns a multiple result tibble
    ## For comparing results over many runs
    multires_tbl <- tibble(
        ## Metadata
        mixing_fraction  = numeric(),
        rms              = numeric(),
        num_regions      = numeric(),
        method           = factor(),
        ## Actual fractions
        hct_actual       = numeric(),
        hela_actual      = numeric(),
        mcf7_actual      = numeric(),
        k562_actual      = numeric(),
        esc_actual       = numeric(),
        kasumi_actual    = numeric(),
        cd4_actual       = numeric(),
        jurkat_actual    = numeric(),
        bj5ta_actual     = numeric(),
        lcl_actual       = numeric(),
        ## Estimated fractions
        hct_estimated    = numeric(),
        hela_estimated   = numeric(),
        mcf7_estimated   = numeric(),
        k562_estimated   = numeric(),
        esc_estimated    = numeric(),
        kasumi_estimated = numeric(),
        cd4_estimated    = numeric(),
        jurkat_estimated = numeric(),
        bj5ta_estimated  = numeric(),
        lcl_estimated    = numeric()
    )
    return(multires_tbl)
}

#######################
## Unified Interface ##
#######################

## This doesn't make sense except for in the context of the homo mixture
## Otherwise we need to keep the same filterset for the hetero, so no aio func
removeSparseRows <- function(counts) {
    filter_regions <- rowSums(counts[2:length(counts)] != 0) > 1
    counts_filtered <- counts[filter_regions, ]
    return(counts_filtered)
}

## How do I do this on both sets of counts cleanly??
subsetROIs <- function(counts, n) {
    filter_set <- find_top_genes_per_celltype(counts, n)
    counts_filtered <- counts[filter_set, ]
}

## So, probably combine the two filtering functions to generate a
## single set of regions? Is this necessary? I don't think so for
## small n since we'll automagically exclude all-zero regions for
## having zero variance

## Please write a docstring for the following R code:
## Unified function to estimate mixing fractions
estimateMixing <- function(counts_homo, counts_hetero, ...) {
    args <- list(...)
    model_choice <- args$model_choice
    param <- args$param
    if (is.null(model_choice)) {
        warning("Model not specified, defaulting to nnls")
        model_choice <- "nnls"
    }
    ## How to only assign parameter if model choice requires it
    if (is.null(param) &&
        model_choice %in% c("nu_svm", "eps_svm", "elastic")) {
        warning(paste0("Model ", model_choice, " requires parameter selection",
                       " Setting to 0.1 by default."))
        param <- 0.1
    }
    if (model_choice == "nu_svm") {
        model <- function(X, Y) {
            nu_svm(X, Y, param)
        }
    } else if (model_choice == "eps_svm") {
        model <- function(X, Y) {
            eps_svm(X, Y, param)
        }
    } else if (model_choice == "invert") {
        model <- function(X, Y) {
            invert_matrix(X, Y)
        }
    } else if (model_choice == "nnls") {
        model <- function(X, Y) {
            decon_nnls(X, Y)
        }
    } else if (model_choice == "mcmc") {
        model <- function(X, Y) {
            decon_mcmc(X, Y)
        }
    }
    else if (model_choice == "lasso") {
        model <- function(X, Y) {
            decon_lasso(X, Y)
        }
    }
    else if (model_choice == "ridge") {
        model <- function(X, Y) {
            decon_ridge(X, Y)
        }
    }
    else if (model_choice == "relasso") {
        model <- function(X, Y) {
            decon_relasso(X, Y)
        }
    }
    else if (model_choice == "elastic") {
        model <- function(X, Y) {
            decon_elastic(X, Y, param)
        }
    }
    coefs <- model(counts_homo, counts_hetero)
    return(coefs)
}

deconvolveSystem <- function(counts_homo, counts_hetero, ...) {
    ## Parse function arguments
    args <- list(...)
    model_choice <- args$model_choice
    param <- args$param
    n <- args$n
    ## Prune the system if requested
    if (is.null(n)) {
        n <- min(100, 0.01*nrow(counts_homo))
        print(paste0("Guessed n = ", n))
    }
    homo_rs <- subsetROIs(counts_homo)
    hetero_rs <- subsetROIs(counts_hetero)
    ## Run the model
    coefs <- estimateMixing(homo_rs, hetero_rs, model_choice, param)
    return(coefs)
}


######################################################################
### decon_utils.R ends here
