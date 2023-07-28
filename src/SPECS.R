### SPECS.R --- Implementation of the SPECS algorithm
##
## Filename: SPECS.R
## Author: Zach Maas
## Created: Wed Jun 21 12:18:47 2023 (-0600)
##
######################################################################
##
### Commentary:
##
## The SPECS algorithm is poorly documented and written in Python, so
## provided here is a reimplementation in R. The original code can be
## found at
## https://github.com/celineeveraert/SPECS/blob/master/SPECs_onlyp_example.py
##
## I have attempted here to provide a cleaner, more readable and
## documented interface than that provided by the original authors.
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

library(tidyverse)

## Function to calculate TPM values
calculate_tpm <- function(countdata) {
    tpm <- countdata %>%
        select(-c(GeneName, Length)) %>%
        mutate(across(everything(), as.numeric)) %>%
        column_to_rownames("GeneID") %>%
        as.matrix() %>%
        t() %>%
        as.data.frame() %>%
        mutate(across(everything(), function(x) x / (countdata$Length / 1000))) %>%
        as.matrix() %>%
        t() %>%
        as.data.frame() %>%
        mutate(across(everything(), function(x) x / sum(x) * 1e6)) %>%
        as.matrix() %>%
        t() %>%
        as.data.frame() %>%
        drop_na() %>%
        rownames_to_column("GeneID")

    return(tpm)
}

## Function to determine unique tissues
get_unique_tissues <- function(metadata) {
    tissues_unique <- metadata %>%
        group_by(tissue, disease) %>%
        summarize(n = n()) %>%
        filter(n > 5) %>%
        ungroup() %>%
        mutate(tissue_disease = paste(tissue, disease, sep = "_")) %>%
        pull(tissue_disease)

    return(tissues_unique)
}

## Function to calculate weight factor
calculate_weight <- function(tissues_nondisease_unique) {
    weight <- 1 / (length(tissues_nondisease_unique) - 1)
    return(weight)
}

## Function to calculate pd estimate
calculate_pd_estimate <- function(countdata_d, countdata_k, weight) {
    countdata_d$I_kidj <- apply(countdata_d[, 1:(ncol(countdata_d) - 3)],
                                1,
                                function(x) sum(countdata_k[, colnames(countdata_d[, 1:(ncol(countdata_d) - 3)])] < x))
    countdata_d$I_kidj_norm <- countdata_d$I_kidj / (ncol(countdata_k) * ncol(countdata_d[, 1:(ncol(countdata_d) - 3)]))
    countdata_d$ptot_dis <- countdata_d$ptot_dis + countdata_d$I_kidj_norm * weight
    return(countdata_d$ptot_dis)
}

## Function to calculate SPECS score
calculate_specs <- function(datafile, metafile) {
    ## Read in count data and metadata
    countdata <- read_delim(datafile, delim = "\t") %>%
        rename(GeneID = `""`) %>%
        select(-`""`) %>%
        as_tibble()

    metadata <- read_delim(metafile, delim = "\t") %>%
        as_tibble()

    ## Calculate TPM values
    tpm <- calculate_tpm(countdata)

    ## Determine unique tissues, with and without disease
    metadata_all <- metadata %>%
        filter(sample_name %in% colnames(tpm)) %>%
        filter(samp_qc_score %>% as.numeric() < 4)

    metadata_nondisease <- metadata_all %>%
        filter(disease == "0")

    metadata_disease <- metadata_all %>%
        filter(disease == "1")

    tissues_nondisease_unique <- get_unique_tissues(metadata_nondisease)

     ## Calculate weight factor
     weight <- calculate_weight(tissues_nondisease_unique)

     nsamps <- length(tissues_nondisease_unique)
     nrna <- nrow(tpm)

     ## Read names transcripts
     mw_RNA <- matrix(nrow = nrna, ncol = length(tissues_nondisease_unique))

     ## Run through next loop for each tissue type (d)
     for (tissue in tissues_nondisease_unique) {
         tissues_noncurrent <- setdiff(tissues_nondisease_unique, tissue)
         tiss <- str_split(tissue, "_", simplify = TRUE)[1]
         current_samples <- metadata_all %>%
             filter(tissue == tiss) %>%
             filter(disease == "0") %>%
             pull(sample_name)
         countdata_d <- tpm %>%
             select(current_samples) %>%
             as.matrix()
         colnames(countdata_d) <- current_samples
         countdata_d <- cbind(countdata_d,
                              I_kidj = rep(0, nrow(countdata_d)),
                              ptot_dis = rep(0, nrow(countdata_d)),
                              I_kidj_norm = rep(0, nrow(countdata_d)))
                                        # Compare with each tissue type not d (k)
         for (comp_type in tissues_noncurrent) {
             comp_tiss <- str_split(comp_type, "_", simplify = TRUE)[1]
             comp_samples <- metadata_all %>%
                 filter(tissue == comp_tiss) %>%
                 filter(disease == "0") %>%
                 pull(sample_name)
             countdata_k <- tpm %>%
                 select(comp_samples) %>%
                 as.matrix()
             colnames(countdata_k) <- comp_samples
             countdata_d$ptot_dis <- calculate_pd_estimate(countdata_d, countdata_k, weight)
         }
         mw_RNA[, which(tissues_nondisease_unique == tissue)] <- countdata_d$ptot_dis
     }

     colnames(mw_RNA) <- tissues_nondisease_unique

     ## Write out data
     write_tsv(
         mw_RNA,
         file.path(dirname(metafile), "results", "nondisease_genes_filt_qc123_specs_pall.txt")
     )

     mw_max_type <- data.frame(
         maxval = apply(mw_RNA, 1, max),
         maxtype = apply(mw_RNA, 1, function(x) colnames(mw_RNA)[which(x == max(x))])
     )
     mw_max_type <- mw_max_type[order(-mw_max_type$maxval, mw_max_type$maxtype), ]
     write_tsv(
         mw_max_type,
         file.path(dirname(metafile), "results", "nondisease_genes_filt_qc123_specs_maxval.txt"),
         row.names = FALSE,
         col.names = FALSE
     )

     mw_min_type <- data.frame(
         minval = apply(mw_RNA, 1, min),
         mintype = apply(mw_RNA, 1, function(x) colnames(mw_RNA)[which(x == min(x))])
     )
     mw_min_type <- mw_min_type[order(-mw_min_type$minval, mw_min_type$mintype), ]
     write_tsv(
         mw_min_type,
         file.path(dirname(metafile), "results", "nondisease_genes_filt_qc123_specs_min.txt"),
         row.names = FALSE,
         col.names = FALSE
     )
}

######################################################################
### SPECS.R ends here
