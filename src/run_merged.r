### run_merged.r ---
##
## Filename: run_merged.r
## Author: Zach Maas
## Created: Fri Jul 21 12:26:02 2023 (-0600)
##
######################################################################
##
### Commentary:
##
## Run consolidated analysis over both genes and enhancers. This is a
## combination of earlier analyses from the project.
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

print("Loading Utilities")
source('./src/decon_utils.R')

read_gene_samples <- function(data_dir) {
    sample_hct <- read_delim(paste0(data_dir   , "Jiang2018multi_counts.txt"),
                             show_col_types=FALSE, skip=1, delim="\t")
    sample_hela <- read_delim(paste0(data_dir  , "Fei2018ndf_counts.txt"),
                              show_col_types=FALSE, skip=1, delim="\t")
    sample_mcf7 <- read_delim(paste0(data_dir  , "Andrysik2017identification_counts.txt"),
                              show_col_types=FALSE, skip=1, delim="\t")
    sample_k562 <- read_delim(paste0(data_dir  , "Dukler2017nascent_counts.txt"),
                              show_col_types=FALSE, skip=1, delim="\t")
    sample_esc <- read_delim(paste0(data_dir   , "Smith2021peppro_counts.txt"),
                             show_col_types=FALSE, skip=1, delim="\t")
    sample_kasumi <- read_delim(paste0(data_dir, "Zhao2016high_counts.txt"),
                                show_col_types=FALSE, skip=1, delim="\t")
    sample_cd4 <- read_delim(paste0(data_dir   , "Danko2018dynamic_counts.txt"),
                             show_col_types=FALSE, skip=1, delim="\t")
    sample_jurkat <- read_delim(paste0(data_dir, "Chu2018chromatin_counts.txt"),
                                show_col_types=FALSE, skip=1, delim="\t")
    sample_bj5ta <- read_delim(paste0(data_dir , "Ikegami2020phosphorylated_counts.txt"),
                               show_col_types=FALSE, skip=1, delim="\t")
    sample_lcl <- read_delim(paste0(data_dir   , "Core2014analysis_counts.txt"),
                             show_col_types=FALSE, skip=1, delim="\t")
    genes_merged <- inner_join(sample_hct, sample_hela,
                                 by="Geneid") %>%
        inner_join(sample_mcf7, by="Geneid") %>%
        inner_join(sample_k562, by="Geneid") %>%
        inner_join(sample_esc, by="Geneid") %>%
        inner_join(sample_kasumi, by="Geneid") %>%
        inner_join(sample_cd4, by="Geneid") %>%
        inner_join(sample_jurkat, by="Geneid") %>%
        inner_join(sample_bj5ta, by="Geneid") %>%
        inner_join(sample_lcl, by="Geneid") %>%
        transmute(GeneID = Geneid,
                  hct    = `/scratch/Users/zama8258/deconvolution_gene_trials/Jiang2018multi.bam`,
                  hela   = `/scratch/Users/zama8258/deconvolution_gene_trials/Fei2018ndf.bam`,
                  mcf7   = `/scratch/Users/zama8258/deconvolution_gene_trials/Andrysik2017identification.bam`,
                  k562   = `/scratch/Users/zama8258/deconvolution_gene_trials/Dukler2017nascent.bam`,
                  esc    = `/scratch/Users/zama8258/deconvolution_gene_trials/Smith2021peppro.bam`,
                  kasumi = `/scratch/Users/zama8258/deconvolution_gene_trials/Zhao2016high.bam`,
                  cd4    = `/scratch/Users/zama8258/deconvolution_gene_trials/Danko2018dynamic.bam`,
                  jurkat = `/scratch/Users/zama8258/deconvolution_gene_trials/Chu2018chromatin.bam`,
                  bj5ta  = `/scratch/Users/zama8258/deconvolution_gene_trials/Ikegami2020phosphorylated.bam`,
                  lcl    = `/scratch/Users/zama8258/deconvolution_gene_trials/Core2014analysis.bam`) #%>%
    return(genes_merged)
}

read_merged_bidirs <- function(data_dir) {
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
    bidirs_merged <- inner_join(sample_hct, sample_hela,
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
    return(bidirs_merged)
}

## Read in bidirs
bidir_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/"
bidirs_merged <- read_merged_bidirs(bidir_dir)
bidirs_homo <- as.matrix(bidirs_merged[2:ncol(bidirs_merged)])
bidirs_homo <- bidirs_homo[,seq_len(9)]

## Filter to just enhancers
roi_info <- read_delim(paste0(bidir_dir, "rois_with_intersect_count.txt"),
                       delim="\t", col_names=c('GeneID', 'promoter_count')) %>%
    transmute(GeneID=GeneID,
              region_type = as_factor(ifelse(promoter_count > 0,
                                             "promoter",
                                             "enhancer"))) %>%
    separate_wider_regex(GeneID, c(chr=".*", ":",
                                   start=".*", "-",
                                   end=".*"),
                         cols_remove=FALSE) %>%
    mutate(start = as.numeric(start),
           end = as.numeric(end),
           length = end - start)

## Read in genes
gene_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/gene_counts/"
genes_merged <- read_gene_samples(gene_dir)
genes_homo <- as.matrix(genes_merged[2:ncol(genes_merged)])
rownames(genes_homo) <- genes_merged$GeneID
## genes_homo <- genes_homo[,c(1,2,3,4,6,7,8)]
genes_homo <- genes_homo[,seq_len(9)]

## Merge together
counts_homo <- rbind(genes_homo, bidirs_homo)
filter_regions <- rowSums(counts_homo != 0) > 1
counts_homo <- counts_homo[filter_regions, ]

## Get coefficients
coef_df <- tibble(
    trial  = integer(),
    hct    = double(),
    hela   = double(),
    mcf7   = double(),
    k562   = double(),
    esc    = double(),
    kasumi = double(),
    cd4    = double(),
    jurkat = double(),
    bj5ta  = double()
    ## lcl    = double()
)
for (iter in seq_len(128)) {
    new_coefs <- read_delim(paste0("data/gene_counts/coefs_",
                                   iter, ".txt"),
                            col_names=c("hct","hela","mcf7",
                                        "k562","esc","kasumi",
                                        "cd4","jurkat","bj5ta"),
                            show_col_types = FALSE) %>%
        mutate(trial=iter)
    new_coefs["esc"] = 0
    new_coefs["bj5ta"] = 0
    new_coefs[1:9] <- new_coefs[1:9] / sum(new_coefs[1:9])
    coef_df <- coef_df %>% add_row(new_coefs)
}
coef_df <- coef_df %>% subset(select=-c(esc, bj5ta))

## Load each heterogenous
## Genes
gene_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/gene_counts/merged_"
genes_random <- read_sample(1, gene_dir)
for (sample_iter in seq(2,128)) {
    sample_random <- read_sample(sample_iter, gene_dir)
    genes_random <- inner_join(genes_random, sample_random, by="Geneid")
}

## Bidirs
bidir_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/random_trials/merged_"
bidirs_random <- read_sample(1, bidir_dir)
for (bidir_iter in seq(2,128)) {
    bidir_random <- read_sample(bidir_iter, bidir_dir)
    bidirs_random <- inner_join(bidirs_random, bidir_random, by="Geneid")
}

merged_random <- rbind(genes_random, bidirs_random)
merged_random <- merged_random[filter_regions, ]

## ## Do deconvolution
## used_celltypes <- c(1,2,3,4,6,7,8)
## counts_homo_used <- counts_homo[,used_celltypes]
## multi_res_tbl <- gen_multires_tbl() %>%
##     subset(select=-c(lcl_estimated, lcl_actual,
##                      esc_estimated, esc_actual,
##                      bj5ta_actual, bj5ta_estimated))
## for (method in c("nnls", "ridge", "lasso", "eps_svm")) {
##     print(paste0("Processing ", method))
##     gene_res <- gen_multires_tbl() %>%
##         subset(select=-c(lcl_estimated, lcl_actual,
##                          esc_estimated, esc_actual,
##                          bj5ta_actual, bj5ta_estimated))
##     pb <- progress_bar$new(format = "Models [:bar] :percent eta: :eta",
##                            total = 128)
##     subset_idxs <- find_top_genes_per_celltype(counts_homo, 10000)
##     for (iter in seq_len(128)) {
##         ## Generate a random mix of ROIs
##         mixing_fraction <- NA
##         coefs_true <- unlist((coef_df %>% subset(trial == iter))[2:ncol(coef_df)])
##         coefs_true <- coefs_true / sum(coefs_true)
##         counts_hetero <- as.matrix(merged_random[,iter+1])
##         ## Subset for speed
##         rs_counts_homo_used <- counts_homo_used[subset_idxs,]
##         rs_counts_hetero <- counts_hetero[subset_idxs,]
##         ## Do Estimation
##         coefs <- estimateMixing(rs_counts_homo_used, rs_counts_hetero,
##                                 model_choice=method, param=0.2)
##         ## coefs <- estimateMixing(rs_counts_homo, rs_counts_hetero,
##         ##                         model_choice=method, param=0.2)
##         rms <- get_rms(coefs, coefs_true)
##         run_data <- as_tibble(t(c(mixing_fraction, rms,
##                                   NA, method,
##                                   coefs_true, coefs)))
##         run_data <- suppressMessages(type_convert(run_data))
##         colnames(run_data) <- colnames(gene_res)
##         gene_res <- bind_rows(gene_res, run_data)
##         pb$tick()
##     }
##     pb$terminate()
##     multi_res_tbl <- bind_rows(multi_res_tbl, gene_res)
## }

num_runs <- 100
rms_df <- tibble(iter=seq_len(num_runs),
                 num_regions=round(10^seq(1,5,length.out=num_runs)))
## Select only things included in this mixture
counts_homo_usedcells <- counts_homo[,c(1,2,3,4,6,7,8)]
multi_res_tbl <- gen_multires_tbl() %>%
    subset(select=-c(lcl_estimated, lcl_actual,
                     esc_estimated, esc_actual,
                     bj5ta_actual, bj5ta_estimated))
for (method in c("nnls", "ridge", "lasso", "eps_svm")) {
    print(paste0("Processing ", method))
    all_rms <- c()
    for (num_regions in rms_df$num_regions) {
        ## Run iterations
        run_res_tbl <- gen_multires_tbl() %>%
            subset(select=-c(lcl_estimated, lcl_actual,
                             esc_estimated, esc_actual,
                             bj5ta_actual, bj5ta_estimated))
        ## num_regions <- 1000
        pb <- progress_bar$new(format = "Models [:bar] :percent eta: :eta",
                               total = 128)
        subset_idxs <- find_top_genes_per_celltype(counts_homo_usedcells, num_regions)
        iter_rms <- c()
        for (iter in seq_len(128)) {
            ## Generate a random mix of ROIs
            mixing_fraction <- mixing_fractions[iter]
            coefs_true <- unlist((coef_df %>% subset(trial == iter))[seq(2,ncol(coef_df))])
            ## coefs_true[5] <- 0
            ## coefs_true[9] <- 0
            ## coefs_true <- coefs_true[c(1,2,3,4,6,7,8)]
            coefs_true <- coefs_true / sum(coefs_true)
            counts_hetero <- as.matrix(merged_random[,iter+1])
            ## Subset for speed
            rs_counts_homo <- counts_homo_usedcells[subset_idxs,]
            rs_counts_hetero <- counts_hetero[subset_idxs,]
            ## Do Estimation
            coefs <- estimateMixing(rs_counts_homo, rs_counts_hetero,
                                    model_choice=method, param=0.2)
            coefs[is.nan(coefs)] <- 0
            rms <- get_rms(coefs, coefs_true)
            iter_rms <- c(iter_rms, rms)
            run_data <- as_tibble(t(c(mixing_fraction, rms,
                                      num_regions, method,
                                      coefs_true, coefs)))
            run_data <- suppressMessages(type_convert(run_data))
            colnames(run_data) <- colnames(run_res_tbl)
            ## run_res_tbl <- bind_rows(run_res_tbl, run_data)
            run_res_tbl <- bind_rows(run_res_tbl, run_data)
            pb$tick()
        }
        pb$terminate()
        all_rms <- c(all_rms, sum(iter_rms))
        multi_res_tbl <- bind_rows(multi_res_tbl, run_res_tbl)
    }
    rms_df[paste0(method, "_", method)] <- all_rms
}
save(multi_res_tbl, file="regions_merged_res.Rdata")

print(multi_res_tbl %>% distinct() %>% group_by(method, num_regions) %>% summarise(agg_rms = sum(rms)), n=40)

## Plot the results
proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
linear_dir <- paste0(proj_dir, "figures/merged_tests/")
celltypes <- c("hct", "hela", "mcf7",
               "k562", "kasumi",
               "cd4", "jurkat")
for (celltype in celltypes) {
    ggplot(data=multi_res_tbl) +
        geom_point(aes(x=.data[[paste0(celltype, "_actual")]],
                       y=.data[[paste0(celltype, "_estimated")]],
                       color=method)) +
        geom_abline(slope=1,intercept=0,color="red") +
        coord_fixed() +
        scale_color_discrete(name="Algorithm",
                             labels=c("SVR",
                                      "LASSO",
                                      "NNLS",
                                      "Ridge")) +
        theme_minimal()
    ggsave(paste0(linear_dir, "genes_", celltype,".pdf"), width=8, height=8)
}

multi_res_tbl %>% group_by(method) %>% summarise(agg_rms = sum(rms))

agg_rms <- multi_res_tbl %>% group_by(method, num_regions) %>% summarise(all_rms = sum(rms))
aggregate_plot <- ggplot(data = agg_rms) +
    geom_line(aes(x=num_regions, y=all_rms, color=method)) +
    geom_hline(aes(yintercept=2.51, color="eps_svm"), linetype="dotted") +
    geom_hline(aes(yintercept=2.06, color="lasso"), linetype="dotted") +
    geom_hline(aes(yintercept=0.294, color="nnls"), linetype="dotted") +
    geom_hline(aes(yintercept=5.51, color="ridge"), linetype="dotted") +
    scale_x_log10() +
    theme_minimal() +
    scale_color_discrete(name="Algorithm",
                         labels=c("SVR",
                                  "LASSO",
                                  "NNLS",
                                  "Ridge")) +
    labs(title="System Error vs Number of Regions",
         x="",
         y="Total RMS Error Over 128 Runs",
         color="Algorithm")
proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
linear_dir <- paste0(proj_dir, "figures/merged_tests/")
ggsave(paste0(linear_dir, "aggregate_error_merged.pdf"), width=8, height=8)

######################################################################
### run_merged.r ends here
