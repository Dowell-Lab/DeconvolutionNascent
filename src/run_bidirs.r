### run_svm.r --- Run the nu-svm model
##
## Filename: run_svm.r
## Author: Zach Maas
## Created: Tue Jul 12 16:37:23 2022 (-0600)
##
######################################################################
##
### Commentary:
##
## This is a main script to run the nu-svm model used by this
## deconvolution algorithm. Eventually this file may also contain
## prefiltering steps as well. We should optimize the objective using
## the L-Curve method for SVM
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

## Simulate data for 6 samples
## num_samples <- 6
## size <- 100
## X <- matrix(runif(n=num_samples*size,min=0,max=10),ncol=num_samples,nrow=size)
## Y <- rowSums(X*(1/num_samples))

## ## Generate a handful of test nu values and get the best value
## ## test_nus <- c(10^(-7:-1), 10^-1*(2:9), 1 - 10^(-1:-7), 1)
## test_nus <- c(1:1000*10^-3)
## best_nu <- get_best_nu(X,Y,test_nus)

## ## Run the model and estimate coefficients
## model <- nu_svm(X, Y, best_nu)
## coefs <- coef(model)[-1] ## Remove the intercept term
## coefs <- ifelse(coefs > 0, coefs, 0) ## Force terms to be positive
## coefs <- coefs / sum(coefs) ## Force terms to sum to one
## coef_rms <- sqrt(sum((coefs - (1/6))^2))

## Import some real data to test on
print("Importing Data")
data_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/"
samples_merged <- read_merged_samples(data_dir)
## samples_merged <- samples_merged %>%
## filter(rowSums(across(where(is.numeric)) != 0) > 6) %>%
## filter(rowSums(across(where(is.numeric))) > 10)

## Generate filter set for later example
filter_regions <- rowSums(samples_merged[2:10] != 0) > 1
samples_merged <- samples_merged[filter_regions, ]

## Read in and process enhancer/promoter data
roi_info <- read_delim(paste0(data_dir, "rois_with_intersect_count.txt"),
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

## Generate distribution graphs
ggplot(data = samples_merged) +
    geom_density(aes(x = hct, color="HCT")) +
    geom_density(aes(x = hela, color="HeLa")) +
    geom_density(aes(x = mcf7, color="MCF7")) +
    geom_density(aes(x = k562, color="K562")) +
    geom_density(aes(x = esc, color="ESC")) +
    geom_density(aes(x = kasumi, color="Kasumi")) +
    geom_density(aes(x = cd4, color="CD4")) +
    geom_density(aes(x = jurkat, color="Jurkat")) +
    geom_density(aes(x = bj5ta, color="BJ5TA")) +
    geom_density(aes(x = lcl, color="LCL")) +
    scale_x_log10()
ggsave("distr.pdf", width=8, height=8)


## These stay consistent
counts_homo <- as.matrix(samples_merged[2:ncol(samples_merged)])
rownames(counts_homo) <- samples_merged$GeneID

print("Running Titrations")
## These can vary
num_points <- 100
## mixing_fractions <- c(1:num_points*(1 / num_points))
mixing_fractions <- seq(0, 1, length.out=num_points)
res_tbl <- gen_res_tbl()
pb <- progress_bar$new(format = "Models [:bar] :percent eta: :eta",
                       total = num_points)
for (mixing_fraction in mixing_fractions) {
    coefs_true <- seq(mixing_fraction, 1, length.out=9)
    coefs_true <- coefs_true / sum(coefs_true)
    counts_hetero <- as.integer(rowSums(sweep(counts_homo, MARGIN=2, coefs_true, `*`)))
    ## Add lots of noise
    counts_hetero <- rpois(length(counts_hetero), lambda=counts_hetero)
    ## Scale things
    ## counts_hetero <- scale(counts_hetero)
    ## counts_homo <- scale(counts_homo)
    ## Matrix Inversion
    ## coefs <- invert_matrix(counts_homo, counts_hetero)
    ## NNLS
    coefs <- decon_nnls(counts_homo, counts_hetero)
    ## coefs <- decon_lasso(counts_homo, counts_hetero)
    ## Epsilon SVR
    ## coefs <- eps_svm(counts_homo, counts_hetero, 0.25)
    ## Nu SVR
    ## coefs <- nu_svm(counts_homo, counts_hetero, 0.25)
    ## MCMC
    ## coefs <- decon_mcmc(counts_homo, counts_hetero)
    ## Calculate rms
    rms <- get_rms(coefs, coefs_true)
    ## res_tbl <- res_tbl %>%
    ##     add_row(c(mixing_fraction, rms, coefs_true, coefs))
    run_data <- as_tibble(t(c(mixing_fraction, rms,
                              coefs_true, coefs)))
    colnames(run_data) <- colnames(res_tbl)
    res_tbl <- bind_rows(res_tbl, run_data)
    pb$tick()
}
pb$terminate()
ggplot(data = res_tbl) +
    geom_line(aes(x = mixing_fraction, y = hct_actual      , linetype = "actual", color = "HCT116")) +
    geom_line(aes(x = mixing_fraction, y = hela_actual     , linetype = "actual", color = "HeLa")) +
    geom_line(aes(x = mixing_fraction, y = mcf7_actual     , linetype = "actual", color = "MCF7")) +
    geom_line(aes(x = mixing_fraction, y = k562_actual     , linetype = "actual", color = "K562")) +
    geom_line(aes(x = mixing_fraction, y = esc_actual      , linetype = "actual", color = "ESC")) +
    geom_line(aes(x = mixing_fraction, y = kasumi_actual   , linetype = "actual", color = "Kasumi")) +
    geom_line(aes(x = mixing_fraction, y = cd4_actual      , linetype = "actual", color = "CD4")) +
    geom_line(aes(x = mixing_fraction, y = jurkat_actual   , linetype = "actual", color = "Jurkat")) +
    geom_line(aes(x = mixing_fraction, y = bj5ta_actual    , linetype = "actual", color = "BJ5TA")) +
    geom_line(aes(x = mixing_fraction, y = hct_estimated   , linetype = "estimated", color = "HCT116")) +
    geom_line(aes(x = mixing_fraction, y = hela_estimated  , linetype = "estimated", color = "HeLa")) +
    geom_line(aes(x = mixing_fraction, y = mcf7_estimated  , linetype = "estimated", color = "MCF7")) +
    geom_line(aes(x = mixing_fraction, y = k562_estimated  , linetype = "estimated", color = "K562")) +
    geom_line(aes(x = mixing_fraction, y = esc_estimated   , linetype = "estimated", color = "ESC")) +
    geom_line(aes(x = mixing_fraction, y = kasumi_estimated, linetype = "estimated", color = "Kasumi")) +
    geom_line(aes(x = mixing_fraction, y = cd4_estimated   , linetype = "estimated", color = "CD4")) +
    geom_line(aes(x = mixing_fraction, y = jurkat_estimated, linetype = "estimated", color = "Jurkat")) +
    geom_line(aes(x = mixing_fraction, y = bj5ta_estimated , linetype = "estimated", color = "BJ5TA")) +
    labs(title = "NNLS Deconvolution Accuracy with Poisson Noise",
         x = "Progress to Equal Mixing",
         y = "Estimated Mixing Fraction",
         color = "Sample",
         linetype = "Condition") +
    scale_linetype_manual(values=c("dashed","solid"))
ggsave("mixing_accuracy_naive.pdf", width=8, height=8)

ggplot(data = res_tbl) +
    geom_point(aes(x = mixing_fraction, y = rms))

print("Importing subset data")
## Now do it on data from random subsampling of real data

## Automate this?
variation_dirs <- c("/with_all",
                    "/without_esc_bj5ta",
                    "/without_esc",
                    "/with_esc",
                    "/add_lcl")
subset_excl <- list("/with_all"=c("bj5ta", "lcl"),
                    "/add_lcl"=c("esc", "bj5ta"),
                    "/without_esc"=c("esc", "lcl"),
                    "/with_esc"=c("bj5ta", "lcl"),
                    "/without_esc_bj5ta"=c("esc", "bj5ta", "lcl"))
for (sub_dir in variation_dirs) {
    print(sub_dir)
    root_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/simulated_counts"
    data_dir <- paste0(root_dir, sub_dir, "/merged_")

    samples_sub <- read_sample(1, data_dir)
    for (sample_iter in seq(2,99)) {
        sample_sub <- read_sample(sample_iter, data_dir)
        samples_sub <- inner_join(samples_sub, sample_sub, by="Geneid")
    }
    ## Filter as before
    samples_sub <- samples_sub[filter_regions, ]

    print("Running subset simulation")
    all_rms <- c()
    ## for (num_regions in round(10^seq(1,5,length.out=100))) {
    num_regions <- 1000
    print(paste0("n=",num_regions))
    subset_idxs <- find_top_genes_per_celltype(counts_homo, num_regions)
    ## mixing_fractions <- c(1:num_points*(1 / num_points))
    res_tbl_2 <- gen_res_tbl()
    pb <- progress_bar$new(format = "Models [:bar] :percent eta: :eta",
                           total = num_points)
    i <- 0
    for (iter in seq(2,99)) {
        i <- i + 1
        mixing_fraction <- mixing_fractions[iter]
        if (sub_dir == "/add_lcl") {
            coefs_true <- seq((1/99)*i, 1, length.out=10) /
                sum(seq((1/99)*i, 1, length.out=10))
            ## coefs_true <- seq(mixing_fraction, 1, length.out=10)
        } else {
            coefs_true <- seq((1/99)*i, 1, length.out=9) /
                sum(seq((1/99)*i, 1, length.out=9))
            ## coefs_true <- seq(mixing_fraction, 1, length.out=9)
            coefs_true <- c(coefs_true, 0)
        }
        names(coefs_true) <- colnames(counts_homo)
        coefs_true[names(coefs_true) %in% unlist(subset_excl[sub_dir])] <- 0
        ## coefs_true[5] <- 0
        ## coefs_true[9] <- 0
        coefs_true <- coefs_true / sum(coefs_true)
        counts_hetero <- as.matrix(samples_sub[iter])
        ## Subset for speed
        rs_counts_homo <- counts_homo[subset_idxs,]
        rs_counts_hetero <- counts_hetero[subset_idxs,]
        ## Rescale
        ## counts_hetero <- scale(counts_hetero, center = TRUE, scale = TRUE)
        ## counts_homo <- scale(counts_homo, center = TRUE, scale = TRUE)
        coefs <- estimateMixing(rs_counts_homo, rs_counts_hetero, model_choice="nnls")
        rms <- get_rms(coefs, coefs_true)
        run_data <- as_tibble(t(c(mixing_fraction, rms,
                                  coefs_true, coefs)))
        colnames(run_data) <- colnames(res_tbl_2)
        res_tbl_2 <- bind_rows(res_tbl_2, run_data)
        pb$tick()
    }
    pb$terminate()

    ggplot(data = res_tbl_2) +
        geom_line(aes(x = mixing_fraction, y = hct_actual      , linetype = "actual", color = "HCT116")) +
        geom_line(aes(x = mixing_fraction, y = hela_actual     , linetype = "actual", color = "HeLa")) +
        geom_line(aes(x = mixing_fraction, y = mcf7_actual     , linetype = "actual", color = "MCF7")) +
        geom_line(aes(x = mixing_fraction, y = k562_actual     , linetype = "actual", color = "K562")) +
        geom_line(aes(x = mixing_fraction, y = esc_actual      , linetype = "actual", color = "ESC")) +
        geom_line(aes(x = mixing_fraction, y = kasumi_actual   , linetype = "actual", color = "Kasumi")) +
        geom_line(aes(x = mixing_fraction, y = cd4_actual      , linetype = "actual", color = "CD4")) +
        geom_line(aes(x = mixing_fraction, y = jurkat_actual   , linetype = "actual", color = "Jurkat")) +
        geom_line(aes(x = mixing_fraction, y = bj5ta_actual    , linetype = "actual", color = "BJ5TA")) +
        geom_line(aes(x = mixing_fraction, y = lcl_actual      , linetype = "actual", color = "LCL")) +
        geom_line(aes(x = mixing_fraction, y = hct_estimated   , linetype = "estimated", color = "HCT116")) +
        geom_line(aes(x = mixing_fraction, y = hela_estimated  , linetype = "estimated", color = "HeLa")) +
        geom_line(aes(x = mixing_fraction, y = mcf7_estimated  , linetype = "estimated", color = "MCF7")) +
        geom_line(aes(x = mixing_fraction, y = k562_estimated  , linetype = "estimated", color = "K562")) +
        geom_line(aes(x = mixing_fraction, y = esc_estimated   , linetype = "estimated", color = "ESC")) +
        geom_line(aes(x = mixing_fraction, y = kasumi_estimated, linetype = "estimated", color = "Kasumi")) +
        geom_line(aes(x = mixing_fraction, y = cd4_estimated   , linetype = "estimated", color = "CD4")) +
        geom_line(aes(x = mixing_fraction, y = jurkat_estimated, linetype = "estimated", color = "Jurkat")) +
        geom_line(aes(x = mixing_fraction, y = bj5ta_estimated , linetype = "estimated", color = "BJ5TA")) +
        geom_line(aes(x = mixing_fraction, y = lcl_estimated   , linetype = "estimated", color = "LCL")) +
        labs(title = "",
             ## paste0("Deconvolution with n=", num_regions,
             ##        " ROIs per celltype, RMS=", sum(res_tbl_2$rms)),
             x = "Progress to Equal Mixing",
             y = "Estimated Mixing Fraction",
             color = "Sample",
             linetype = "Condition") +
        scale_linetype_manual(values=c("dashed","solid")) +
        theme_minimal() +
        theme(text=element_text(family="Helvetica", size=18))
    proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
    mix_dir <- paste0(proj_dir, "figures/differentiation_tests")
    ggsave(paste0(mix_dir, sub_dir, ".pdf"), width=8, height=8)
    ## all_rms <- c(all_rms, sum(res_tbl_2$rms))
    ## }
}

## Random linearity tests to match literature
## Load in data
data_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/random_trials/merged_"
samples_random <- read_sample(1, data_dir)
for (sample_iter in seq(2,128)) {
    sample_random <- read_sample(sample_iter, data_dir)
    samples_random <- inner_join(samples_random, sample_random, by="Geneid")
}
## Filter as before
samples_random <- samples_random[filter_regions, ]

## Load in true coefficients
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
)
for (iter in seq_len(128)) {
    new_coefs <- read_delim(paste0("data/random_trials/coefs_",
                                   iter, ".txt"),
                            col_names=c("hct","hela","mcf7",
                                        "k562","esc","kasumi",
                                        "cd4","jurkat","bj5ta"),
                            show_col_types = FALSE) %>%
        mutate(trial=iter)
    coef_df <- coef_df %>% add_row(new_coefs)
}

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
            coefs_true <- unlist((coef_df %>% subset(trial == iter))[2:10])
            ## coefs_true[5] <- 0
            ## coefs_true[9] <- 0
            coefs_true <- coefs_true[c(1,2,3,4,6,7,8)]
            coefs_true <- coefs_true / sum(coefs_true)
            counts_hetero <- as.matrix(samples_random[,iter+1])
            ## Subset for speed
            rs_counts_homo <- counts_homo_usedcells[subset_idxs,]
            rs_counts_hetero <- counts_hetero[subset_idxs,]
            ## Do Estimation
            coefs <- estimateMixing(rs_counts_homo, rs_counts_hetero,
                                    model_choice=method, param=0.2)
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

        ## Plot each of these and arrange together
        do_plot <- FALSE
        if (do_plot) {
            proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
            linear_dir <- paste0(proj_dir, "figures/linearity_tests/")
            celltypes <- c("hct", "hela", "mcf7",
                           "k562", "esc", "kasumi",
                           "cd4", "jurkat", "bj5ta")
            for (celltype in celltypes) {
                ggplot(data=res_tbl_3) +
                    geom_point(aes(x=.data[[paste0(celltype, "_actual")]],
                                   y=.data[[paste0(celltype, "_estimated")]])) +
                    geom_abline(slope=1,intercept=0,color="red") +
                    coord_fixed() +
                    theme_minimal()
                ggsave(paste0(linear_dir, "linearity_", celltype,".pdf"), width=8, height=8)
            }
        }
    }
    rms_df[paste0(method, "_", method)] <- all_rms
}
save(multi_res_tbl, file="regions_bidirs_res.Rdata")
## save(multi_res_tbl, file="multi_res_percentile_celltypes.Rdata")

## Old multi_res_tbl without percentile filtering
## load('multi_res.Rdata')

## Add enhancer/promoter information (done here since the long-slow
## run is already done and I don't want to waste hours on compute)
roi_df <- tibble(num_regions=integer(),
                 enhancer_count=integer(),
                 promoter_count=integer())
regions <- unique(multi_res_tbl$num_regions)
## Testing doing a denser grid for plotting the E/P ratio
regions <- round(10^seq(1,5,length.out=1000))
for (num_regions in regions) {
    subset_idxs <- find_top_genes_per_celltype(counts_homo, num_regions)
    roi_types <- roi_info[subset_idxs,] %>% count(region_type)
    roi_df <- add_row(roi_df, tibble(num_regions=c(num_regions),
                                     enhancer_count=c(c(roi_types[1,2])$n),
                                     promoter_count=c(c(roi_types[2,2])$n)))
}
roi_df <- mutate_all(roi_df, ~coalesce(.,0))

## multi_res_tbl <- multi_res_tbl %>% left_join(roi_df, by="num_regions")
region_counts <- roi_info %>%
    group_by(region_type) %>%
    summarise(count=length(region_type))
total_enhancer <- (region_counts %>% filter(region_type == "enhancer"))$count
total_promoter <- (region_counts %>% filter(region_type == "promoter"))$count

roi_df<- roi_df %>%
    mutate(roi_count = enhancer_count+promoter_count,
           promoter_frac = promoter_count / roi_count,
           enhancer_frac = enhancer_count / roi_count,
           roi_ratio = enhancer_count / promoter_count,
           pval = phyper(enhancer_count, total_enhancer,
                         total_promoter, roi_count,
                         lower.tail = TRUE),
           padj = pmin(pval * nrow(unique(roi_df)), 1),
           significant = padj < 0.05)

##TODO Hypergeometric test
## pval = phyper(enhancer_count, total_enhancer, total_promoter, roi_count)
## padj = pval / nrow(unique(roi_df))

proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
roi_dir <- paste0(proj_dir, "figures/roi_tests/")

roi_signif <- ggplot(roi_df) +
    geom_line(aes(x=num_regions,y=-log10(padj))) +
    geom_hline(yintercept=-log10(0.05), color="red", linetype="dotted") +
    scale_x_log10() +
    theme_minimal() +
    labs(x="Number of Regions Per Celltype",
         y="-log10(p-adjusted)",
         title="Enhancer Enrichment (Hypergeometric Test)")

ggplot(data = roi_df) +
    geom_line(aes(x=num_regions, y=enhancer_frac, color="Enhancer")) +
    geom_line(aes(x=num_regions, y=promoter_frac, color="Promoter")) +
    scale_x_log10() +
    theme_bw() +
    labs(x="Number of Regions Per-Sample in Subset",
         y="Proportion of ROIs in Group",
         title="Enhancers Dominate Sample-Determining ROIs",
         color="Region Type") +
    theme_minimal() +
    theme(text=element_text(family="Helvetica", size=16))
ggsave(paste0(roi_dir, "ROI_Fraction.pdf"), width=8, height=8)

ggplot(data = roi_df) +
    geom_line(aes(x=num_regions, y=roi_ratio)) +
    scale_x_log10() +
    theme_minimal() +
    lims(y=c(0,max(roi_df$roi_ratio[!is.infinite(roi_df$roi_ratio)]))) +
    labs(x="Number of Regions Per-Sample in Subset",
         y="Enhancer:Promoter Ratio of ROIs",
         title="ROI Selection Picks Many-Fold More Enhancers than Promoters") +
    theme_bw() +
    theme(text=element_text(family="Helvetica"))
ggsave(paste0(roi_dir, "ROI_Ratio.pdf"), width=8, height=8)

## Plot linear group plots for each number of regions
agg_num <- multi_res_tbl %>%
    group_by(num_regions) %>%
    nest()

for (i in seq_len(100)) {
    iter_df <- unnest(agg_num[i,2], cols=c(data))
    proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
    linear_dir <- paste0(proj_dir, "figures/linearity_tests/")
    celltypes <- c("hct", "hela", "mcf7",
                   "k562", "esc", "kasumi",
                   "cd4", "jurkat", "bj5ta")
    for (celltype in celltypes) {
        ggplot(data=iter_df) +
            geom_point(aes(x=.data[[paste0(celltype, "_actual")]],
                           y=.data[[paste0(celltype, "_estimated")]],
                           color=method)) +
            geom_smooth(method="lm",
                        aes(x=.data[[paste0(celltype, "_actual")]],
                            y=.data[[paste0(celltype, "_estimated")]],
                            color=method),
                        formula = 'y~x',
                        size=0.5,
                        alpha=0.5,
                        se=FALSE) +
            labs("Actual Proportion", "Estimated Proportion",
                 title=paste0("Linearity of ", celltype, " Estimates, N=",
                              c(agg_num[i,1])$num_regions),
                 color="Algorithm") +
            geom_abline(slope=1,intercept=0,color="red",linetype="dotted") +
            scale_color_discrete(name="Algorithm",
                                 labels=c("SVR",
                                          "LASSO",
                                          "NNLS",
                                          "Ridge")) +
            ## coord_fixed() +
            xlim(0,0.5) +
            ylim(0,0.5) +
            theme_minimal()
        ggsave(paste0(linear_dir, "/roi_filter/linearity_", celltype, "_roi:",
                      c(agg_num[i,1])$num_regions, ".pdf"), width=8, height=8)
    }
}

## Plot aggregated RMS plot
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
linear_dir <- paste0(proj_dir, "figures/linearity_tests/")
ggsave(paste0(linear_dir, "aggregate_error_percentile.pdf"), width=8, height=8)

## Subsetted version
sub_agg_rms <- sub_multi_res_tbl %>% group_by(method, num_regions) %>% summarise(all_rms = sum(rms))
aggregate_plot <- ggplot(data = sub_agg_rms) +
    geom_line(aes(x=num_regions, y=all_rms, color=method)) +
    scale_x_log10() +
    theme_minimal() +
    labs(title="System Error vs Number of Regions",
         x="",
         y="Total RMS Error Over 128 Runs",
         color="Algorithm")
proj_dir <- "/Users/zachmaas/Dropbox/phd/research/dna_lab/deconascent/"
linear_dir <- paste0(proj_dir, "figures/linearity_tests/")
ggsave(paste0(linear_dir, "aggregate_error_sub.pdf"), width=8, height=8)

## Plot padj plot with aggregate RMS
roi_signif / aggregate_plot  +
    plot_layout(widths=c(1,1), heights=c(1,0.5)) +
    plot_annotation(tag_levels="A") +
    labs(x="Number of Regions Per Celltype")
ggsave(paste0(linear_dir, "aggregate_error_merged.pdf"), width=8, height=8)

## Testing quality of selecting top genes
topgene_simulations <- tibble(num_genes=round(10^seq(1,5,length.out=100)),
                              rms = all_rms)
ggplot(data = topgene_simulations) + geom_line(aes(x=num_genes, y=rms)) +
    labs(title="NNLS Deconvolution Subsampling Top Genes",
         x="Genes selected per sample", y="Total RMS across 100 step titration") +
    theme_minimal() +
    scale_x_log10()
ggsave(paste0(mix_dir, "topgene_accuracy.pdf"), width=8, height=8)
ggsave(paste0(mix_dir, "topgene_accuracy.png"), width=8, height=8)

data_cor <- cor(samples_merged[2:10])
pdf(file="correlation.pdf")
corrplot(data_cor, order = "hclust",
         method="color", tl.col = "black", tl.srt = 45)
dev.off()

## Number of nonzero columns
colSums(samples_merged[2:10] != 0) / nrow(samples_merged)

## save(counts_homo, samples_sub, mixing_fractions,
##      file="test_data.Rdata", compress='xz')

without_esc <- counts_homo[, colnames(counts_homo)!=c("esc")]
without_bj5ta <- counts_homo[, colnames(counts_homo)!=c("bj5ta")]
differentiated <- without_bj5ta[, colnames(without_bj5ta)!=c("esc")]
print("Condition Numbers:")
kappa(differentiated)
kappa(without_esc)
kappa(without_bj5ta)
kappa(counts_homo)
print("Matrix Rank (approximates information content):")
rankMatrix(differentiated)
rankMatrix(without_esc)
rankMatrix(without_bj5ta)
rankMatrix(counts_homo)
## print("Von Neumann Information Content")
## vonNeumann(differentiated)
## vonNeumann(without_esc)
## vonNeumann(without_bj5ta)
## vonNeumann(counts_homo)

pca <- prcomp(t(counts_homo), scale=TRUE)

## Attempt to generate a pseudo-ESC celltype without specific markers
homo_binary <- as.matrix((counts_homo>0)+0)
esc_pure <- homo_binary[,5]
differentiated <- homo_binary[, colnames(homo_binary)!="esc"]
shared_markers <- apply(differentiated, 1, function(x) { Reduce("|", x)} )
## for (celltype in colnames(homo_binary)) {

## }

library('topicmodels')
topics <- LDA(t(counts_homo), k=9)
print(tidy(topics, matrix='beta'), n=300)

top_rois <- tidy(topics) %>%
    group_by(topic) %>%
    slice_max(beta, n = 10) %>%
    ungroup() %>%
    arrange(topic, -beta)

beta_wide <- tidy(topics) %>%
    mutate(topic = paste0("topic", topic)) %>%
    pivot_wider(names_from = topic, values_from = beta) %>%
    mutate(log_ratio = log2(topic2 / topic1)) %>%
    arrange(desc(log_ratio))


for (celltype in colnames(counts_homo)) {
    print(sum(counts_homo[,colnames(counts_homo)==celltype] != 0))
}

subset_idxs <- find_top_genes_per_celltype(counts_homo, 1000)
homo_filtered <- counts_homo[subset_idxs,]
hetero_filtered <- counts_hetero[subset_idxs,]

info_content <- cor(counts_homo)
for (i in seq_len(10)) {
    for (j in seq_len(10)) {
        ## info_content[i,j] <- round(mutinformation(counts_homo[,i], counts_homo[,j]), 3)
        ## info_content[i,j] <- round(jsd(counts_homo[,i],
        ##                               counts_homo[,j]), 3)
        info_content[i,j] <- round(mean(jsd(counts_homo[,i],
                                            counts_homo[,j])), 3)
        ## if (i == j) {
        ##     info_content[i,j] <- 1
        ## }
    }
}

jsd <- function(x,y) {
    if (length(x) != length(y)) {
        print("Incompatible matrix dimension")
    }
    norms <- matrix(0,nrow=length(x),ncol=3)
    colnames(norms) <- c("x","y","m")
    norms[,"x"] <- t(x) / sum(x)
    norms[,"y"] <- t(y) / sum(y)
    norms[,"m"] <- 0.5*(norms[,"x"]+norms[,"y"])
    norm_sub <- apply(norms, 1, function(row) any(row != 0))
    norms <- norms[norm_sub,]
    probs <- matrix(0,nrow=nrow(norms),ncol=2)
    colnames(probs) <- c("logx","logy")
    probs[,"logx"] <- norms[,"x"] * log2(norms[,"x"]/norms[,"m"])
    probs[,"logy"] <- norms[,"y"] * log2(norms[,"y"]/norms[,"m"])
    probs_sub <- apply(probs, 1, function(row) all(!is.na(row)))
    probs <- probs[probs_sub,]
    divergence_x <- sum(probs[,"logx"])
    divergence_y <- sum(probs[,"logy"])
    divergence <- 0.5*divergence_x + 0.5*divergence_y
    return(divergence)
}

kl <- function(x,y) {
    norms <- matrix(0,nrow=length(x),ncol=3)
    colnames(norms) <- c("x","y","m")
    norms[,"x"] <- t(x) / sum(x)
    norms[,"y"] <- t(y) / sum(y)
    norm_sub <- apply(norms, 1, function(row) any(row != 0))
    norms <- norms[norm_sub,]
    probs <- matrix(0,nrow=nrow(norms),ncol=1)
    colnames(probs) <- c("logx")
    probs[,"logx"] <- norms[,"x"] * log2(norms[,"x"]/norms[,"y"])
    probs_sub <- apply(probs, 1, function(row) all(!is.na(row)))
    divergence <- sum(probs[,"logx"])
    return(divergence)
}

######################################################################
### run_svm.r ends here
