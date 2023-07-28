### trisomy_test.R --- Test simulated mosaicism
##
## Filename: trisomy_test.R
## Author: Zach Maas
## Created: Wed Nov 30 12:39:52 2022 (-0700)
##
######################################################################
##
### Commentary:
##
## This file contains code to test the deconvolution algorithm I've
## been using on mixtures of disomic and trisomic data
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

source('./src/decon_utils.R')

data_dir <- "~/Dropbox/phd/research/dna_lab/deconascent/data/trisomy/"
files <- c("Ethan37_rep1.sorted.unstranded.bidir_counts.txt",
           "Ethan37_rep2.sorted.unstranded.bidir_counts.txt",
           "Ethan42_rep1.sorted.unstranded.bidir_counts.txt",
           "Ethan42_rep2.sorted.unstranded.bidir_counts.txt",
           "Eric37_rep1.sorted.unstranded.bidir_counts.txt",
           "Eric37_rep2.sorted.unstranded.bidir_counts.txt",
           "Eric42_rep1.sorted.unstranded.bidir_counts.txt",
           "Eric42_rep2.sorted.unstranded.bidir_counts.txt")
samples_t21 <- read_delim(paste0(data_dir, files[1]), show_col_types=FALSE) %>%
    subset(select=-c(Source))
for (curr_file in files[2:8]) {
    sample_t21 <- read_delim(paste0(data_dir, curr_file), show_col_types=FALSE) %>%
        subset(select=-c(Source))
    samples_t21 <- samples_t21 %>% inner_join(sample_t21, by="GeneID")
}

counts_homo <- samples_t21 %>%
    subset(select=c("Ethan42.rep1.sorted.sorted.bam",
                    "Eric42.rep1.sorted.sorted.bam")) %>%
    as.matrix()

counts_homo_mix <- samples_t21 %>%
    subset(select=c("Ethan42.rep2.sorted.sorted.bam",
                    "Eric42.rep2.sorted.sorted.bam")) %>%
    as.matrix()

num_points <- 100
## mixing_fractions <- c(1:num_points*(1 / num_points))
mixing_fractions <- seq(0, 1, length.out=num_points)
res_tbl <- gen_res_tbl()
res_tbl <- tibble(
    mixing_fraction  = numeric(),
    rms              = numeric(),
    eric_actual       = numeric(),
    ethan_actual      = numeric(),
    eric_estimated = numeric(),
    ethan_estimated  = numeric()
)
pb <- progress_bar$new(format = "Models [:bar] :percent eta: :eta",
                       total = num_points)
for (mixing_fraction in mixing_fractions) {
    coefs_true <- seq(mixing_fraction, 1, length.out=2)
    coefs_true <- coefs_true / sum(coefs_true)
    counts_hetero <- as.integer(rowSums(sweep(counts_homo_mix, MARGIN=2, coefs_true, `*`)))
    ## Add lots of noise
    counts_hetero <- rpois(length(counts_hetero), lambda=counts_hetero)
    ## Scale things
    ## counts_hetero <- scale(counts_hetero)
    ## counts_homo <- scale(counts_homo)
    ## Matrix Inversion
    ## coefs <- invert_matrix(counts_homo, counts_hetero)
    ## NNLS
    coefs <- decon_nnls(counts_homo, counts_hetero)
    ## Epsilon SVR
    ## coefs <- eps_svm(counts_homo, counts_hetero, 0.25)
    ## Nu SVR
    ## coefs <- nu_svm(counts_homo, counts_hetero, 0.25)
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
    geom_line(aes(x = mixing_fraction, y = eric_actual      , linetype = "actual", color = "Eric")) +
    geom_line(aes(x = mixing_fraction, y = ethan_actual     , linetype = "actual", color = "Ethan")) +
    geom_line(aes(x = mixing_fraction, y = eric_estimated   , linetype = "estimated", color = "Eric")) +
    geom_line(aes(x = mixing_fraction, y = ethan_estimated  , linetype = "estimated", color = "Ethan")) +
    labs(title = "NNLS Deconvolution Accuracy",
         x = "Progress to Equal Mixing",
         y = "Estimated Mixing Fraction",
         color = "Sample",
         linetype = "Condition") +
    scale_linetype_manual(values=c("dashed","solid"))
ggsave("mixing_t21.pdf", width=8, height=8)

t21_cor <- cor(samples_t21[2:8])
pdf(file="cor_t21.pdf")
corrplot(t21_cor, order = "hclust",
         method="color", tl.col = "black", tl.srt = 45)
dev.off()

######################################################################
### trisomy_test.R ends here
