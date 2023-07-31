# Code to accompany the manuscript "Deconvolution of Nascent Sequencing Data Using Transcriptional Regulatory Elements"

This repository contains code to accompany the manuscript "Nascent Sequencing Data Facilitates Highly Accurate Deconvolution of Mixtures of Fully Differentiated Celltypes". All code for this project can be found in the `/src` directory and a running lab notebook from during the project is located in the `notes.md` file, documenting the process of developing this work.

Data is available in the [./data](data) while supplementary documentation for this work are provided in the [./doc](doc) directory.

## Project Dependencies

This work was developed in R, in addition to the use of published pipelines (more information on these can be found in the [./doc](doc) directory). The following R packages were used in this project:
- tidyverse
- ggthemes
- e1071
- LiblineaR
- MASS
- nnls
- progress
- corrplot
- rstan
- glmnet
- expm
- Matrix
- expm
- memoise
- gtools
- scales
- patchwork

Shared implementation details [./src/decon_utils.R](src/decon_utils.R) and which implements a standard interface used for other scripts in this project.

Scripts used in the manuscript
- `bidir_overlap_stats.R` - Generates overlap statistics for bidirectional calls
- `decon_utils.R` - Shared implementation details
- `gen_gene_counts.sbatch` - Generates gene counts from bam files
- `gen_subset.sbatch` - Generates titration trials over bidirectionals
- `gen_subset_random.sbatch` - Generates random mixtures from bidirectionals
- `gen_subset_genes.sbatch` - Generates random mixtures from genes
- `run_bidirs.r` - Runs deconvolution over bidirectionals, plus other early exploratory analysis
- `run_genes.r` - Runs deconvolution over genes
- `run_merged.r` - Runs deconvolution over merged bidirectionals and genes
- `submit_all_subset.sh` - Submits titration subset jobs
- `submit_gene_subset.sh` - Submits gene subset jobs
- `submit_random_subset.sh` - Submits random subset jobs

Scripts from earlier, unused analysis:
- `SPECS.R`
- `analyze_collated_results.R`
- `check_readcount_simulations.R`
- `collate_iteration_results.py`
- `gen_giggle_db.sbatch`
- `gen_giggle_table.py`
- `giggle_parse.awk`
- `install_deps.R`
- `linreg.rds`
- `linreg.stan`
- `mcmc_decon.R`
- `robust.stan`
- `run_svm_iter.R`
- `search_giggle_db.sbatch`
- `parse_giggle_res.py`
- `specs_from_scratch.R`
- `specs_gene_filt_qc123.py`
- `test_iter_run.sh`
- `trisomy_test.R`
- `run_iters.sh`
