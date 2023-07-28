#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH --partition=sandbox
#SBATCH --nodelist=fijinode-48
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --ntasks=64
#SBATCH --mem=500gb
#SBATCH --time=2:00:00
#SBATCH --mail-user=zama8258@colorado.edu

ulimit -n 12288
set -euo pipefail

module load bedtools
module load samtools

data_dir=/scratch/Users/zama8258/deconvolution_giggle_db
bgzip=/Users/zama8258/giggle/lib/htslib/bgzip
giggle=/Users/zama8258/giggle/bin/giggle
fd=~/usr/bin/fd
GIGGLE_ROOT=/Users/zama8258/giggle/
mkdir -p /scratch/Users/zama8258/deconvolution_giggle_db
"$fd" '.sorted_split_bidir_cov_filtered.bed' /Shares/dbnascent \
		| parallel --progress -j64 rsync {} \
							 /scratch/Users/zama8258/deconvolution_giggle_db
parallel -j64 --progress bedtools sort -i {} \
				 '|' "$bgzip" -c '>' {}.gz ::: "$data_dir"/*.bed
"$giggle" index -i "$data_dir/*.gz" -o "$data_dir"_b -f -s
