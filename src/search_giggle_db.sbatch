#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euxo pipefail

function finish {
		if command -v apprise &> /dev/null
		then
				apprise -t "Fiji Run Complete" -b "Giggle Search"
				exit
		fi
}
trap finish EXIT

data_dir=/scratch/Users/zama8258/deconvolution_giggle_db
test_dir=/scratch/Users/zama8258/giggle_test
bgzip=/Users/zama8258/giggle/lib/htslib/bgzip
giggle=/Users/zama8258/giggle/bin/giggle
# "$giggle" search -i "$data_dir"_b -q \
		# 					/scratch/Users/zama8258/deconascent/src/regions.bed.gz \
		# 					-o > "$data_dir"/overlaps.txt
<"$data_dir"/overlaps.txt sed \
 -e 's?#/scratch/Users/zama8258/deconvolution_giggle_db/??' \
 -e 's/.sorted_split_bidir_cov_filtered.bed.gz//' | \
		grep -v overlaps:0 | grep -v overlaps.bed.gz > "$data_dir"/overlaps_cleaned.txt
xz --threads=16 --memory=32G --keep \
	 -6 -z "$data_dir"/overlaps_cleaned.txt
