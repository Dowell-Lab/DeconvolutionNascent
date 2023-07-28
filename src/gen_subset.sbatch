#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=2gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euxo pipefail

# Load Modulea
module load samtools
module load subread

# Set up directory
workDir=/scratch/Users/zama8258/deconvolution_titration
mkdir -p "$workDir"
pushd "$workDir" || exit

# Shrink names for readability
db="/Shares/dbnascent"
filetype="crams"
suffix=".sorted.cram"

# For each file
# Subsample to size provided on stdin
refSeq=/scratch/Shares/dowell/genomes/hg38/hg38.fa
parallel -j 10 --link --progress samtools view -C -T "$refSeq" \
				 -s {2} -o "$workDir"/sample_"$1"_{3}.cram {1} ::: \
				 "$db"/Jiang2018multi/"$filetype"/SRR6789175"$suffix" \
				 "$db"/Fei2018ndf/"$filetype"/SRR7010982"$suffix" \
				 "$db"/Andrysik2017identification/"$filetype"/SRR4090102"$suffix" \
				 "$db"/Dukler2017nascent/"$filetype"/SRR5364303"$suffix" \
				 "$db"/Smith2021peppro/"$filetype"/SRR10669536"$suffix" \
				 "$db"/Zhao2016high/"$filetype"/SRR3713700"$suffix" \
				 "$db"/Danko2018dynamic/"$filetype"/SRR6780907"$suffix" \
				 "$db"/Chu2018chromatin/"$filetype"/SRR7616132"$suffix" \
				 "$db"/Ikegami2020phosphorylated/"$filetype"/SRR10601203"$suffix" \
				 "$db"/Core2014analysis/"$filetype"/SRR1552485"$suffix" \
				 ::: \
				 $(Rscript -e "cat(seq((1/99)*""$1"", 1, length.out=10) / sum(seq((1/99)*""$1"", 1, length.out=10)))") \
				 ::: \
				 "$(seq 1 10)"
samtools merge --threads 10 -o "merged_""$1"".bam" \
				 "$workDir"/sample_"$1"_4.cram \
				 "$workDir"/sample_"$1"_6.cram \
				 "$workDir"/sample_"$1"_7.cram \
				 "$workDir"/sample_"$1"_8.cram \
				 "$workDir"/sample_"$1"_9.cram \
				 "$workDir"/sample_"$1"_10.cram
# "$workDir"/sample_"$1"_1.cram \
		# "$workDir"/sample_"$1"_2.cram \
		# "$workDir"/sample_"$1"_3.cram \
		# "$workDir"/sample_"$1"_5.cram \
		# "$workDir"/sample_"$1"_*.cram

# Generate SAF for featureCounts
db=/scratch/Shares/dowell/dbnascent/out/meta_analysis/mumerge/archive
refFile="$db"/bidirectionals_dreg_tfit/hg38_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.bed
safFile="$workDir"/bidirs_"$1".saf
awk -v OFS='\t' '{print $1":"$2"-"$3, $1, $2, $3, "+"}' "$refFile" > "$safFile"
# Run featureCounts
featureCounts \
		-T 9 \
		-s 0 \
		-F 'SAF' \
		-a "$safFile" \
		-o "$workDir"/"merged_""$1""_counts.txt" \
		"$workDir"/"merged_""$1"".bam"
