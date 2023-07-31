#!/bin/bash
#SBATCH --output=/scratch/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/zama8258/e_and_o/%x_%j.err
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=4gb
#SBATCH --mail-user=zama8258@colorado.edu

set -euxo pipefail

# Load Modulea
module load samtools
module load subread

# Variables
workDir=/scratch/Users/zama8258/deconvolution_gene_trials_300

# Shrink names for readability
db="/Shares/dbnascent"
filetype="crams"
suffix=".sorted.cram"

# Convert to bam for counting
refSeq=/scratch/Shares/dowell/genomes/hg38/hg38.fa
# parallel -j 8 --progress --link samtools view -b -T "$refSeq" \
		# 				 -o "$workDir"/{2}.bam {1} ::: \
		# 				 "$db"/Jiang2018multi/"$filetype"/SRR6789175"$suffix" \
		# 				 "$db"/Fei2018ndf/"$filetype"/SRR7010982"$suffix" \
		# 				 "$db"/Andrysik2017identification/"$filetype"/SRR4090102"$suffix" \
		# 				 "$db"/Dukler2017nascent/"$filetype"/SRR5364303"$suffix" \
		# 				 "$db"/Smith2021peppro/"$filetype"/SRR10669536"$suffix" \
		# 				 "$db"/Zhao2016high/"$filetype"/SRR3713700"$suffix" \
		# 				 "$db"/Danko2018dynamic/"$filetype"/SRR6780907"$suffix" \
		# 				 "$db"/Chu2018chromatin/"$filetype"/SRR7616132"$suffix" \
		# 				 "$db"/Ikegami2020phosphorylated/"$filetype"/SRR10601203"$suffix" \
		# 				 "$db"/Core2014analysis/"$filetype"/SRR1552485"$suffix" \
		# 				 ::: \
		# 				 Jiang2018multi \
		# 				 Fei2018ndf \
		# 				 Andrysik2017identification \
		# 				 Dukler2017nascent \
		# 				 Smith2021peppro \
		# 				 Zhao2016high \
		# 				 Danko2018dynamic \
		# 				 Chu2018chromatin \
		# 				 Ikegami2020phosphorylated \
		# 				 Core2014analysis

# Do counting
# Generate SAF
db=/scratch/Users/zama8258/deconascent
refFile="$db"/data/hg38_refseq.bed
safFile="$workDir"/hg38_refseq.saf
# awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$refFile" > "$safFile"
awk -v OFS='\t' '{ if ($3-$2 > 300) {if ($6 == "+") print $4, $1, $2+300, $3, $6; else print $4, $1, $2, $3-300, $6}}' "$refFile" > "$safFile"

# Run featureCounts on each individual file
# Only want to run once, so we can separate this out?
parallel -j 8 --progress --link featureCounts \
				 -T 1 \
				 -s 1 \
				 -F 'SAF' \
				 -a "$safFile" \
				 -o "$workDir"/{1}"_counts.txt" \
				 "$workDir"/{1}.bam \
				 ::: \
				 Jiang2018multi \
				 Fei2018ndf \
				 Andrysik2017identification \
				 Dukler2017nascent \
				 Smith2021peppro \
				 Zhao2016high \
				 Danko2018dynamic \
				 Chu2018chromatin \
				 Ikegami2020phosphorylated \
				 Core2014analysis
