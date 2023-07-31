echo "Submitting Jobs"
parallel --progress sbatch gen_subset.sbatch {}  ::: $(seq 1 99)
