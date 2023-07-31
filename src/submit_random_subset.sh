# for iter in $(seq 2 128)
# do
# 		echo "Submitting $iter"
# 		sbatch gen_subset_random.sbatch "$iter"
# done
parallel sbatch gen_subset_random.sbatch ::: $(seq 1 128)
