#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G # < if you are threading, ask for whole memory and not memory per CPU
##SBATCH --time=20:00:00     # 20 hrs
#SBATCH --output=Rename_Contigs_10.10.22.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="RenameContigs"
#SBATCH -p aronsonlab

for i in *_contigs.fasta;
do
    awk -i inplace '/^>/{print ">c_" ++i; next}{print}' ${i}
done

# Helps with downstream processing specifically with gene assignment per contig and the naming conventions used

# changes long name to "c_1, c_2, c_3 etc
