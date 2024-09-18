#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=700G
#SBATCH --time=5-00:00:00     # 5 days, 0 hrs
#SBATCH --output=metaBAT_contig_binning.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Binning contigs from assembled mgms with metaBAT"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu


#today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

conda activate MAG_Binning # used conda env since we did not have correct metabat version installed on CC at the time
#module load metabat/0.32.4 # load the module on your computing cluster system (ignore this if running locally)

for i in *_contigs.fasta;
do
    f=$(basename $i)
    SAMPLE=${f%_contigs*}
    echo ${SAMPLE} " -- binning contigs in this metagenome"
    
    runMetaBat.sh -t 4 -m 1500 ${SAMPLE}_contigs.fasta ${SAMPLE}_contigs_sorted.bam

done

for i in *_contigs.fasta;
do
    f=$(basename $i)
    SAMPLE=${f%_contigs*}
    
    find . -name "${SAMPLE}*fa" -exec mv -t ./${SAMPLE}_bins_metabat {} +

done

# Rename metabat bin directories
for dir in *metabat-bins1500
do
    file=$(basename $dir)
    SAMPLE=${file%_contigs*}
    mv $dir ${SAMPLE}_bins

done

# Rename actual bins so you know which sample they came from
#for dir in *_good_bins;
#do
#    file=$(basename $dir)
#    SAMPLE=${file%_good*}
#
#    for bin in ${dir}/*;
#    do
#
#        bin_name=$(basename $bin)
        #SAMPLE=${file%_good*}
#        rename bin ${SAMPLE}_bin ./${dir}/${bin_name}
#
#    done
#
#done

conda deactivate
### * bam files must be sorted first
### ${SAMPLE}_depth.txt -- includes mean and variance of base coverage depth

## Running locally with metabat installed....
# jgi_summarize_bam_contig_depths --outputDepth ${SAMPLE}_depth.txt ${SAMPLE}_bam.bam
# metabat2 -i ${SAMPLE}_assembly.fa.gz -a ${SAMPLE}_depth.txt -o ${SAMPLE}_bins/bin -v
# for f in *.fa; do mv -i -- "$f" "${f//contigs.fasta.metabat-bins-_-t_4_-m_1500/Bin}"; done

