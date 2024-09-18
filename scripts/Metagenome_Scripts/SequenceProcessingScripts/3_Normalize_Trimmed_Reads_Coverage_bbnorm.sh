#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G
#SBATCH --time=10-00:00:00     # 7 days, 0 hrs
#SBATCH --output=Normalize_Read_Coverage.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Normalize Read Coverage"
#SBATCH -p aronsonlab

# you can use any of the following: intel, batch, highmem, gpu

Path="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_7.12.2022/Sequences/"
# ^ replace with the path where your trimmed MGM sequence files are stored

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load BBMap/38.95

if [[ ! -d ${Path}/Normalized_Seqs ]]; then
    mkdir ${Path}/Normalized_Seqs
fi

# Normalized Trimmed Metagenome Reads

for i in *_R1_clean.fastq;
do
    SAMPLE=$(echo ${i} | sed "s/_R1_clean.fastq//")
    echo "Normalizing trimmed reads from" ${SAMPLE}_R1_clean.fastq ${SAMPLE}_R2_clean.fastq
    
    if [[ ! -f ${Path}/Normalized_Seqs/${SAMPLE}_R1_norm.fastq ]] && [[ ! -f ${Path}/Normalized_Seqs/${SAMPLE}_R2_norm.fastq ]]; then
        bbnorm.sh in1=${SAMPLE}_R1_clean.fastq in2=${SAMPLE}_R2_clean.fastq out1=${Path}/Normalized_Seqs/${SAMPLE}_R1_norm.fastq out2=${Path}/Normalized_Seqs/${SAMPLE}_R2_norm.fastq target=100 min=5
        # This will run 2-pass normalization to produce an output file of reads with an average depth of 100x. Reads with an apparent depth of under 5x will be presumed to be errors and discarded.
    
    
        # cp ${SAMPLE}_R1_norm.fastq ${SAMPLE}_R2_norm.fastq ${Path}/Normalized_Seqs/

    fi
        
done

# ** Normalization happens AFTER trimming reads

#bbnorm.sh in=reads.fq out=normalized.fq target=100 min=5



