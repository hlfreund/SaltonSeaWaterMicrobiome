#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --time=5-00:00:00     # 5 days, 0 hrs
#SBATCH --mem=600G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --output=BWA_MEM_map_reads_to_contigs.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Mapping reads to assembled contigs with BWA-MEM"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu

#cd /bigdata/aronsonlab/shared/HannahFreund/HannahFTemp/ # change to your desired directory

#Path2Seqs="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SSea_POW_Mgm_Analysis/Checkm_Test/All_Contigs/Normalized_ErrorCorrected_Samples" # initiate your path for later

#today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load bwa-mem2/2.2.1
module load samtools/1.14
module load bwa/0.7.17

# Make directory to store results
if [[ ! -d ./BWA_Map_Reads_to_Contigs ]]; then
    mkdir ./BWA_Map_Reads_to_Contigs
fi


for i in *_contigs.fasta;
do
    file=$(basename $i)
    SAMPLE=${file%_contigs*}
    echo ${SAMPLE} " -- mapping reads to assembled contigs with BWA-MEM"
    

    ### Map NON-NORMALIZED reads back using local alignment
    if [[ ! -f ${SAMPLE}_contigs_aln.sam ]]; then
    
        ## First, build the index aka reference for mapping reads to from the contigs fasta file
        bwa index ${file}
        
        ### Map trimmed, non-error corrected, non-normalized reads back using local alignment
        bwa mem ${file} ${SAMPLE}_R1_clean.fastq ${SAMPLE}_R2_clean.fastq -t 8 > ${SAMPLE}_contigs_aln.sam
        #bwa-mem2 mem -t 8 ${SAMPLE}_c_assembly ${SAMPLE}_L001_R1_clean.fastq.gz ${SAMPLE}_L001_R2_clean.fastq.gz > ${SAMPLE}_contigs_aln-pe.sam

    fi
       
done

for i in *_contigs_aln.sam;
do
    file=$(basename $i)
    SAMPLE=${file%_contigs*}
  
    if [[ ! -f ./${SAMPLE}_contigs_unsort.bam ]]; then
        
        ## Convert SAM file to BAM file with samtools
        samtools view -S -b ${SAMPLE}_contigs_aln.sam > ${SAMPLE}_contigs_unsort.bam  ## views & converts SAM to BAM file
        samtools sort ${SAMPLE}_contigs_unsort.bam -o ${SAMPLE}_contigs_sorted.bam ## sorts BAM file; sort alignments by leftmost coordinates, or by read name when -n is used
        samtools index  ${SAMPLE}_contigs_sorted.bam ## indexes BAM file
        samtools flagstat -@ 8 -O tsv ${SAMPLE}_contigs_sorted.bam > ${SAMPLE}_contigs_stats.tsv
        samtools coverage  ${SAMPLE}_contigs_sorted.bam -o ${SAMPLE}_contigs_coverage.tsv
        samtools depth  ${SAMPLE}_contigs_sorted.bam -o ${SAMPLE}_contigs_depth.tsv
        
        mv ${SAMPLE}_contigs_unsort.bam ${SAMPLE}_contigs_sorted.bam ${SAMPLE}_contigs_sorted_stats.tsv ${SAMPLE}_contigs_sorted_coverage.tsv ${SAMPLE}_contigs_sorted_depth.tsv ./BWA_Map_Reads_to_Contigs/

    fi
    
    
    
    if [[ ! -d ./BWA_Map_Reads_to_Contigs/SamTools_Results ]]; then
        mkdir ./BWA_Map_Reads_to_Contigs/SamTools_Results
        cp ${SAMPLE}_contigs_stats.tsv ${SAMPLE}_contigs_coverage.tsv ${SAMPLE}_contigs_depth.tsv ./BWA_Map_Reads_to_Contigs/SamTools_Results/
    fi
    
    
done

# * You can only index BAM files on position, and only when the data is sorted by position to begin

