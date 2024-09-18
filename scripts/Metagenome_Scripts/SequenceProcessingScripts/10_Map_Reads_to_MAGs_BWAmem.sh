#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --time=7-00:00:00     # 7 days, 0 hrs
#SBATCH --mem=600G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --output=BWA_MEM_map_reads_to_genes.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Mapping reads to MAG bins with BWA-MEM"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu

module load bwa-mem2/2.2.1
module load samtools/1.14
module load bwa/0.7.17

# Make directory to store results
if [[ ! -d ./BWA_Map_Reads_to_Bins ]]; then
    mkdir ./BWA_Map_Reads_to_Bins
fi


for i in ./Good_Bins/*.fa;
do
    file=$(basename $i)
    SAMPLE=${file%.fa*}
    #echo ${SAMPLE} " -- map reads to MAG bins"

    ### Map UN-NORMALIZED reads back using local alignment
    if [[ ! -f ./BWA_Map_Reads_to_Bins/${SAMPLE}_bin_aln.sam ]]; then
    
        ## First, build the index aka reference for mapping reads to from the contigs fasta file
        bwa index ${i}
        
        ### Map trimmed, non-error corrected, un-normalized reads back using local alignment
        if [[ ${SAMPLE} == *"bin."* ]]; then
            # ensures that specific genes in bins are mapped to original contigs fasta file
            OGSample=${SAMPLE%_bin*}
            
            bwa mem ${i} ${OGSample}_R1_clean.fastq ${OGSample}_R2_clean.fastq -t 8 > ./BWA_Map_Reads_to_Bins/${SAMPLE}_bin_aln.sam
        fi

    fi
       
done

for i in ./BWA_Map_Reads_to_Bins/*_bin_aln.sam;
do
    file=$(basename $i)
    SAMPLE=${file%_bin*}
  
    if [[ ! -f ./BWA_Map_Reads_to_Bins/${SAMPLE}_bins_unsort.bam ]]; then
        
        ## Convert SAM file to BAM file with samtools
        samtools view -S -b ${i} > ${SAMPLE}_bin_unsort.bam  ## views & converts SAM to BAM file
        samtools sort ${SAMPLE}_bin_unsort.bam -o ${SAMPLE}_bin_sorted.bam ## sorts BAM file; sort alignments by leftmost coordinates, or by read name when -n is used
        samtools index  ${SAMPLE}_bin_sorted.bam ## indexes BAM file
        samtools flagstat -@ 8 -O tsv ${SAMPLE}_bin_sorted.bam > ${SAMPLE}_bin_stats.tsv
        samtools coverage  ${SAMPLE}_bin_sorted.bam -o ${SAMPLE}_bin_coverage.tsv
        samtools depth  ${SAMPLE}_bin_sorted.bam -o ${SAMPLE}_bin_depth.tsv
        
        #mv ${SAMPLE}_bin_unsort.bam ${SAMPLE}_bin_sorted.bam ${SAMPLE}_bin_sorted_stats.tsv ${SAMPLE}_bin_sorted_coverage.tsv ${SAMPLE}_bin_sorted_depth.tsv ./BWA_Map_to_Bins/

    fi
    
    
    
    if [[ ! -d ./BWA_Map_to_Bins/SamTools_Results ]]; then
        mkdir ./BWA_Map_to_Bins/SamTools_Results
        cp ${SAMPLE}_bin_stats.tsv ${SAMPLE}_bin_coverage.tsv ${SAMPLE}_bin_depth.tsv ./BWA_Map_to_Bins/SamTools_Results/
    else
        cp ${SAMPLE}_bin_stats.tsv ${SAMPLE}_bin_coverage.tsv ${SAMPLE}_bin_depth.tsv ./BWA_Map_to_Bins/SamTools_Results/
    fi
    
    
done

# * You can only index BAM files on position, and only when the data is sorted by position to begin

