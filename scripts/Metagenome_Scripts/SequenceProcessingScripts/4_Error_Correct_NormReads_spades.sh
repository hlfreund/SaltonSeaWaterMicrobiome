#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=200G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=800G
#SBATCH --time=5-00:00:00     # 5 days, 0 hrs
#SBATCH --output=error_correction_spades.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Error correcting clean reads with Spades"
#SBATCH -p highmem # partition name

Path="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_7.12.2022/Sequences_Analysis"

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load spades/3.15.4

for i in *_R1_norm.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    if [[ ! -d ./${SAMPLE}_error_corrected ]]; then
        spades.py -1 ${SAMPLE}_R1_norm.fastq -2 ${SAMPLE}_R2_norm.fastq -o ${SAMPLE}_error_corrected -t 16 --meta --only-error-correction
    fi

    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"

done

echo -e "Finished \aerror correction on mgm reads with SPades"


for i in ${Path}/spades_error_corrected/corrected/*.fastq.gz;
do
    rename _norm _norm_EC $i
    cp *.fastq.gz ${Path}/
done

## Spades assembly - https://github.com/ablab/spades#meta
## Notes for Spades...
## spades.py -1 read1.fq -2 read2.fq --merged merged.fq -o spades_test <<< -s is the flag to indicate merged reads as a "library"
## Error corrected sequences wil be found in output_directory/corrected/ – files from read error correction

## Spades Output
#The full list of <output_dir> content is presented below:

#scaffolds.fasta – resulting scaffolds (recommended for use as resulting sequences) *****
#contigs.fasta – resulting contigs
#assembly_graph.fastg – assembly graph
#contigs.paths – contigs paths in the assembly graph
#scaffolds.paths – scaffolds paths in the assembly graph
#before_rr.fasta – contigs before repeat resolution
#corrected/ – files from read error correction
#   configs/ – configuration files for read error correction
#   corrected.yaml – internal configuration file
#   Output files with corrected reads
#params.txt – information about SPAdes parameters in this run
#spades.log – SPAdes log
#dataset.info – internal configuration file
#input_dataset.yaml – internal YAML data set file
#K<##>/– directory containing intermediate files from the run with K=<##>. These files should not be used as assembly results; use resulting contigs/scaffolds in files mentioned above.

## HPCC (Cluster systems) note
##SBATCH --mem-per-cpu=500G # --mem=900gb --> how to request for total allocation of mem rather than mem per cpu; include in job submission command
