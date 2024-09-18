#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=8 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
##SBATCH --mem=800G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=1-00:00:00  # 14 days, 0 hrs
#SBATCH --output=Parse_KOFamScan_Annotations_11.9.22.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Parse_KOFamScan_Annotations_11.9.22"
#SBATCH -p aronsonlab
# you can use any of the following: intel, batch, highmem, gpu


# Make directory to store results
if [[ ! -d ./Parsed_KOFamScan_Results ]]; then
    mkdir ./Parsed_KOFamScan_Results
fi

# input into KoFam_Scan will be .faa files from Prokka or Prodigal (FASTA file of translated CDSs found in the contigs into Protein sequences)

# Parse through functions with Python script written by Mike Lee
for FILE in *_functions_kfs.txt;
do
    f=$(basename $FILE)
    SAMPLE=${f%_functions*}
    
    if [[ ! -f ${SAMPLE}_parsed_fxns.tsv ]]; then
        python Filter_KoFamScan_Results.py -i ${SAMPLE}_functions_kfs.txt -o ${SAMPLE}_parsed_fxns.tsv
        
        cp ${SAMPLE}_parsed_fxns.tsv ./Parsed_KOFamScan_Results
    fi
    
    
    if [[ ! -f ${SAMPLE}_parsed_fxns_mapper.txt ]]; then
        awk -F '\t' '{OFS = FS} NR>=2 {print $1,$2}' ${SAMPLE}_parsed_fxns.tsv > ${SAMPLE}_parsed_fxns_mapper.txt
        cp ${SAMPLE}_parsed_fxns_mapper.txt ./Parsed_KOFamScan_Results
    fi
done

#for FILE in ./*_functions.txt;
#do
#    f=$(basename $FILE)
#    SAMPLE=${f%_kfs*}
    
    # use awk to parse out functions identified by KoFamScan w/ evalue below 1e-3
#    awk -F"\t" 'NR==1{print;next}$6<=0.001' ${f} > ${SAMPLE}_kfs_below_threshold.txt

#done
