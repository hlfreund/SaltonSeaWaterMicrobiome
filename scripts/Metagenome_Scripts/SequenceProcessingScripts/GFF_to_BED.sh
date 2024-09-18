#!/bin/bash -l

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00     # 1 days, 0 hrs
#SBATCH --output=GFF_to_BED.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="GFF_to_BED_file2"
##SBATCH -p intel

# Convert GFF files to BED files

for i in *_genes.gff; # change to your output file names from when you ran CheckM
do
    file=$(basename $i)
    SAMPLE=${file%_genes*}
    
    
    if [[ ! -f ${SAMPLE}_genes.bed ]]; then
        awk -F '\t' '/#/ {next} {OFS = FS} {print $1,$4,$5,$1=$1 "_" (++count[$1])}' ${i} > ${SAMPLE}_genes.bed
    fi
    
    # $1=$1 "_" (++count[$1]) -> count every unique ID and add counter after "_" and original ID name in column 1 ($1)
done
