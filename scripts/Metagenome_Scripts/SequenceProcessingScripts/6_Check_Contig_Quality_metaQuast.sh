#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=700G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
##SBATCH --time=8-00:00:00     # 8 days, 0 hrs
#SBATCH --output=Check_Assembly_Quality_metaquast.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Check_Assembly_Quality_metaquast"
#SBATCH -p aronsonlab
# you can use any of the following: intel, batch, highmem, gpu

#cd /bigdata/aronsonlab/shared/HannahFreund/HannahFTemp/

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

# conda activate MGM_Quality
module load quast/5.1.0rc1 # ignore this and the sbatch commands if running this locally or not on a Slurm system

if [[ ! -d ./ContigAssembly_Quality_Results ]]; then
    mkdir ./ContigAssembly_Quality_Results
fi

for i in *_contigs.fasta;
do
    file=$(basename $i)
    SAMPLE=${file%_contigs*}
    echo "Checking ${SAMPLE} assembled CONTIGS with metaQUAST"
    
    metaquast.py -t 4 -f -o ./${SAMPLE}_metaQUAST_contigs --gene-finding --mgm ${file}
    #metaquast.py -t 4 -f -o ./${SAMPLE}_metaQUAST_contigs --glimmer ${file}

    # if using sam file...
    #metaquast.py -t 4 --glimmer -o ./${SAMPLE}_metaQUAST_wsam --sam ${SAMPLE}_contigs_aln_pe.sam ${file}
    # if using bam file...
    # python metaquast.py -t 4 -f -o ./${SAMPLE}_metaQUAST_wbam_${today} --gene-finding --bam ${SAMPLE}_sorted.bam ${SAMPLE}_contigs.fasta
    # metaquast.py -t 4 -f -o ./${SAMPLE}_metaQUAST_wbam --gene-finding --bam ${SAMPLE}_sorted.bam ${file}
    
done

mv *_metaQUAST_contigs ./ContigAssembly_Quality_Results

# metaQUAST notes...
# http://quast.sourceforge.net/docs/manual.html#faq_q15 | https://github.com/ablab/quast/
# --gene-finding == -f -- Enables gene finding. Affects performance, thus disabled by default. By default, we assume that the genome is prokaryotic, and apply GeneMarkS for gene finding. If the genome is eukaryotic, use --eukaryote option to enable GeneMark-ES instead. If the genome is fungal, use --fungus option to run GeneMark-ES in a special mode for analyzing fungal genomes. If it is a metagenome (you are running metaquast.py), MetaGeneMark is used. You can also force MetaGeneMark predictions with --mgm option described below.
# must have a sam or bam file for each contig or scaffold file included in command or else it won't run

## Usage
# ./quast.py test_data/contigs_1.fasta \
           test_data/contigs_2.fasta \
        -r test_data/reference.fasta.gz \
        -g test_data/genes.txt \
        -1 test_data/reads1.fastq.gz -2 test_data/reads2.fastq.gz \
        -o quast_test_output
## Output
# report.txt      summary table
# report.tsv      tab-separated version, for parsing, or for spreadsheets (Google Docs, Excel, etc)
# report.tex      Latex version
# report.pdf      PDF version, includes all tables and plots for some statistics
# report.html     everything in an interactive HTML file
# icarus.html     Icarus main menu with links to interactive viewers
# contigs_reports/        [only if a reference genome is provided]
#   misassemblies_report  detailed report on misassemblies
#   unaligned_report      detailed report on unaligned and partially unaligned contigs
# k_mer_stats/            [only if --k-mer-stats is specified]
#   kmers_report          detailed report on k-mer-based metrics
# reads_stats/            [only if reads are provided]
#   reads_report          detailed report on mapped reads statistics

