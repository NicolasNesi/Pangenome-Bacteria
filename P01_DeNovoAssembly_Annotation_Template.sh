#!/usr/bin/env bash
#SBATCH --job-name=P01_DeNovoAssembly_Annotation              # Job name
#SBATCH --output=P01_DeNovoAssembly_Annotation.%a.%A.out      # File to which stdout will be written
#SBATCH --error=P01_DeNovoAssembly_Annotation.%a.%A.err       # File to which stderr will be written
#SBATCH --partition=tcourt                                    # Partition
#SBATCH --cpus-per-task=1                                     # Number of cpus ask per task
#SBATCH --time=01-00:00                                       # Runtime in DD-HH:MM
#SBATCH --mem=100G                                            # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                                       # BEGIN,END,ALL
#SBATCH --mail-user=mymail                                    # Email address

# ---------------------------------
# Version v1.2
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, INSERM UMR 1311 DYNAMICURE
# Date: 25/11/2024
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
QSCORE="${4}"
GENUS="${5}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="${FOLDER}/Output_Illumina/${RUN}"
SCRIPT="/home/2019013/PARTAGE/Various_Scripts_for_various_purposes"

# Databases
BAKTADB="/dlocal/home/2019013/Databases/Bakta/db"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "${FOLDER}/Scripts/List_samples_$RUN.txt")"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
# safety for failing steps
set -euo pipefail

cp "${PATH1}"/Fastq/"${SAMPLE}"_R*.fastq.gz "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/
gunzip *.fastq.gz

# Generate setting file of the analysis:
printf "Analysis run the: $RUNDATE \n
Working directory: $PWD \n
Paths variables are: group $GROUP project $PROJECT and analyses of the run $RUN \n
Filtering using fastq with a minimum qscore of $QSCORE \n
De novo assembly using Unicycler \n
genomes completeness and contamination with CheckM and BUSCO \n
Annotation using Bakta for genus $GENUS \n"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
eval "$(conda shell.bash hook)"
conda activate orthoproka_env
# Include:
# BandageNG=2022.9-dev
# bowtie2=2.5.4
# checkm-genome=1.2.3
# fastp=0.23.4
# unicycler=0.5.1
# SPAdes=4.0.0
# samtools=1.20
# ---------------------------------

# ---------------------------------
# Short read QC
# ---------------------------------
echo "QC and filtering of short reads for sample: $SAMPLE"

fastp \
--in1 "${SAMPLE}_R1.fastq" \
--in2 "${SAMPLE}_R2.fastq" \
--detect_adapter_for_pe \
--length_required 50 \
--qualified_quality_phred "${QSCORE}" \
--low_complexity_filter \
--thread 16 \
--json "${SAMPLE}"_ReportFastp.json \
--html "${SAMPLE}"_ReportFastp.html \
--report_title "$SAMPLE fastq report" \
--out1 "${SAMPLE}_CleanedReads_R1.fastq" \
--out2 "${SAMPLE}_CleanedReads_R2.fastq" \
--unpaired1 "${SAMPLE}_CleanedReads_unpaired.fastq" \
--unpaired2 "${SAMPLE}_CleanedReads_unpaired.fastq" 2>&1 | tee -a "${SAMPLE}_fastp_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# De novo assembly from Illumina reads
# ---------------------------------
echo "De novo assembly with Unicycler for sample: $SAMPLE"
 
unicycler \
--short1 "${SAMPLE}_CleanedReads_R1.fastq" \
--short2 "${SAMPLE}_CleanedReads_R2.fastq" \
--unpaired "${SAMPLE}_CleanedReads_unpaired.fastq" \
--out "${SAMPLE}"_unicycler_assembly/ \
--min_fasta_length 500 \
--depth_filter 0.25 \
--keep 0 \
--threads "${THREADS}" 2>&1 | tee -a "${SAMPLE}_Unicycler_run${RUNDATE}.log"

mv "${SAMPLE}"_unicycler_assembly/assembly.fasta "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta"
mv "${SAMPLE}"_unicycler_assembly/assembly.gfa "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.gfa"
mv "${SAMPLE}"_unicycler_assembly/unicycler.log "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_unicycler.log"

echo "Calculate statistics of Unicycler assembly for sample: $SAMPLE"

Bandage info "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.gfa" > "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_Bandage_statistics.temp"
grep -E 'Node count:|Edge count:|Smallest edge overlap \(bp\):|Largest edge overlap \(bp\):|Total length \(bp\):|Total length no overlaps \(bp\):|Dead ends:|Percentage dead ends:|Connected components:|Largest component \(bp\):|Total length orphaned nodes \(bp\):|N50 \(bp\):|Shortest node \(bp\):|Lower quartile node \(bp\):|Median node \(bp\):|Upper quartile node \(bp\):|Longest node \(bp\):|Median depth:|Estimated sequence length \(bp\):' "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_Bandage_statistics.temp" | awk '{print $NF}' > "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_Bandage_statistics.temp2"
rm "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_Bandage_statistics.temp"
# ---------------------------------

# ---------------------------------
# Quality average coverage
# ---------------------------------
echo "Mapping back reads against assembly for sample: $SAMPLE"

bowtie2-build "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" ./"${SAMPLE}"_assembly

bowtie2 \
-1 "${SAMPLE}_CleanedReads_R1.fastq" \
-2 "${SAMPLE}_CleanedReads_R2.fastq" \
-U "${SAMPLE}_CleanedReads_unpaired.fastq" \
-x ./"${SAMPLE}"_assembly \
--very-sensitive \
--dovetail \
--time \
--threads "${THREADS}" \
-S "${SAMPLE}_assembly.sam" 2>&1 | tee -a "${SAMPLE}_Bowtie2Assembly_run${RUNDATE}.log"

echo "Generate coverage statistics"

samtools view --threads 96 -S --bam "${SAMPLE}_assembly.sam" -o "${SAMPLE}_assembly.bam"
samtools sort --threads 96 "${SAMPLE}_assembly.bam" -o "${SAMPLE}_assembly.sorted.bam"
samtools index "${SAMPLE}_assembly.sorted.bam"
samtools coverage "${SAMPLE}_assembly.sorted.bam" > "${SAMPLE}_assembly_coverage.txt"

# Klebsiella genome OK: - contigs <300
#                       - dead ends <200
#                       - N50 >50,000
#                       - average coverage >30x
#                       - genome size 4,969,898 to 6,132,846 bp
#                       - completeness >95%
#                       - contamination <5%
# ---------------------------------

# ---------------------------------
# Genome annotation
# ---------------------------------
conda activate bakta_env
# Include:
# bakta=1.9.4

echo "Annotation of sample: $SAMPLE"

bakta \
--db "${BAKTADB}" \
"${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" \
--prefix "${SAMPLE}" \
--genus "${GENUS}" \
--output "${SAMPLE}"_bakta/ 2>&1 | tee -a "${SAMPLE}_bakta_run${RUNDATE}.log"

python "${SCRIPT}"/convert_gbff_to_gbk.py -input "${SAMPLE}"_bakta/"${SAMPLE}".gbff -output "${SAMPLE}"_bakta/"${SAMPLE}".gbk
# ---------------------------------

# ---------------------------------
# Genome completeness contamination
# ---------------------------------
echo "Genome completeness and contamination of sample: $SAMPLE"

conda activate orthoproka_env

checkm lineage_wf \
--extension fasta \
"${SAMPLE}"_unicycler_assembly/ \
--tab_table \
--threads "${THREADS}" \
"${SAMPLE}"_checkm/

conda activate BUSCO_env
# Include:
# busco=5.7.0

busco \
--in "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" \
--mode genome \
--lineage_dataset /dlocal/home/2019013/Databases/BUSCO/lineages/enterobacterales_odb10/ \
--out "${SAMPLE}"_BUSCO/ \
--cpu "${THREADS}"

python3 /soft/2019013/Conda_env/BUSCO_env/bin/generate_plot.py --working_directory "${SAMPLE}"_BUSCO/
mv "${SAMPLE}"_BUSCO/busco_figure.png "${SAMPLE}"_BUSCO/"${SAMPLE}"_busco_figure.png
# ---------------------------------

# ---------------------------------
# Output directories and Move data
# ---------------------------------
if [[ ! -d "${FOLDER}/CleanedReads_Illumina" ]]; then
   mkdir -p "${FOLDER}/CleanedReads_Illumina"
fi

mv "${SAMPLE}"_CleanedReads_*.fastq "${FOLDER}/CleanedReads_Illumina/"

if [[ ! -d "${FOLDER}/Assembly_Unicycler" ]]; then
   mkdir -p "${FOLDER}/Assembly_Unicycler"
fi

mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.gfa" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_unicycler.log" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_Bandage_statistics.temp2" "${FOLDER}/Assembly_Unicycler"

if [[ ! -d "${FOLDER}/CheckM" ]]; then
   mkdir -p "${FOLDER}/CheckM"
fi

mv "${SAMPLE}"_checkm/storage/bin_stats_ext.tsv "${FOLDER}"/CheckM/"${SAMPLE}"_qa_statistics.txt
mv "${SAMPLE}"_checkm/storage/bin_stats.analyze.tsv "${FOLDER}"/CheckM/"${SAMPLE}"_general_statistics.txt

if [[ ! -d "${FOLDER}/BUSCO" ]]; then
   mkdir -p "${FOLDER}/BUSCO"
fi

mv "${SAMPLE}"_BUSCO/short_summary*.txt "${FOLDER}"/BUSCO/"${SAMPLE}"_short_summary.txt
mv "${SAMPLE}"_BUSCO/"${SAMPLE}"_busco_figure.png "${FOLDER}"/BUSCO/"${SAMPLE}"_busco_figure.png

if [[ ! -d "${FOLDER}/Bakta" ]]; then
   mkdir -p "${FOLDER}/Bakta"
fi

mv "${SAMPLE}"_bakta "${FOLDER}/Bakta"

if [[ ! -d "${FOLDER}/Bowtie2" ]]; then
   mkdir -p "${FOLDER}/Bowtie2"
fi

mv "${SAMPLE}_assembly.sam" "${FOLDER}/Bowtie2"
mv "${SAMPLE}_assembly.bam" "${FOLDER}/Bowtie2"
mv "${SAMPLE}"_assembly.sorted.ba* "${FOLDER}/Bowtie2"
mv "${SAMPLE}_assembly_coverage.txt" "${FOLDER}/Bowtie2"
mv "${SAMPLE}_Bowtie2Assembly_run${RUNDATE}.log" "${FOLDER}/Bowtie2"
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
ls -R > "${FOLDER}/Scripts/List_files_remaining_in_LocalWorkingDir_for_${SAMPLE}_run${RUNDATE}.txt"
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P01_DeNovoAssembly_Annotation."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}/Scripts/Sbatch_log"
rm -f "${FOLDER}"/Scripts/P01_DeNovoAssembly_Annotation.*.*.out 
# ---------------------------------
