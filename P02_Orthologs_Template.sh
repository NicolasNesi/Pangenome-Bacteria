#!/usr/bin/env bash
#SBATCH --job-name=P2_Orthologs             # Job name
#SBATCH --output=P2_Orthologs.%a.%A.out     # File to which stdout will be written
#SBATCH --error=P2_Orthologs.%a.%A.err      # File to which stderr will be written
#SBATCH --partition=tlong                   # Partition
#SBATCH --cpus-per-task=24             # Number of cpus ask per task
#SBATCH --time=12-00:00                # Runtime in DD-HH:MM
#SBATCH --mem=500G                     # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                # BEGIN,END,ALL
#SBATCH --mail-user=mymail             # Email address

# ---------------------------------
# Version v1.1
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
MODE="${4}"
ALIGN="${5}"
ALIGNER="${6}"
CORETHRES="${7}"
OUTGROUP="${8}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="${FOLDER}/Bakta"
SOFT="/soft/2019013/Logiciels"
SCRIPTS="${FOLDER}/Scripts"

RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Create directories and copy input data
# ---------------------------------
# safety for failing steps
set -euo pipefail

mkdir -p "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/Input_gff
cp "${PATH1}"/*/*.gff3 "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/Input_gff/

# Genereate setting file of the analysis:
printf "Analysis run the: $RUNDATE\n
Working directory: $PWD\n
Paths variables are: group $GROUP project $PROJECT and analyses of the run $RUN\n
Average Nucleotide Identity performed with pyani\n
Pangenome orthologs CDSs obtained using Panaroo with the mode: $MODE and alignment for $ALIG\n"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
eval "$(conda shell.bash hook)"
conda activate pyani_env
# include:
# pyani=0.2.13.1
# csvtk=0.30.0
# ---------------------------------

# ---------------------------------
# Bandage Statistics
# ---------------------------------
echo "Create Bandage tsv table of all samples"

cd "${FOLDER}"/Assembly_Unicycler

for f in *.temp2; do awk 'NR==1 {split(FILENAME,a,"_"); print a[1] "\n" $0; next} { print $0}' "${f}" > "${f%.temp2}".tsv; done

paste -d "\t" "${FOLDER}"/Scripts/Bandage_Template.tsv *_statistics.tsv | csvtk transpose -tT > Bandage_Assembly_Statistics_AllSamples.tsv

rm *.tmp2 *_statistics.tsv

cd "${FOLDER}"/
# ---------------------------------

# ---------------------------------
# Average Nucleotide Identity
# ---------------------------------
echo "Calculate Average Nucleotide identity"

average_nucleotide_identity.py \
--indir "${FOLDER}"/pyani/"${RUN}"/ani_input/ \
--outdir "${FOLDER}"/pyani/"${RUN}"/ANIm_output/ \
--method ANIm \
--verbose \
--graphics \
--gformat pdf,svg \
--write_excel \
--jobprefix "${RUN}" \
--logfile "${FOLDER}"/pyani/"${RUN}"/ANIm.log
# ---------------------------------

# ---------------------------------
conda deactivate
conda activate orthoproka_env
# include:
# panaroo=1.5.1
# checkm-genome=1.2.3
# ---------------------------------

# ---------------------------------
# CheckM statistics
# ---------------------------------
echo "Create CheckM tsv table of all samples"

checkm lineage_wf \
--extension fasta \
"${FOLDER}"/Assembly_Unicycler \
--tab_table \
--threads "${THREADS}" \
"${FOLDER}"/CheckM/ > "${FOLDER}"/CheckM/CheckM_Assembly_Statistics_AllSamples.tsv

mv "${FOLDER}"/CheckM/bins/*.tsv "${FOLDER}"/CheckM/

rm "${FOLDER}"/CheckM/lineage.ms "${FOLDER}"/CheckM/checkm.log
rm -R "${FOLDER}"/CheckM/bins
# ---------------------------------

# ---------------------------------
# Orthologs with Panaroo and Ptolemy
# ---------------------------------
echo "Get Pangenome using Panaroo"

cd "${FOLDER}"/Orthologs/"${RUN}"/Panaroo

mkdir -p results_qc_"${RUN}"
mkdir -p outputs_panaroo
#mkdir -p outputs_ptolemy

echo "Runnning Panaroo QC"

panaroo-qc \
--input Input_gff/*.gff3 \
--graph_type all \
--ref_db /dlocal/home/2019013/Databases/Mash/refseq.genomes+plasmid.k21s1000.msh \
--out_dir results_qc_"${RUN}" \
--threads "${THREADS}"

echo "Running initial estimate of pangenome using $MODE mode, alignment of $ALIGN and a core threshold of $CORETHRES"

panaroo \
--input Input_gff/*.gff3 \
--out_dir outputs_panaroo \
--clean-mode "${MODE}" \
--threshold 0.98 \
--len_dif_percent 0.98 \
--family_threshold 0.7 \
--family_len_dif_percent 0.0 \
--alignment "${ALIGN}" \
--aligner "${ALIGNER}" \
--codons \
--core_threshold "${CORETHRES}" \
--no_clean_edges \
--remove-invalid-genes \
--refind-mode default \
--search_radius 5000 \
--refind_prop_match 0.5 \
--merge_paralogs \
--threads "${THREADS}"

echo "Update GFF files"

panaroo-generate-gffs \
--input Input_gff/*.gff3 \
--out_dir outputs_panaroo \
--threads "${THREADS}"

echo "Filtering of presence/absence gene matrices"

cd "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/outputs_panaroo

panaroo-filter-pa \
--input gene_presence_absence.csv \
--out_dir ./ \
--type pseudo,length,frag

cd "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/outputs_panaroo/unaligned_dna_sequences

ls *.fasta | sed 's/.fasta//g' > "${FOLDER}/Scripts/List_multifasta_$RUN.txt"

#echo "Step2 indexing of sequences with Ptolemy"

# index graph
#java -jar "${SOFT}"/ptolemy.jar extract \
#--genomes "${SCRIPTS}"/ptolemy-input-list_"${RUN}".txt \
#--output-database "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/outputs_ptolemy \
#--use-cds \
#--verbose
#
#python "${SCRIPTS}"/createSA.py "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/outputs_panaroo/ "${FOLDER}"/Orthologs/"${RUN}"/Panaroo/outputs_ptolemy/
#
#java -jar "${SOFT}"/ptolemy.jar canonical-quiver \
#-s outputs_ptolemy/syntenic_anchors.txt \
#--db output_ptolemy \
#--output-directory o"${FOLDER}"/Orthologs/"${RUN}"/Panaroo/output_ptolemy \
#--verbose

# ---------------------------------
conda deactivate
conda activate phylo_env
# include:
# iqtree=2.3.6
# ---------------------------------

# ---------------------------------
# Dated Core Genome Tree
# ---------------------------------
#echo "Inference of dated phylogeny based on core genome for estimate Pangenome size and gene gain loss rates analyses"

#iqtree2 \
#-s core_gene_alignment_filtered.aln \
#--date "${SCRIPTS}"/Date-file_"${RUN}".tsv \
#-o "${OUTGROUP}" \
#--date-no-outgroup \
#--date-ci 100
# ---------------------------------

# ---------------------------------
conda deactivate
conda activate orthoproka_env
# include:
# panaroo=1.5.1
# ---------------------------------

# ---------------------------------
# Infinetely many Genes model
# ---------------------------------
echo "Estimate gene gain and loss rates according to the Infinitely Many Genes (IMG) model"

mkdir img_results
mkdir fmg_results

panaroo-img \
--pa gene_presence_absence.Rtab \
--tree dated_phylogeny.tree \
--model fixed \
--verbose \
--out_dir img_results

echo "Estimate gene gain and loss rates according to the Finitely Many Genes (FMG) model"

panaroo-fmg \
--pa gene_presnece_absence.Rtab \
--tree dated_phylogeny.tree \
--nboot 100 \
--verbose \
--threads "${THREADS}" \
--outfile fmg_results
# ---------------------------------





