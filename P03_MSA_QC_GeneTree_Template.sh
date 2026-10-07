#!/usr/bin/env bash
#SBATCH --job-name=P03_MSA-QC-GeneTree             # Job name
#SBATCH --output=P03_MSA-QC-GeneTree.%a.%A.out     # File to which stdout will be written
#SBATCH --error=P03_MSA-QC-GeneTree.%a.%A.err      # File to which stderr will be written
#SBATCH --partition=2tcourt                        # Partition
#SBATCH --cpus-per-task=48                         # Number of cpus ask per task
#SBATCH --time=00-12:00                            # Runtime in DD-HH:MM
#SBATCH --mem=50G                                  # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                            # BEGIN,END,ALL
#SBATCH --mail-user=mymail                         # Email address


# ---------------------------------
# Version v1.2
# ---------------------------------
# author: Nicolas Nesi
# University Caen Normandy, INSERM UMR 1311 DYNAMICURE
# Date: 26/11/2024
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# Positional parameters
GROUP="${1}"
PROJECT="${2}"
RUN="${3}"
ALIGNMENT1="${4}"
ALIGNMENT2="${5}"
CODE="${6}"
BOOTSTRAP="${7}"

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="$FOLDER/Orthologs/$RUN/Panaroo/outputs_panaroo/unaligned_dna_sequences"
SCRIPTS="$FOLDER/Scripts"
SOFT="/soft/2019013/Logiciels"
UTILS="/home/2019013/PARTAGE/Various_Scripts_for_various_purposes"

# Array job
GENE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "${FOLDER}/Scripts/List_multifasta_$RUN.txt")"
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used:$THREADS"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
eval "$(conda shell.bash hook)"
conda activate phylo_env
# include:
# clustalw=2.1
# mafft=7.520
# muscle=5.1.linux64
# prank=v.170427  #!!Doesn't work with the current version of translatorX!!
# phipack=1.1
# seqkit=2.8.2
# iqtree=2.3.6
# ---------------------------------

# ---------------------------------
# Check for duplicate sequence in multifasta
# ---------------------------------
echo "Check and removing duplicate sequences if exist for gene: $GENE"

rename ".fasta" "_multifasta.fasta" "${PATH1}"/*.fasta

seqkit rmdup --by-name "${PATH1}"/"${GENE}"_multifasta.fasta > "${PATH1}"/"${GENE}"_multifasta_nodup.fasta
# ---------------------------------

# ---------------------------------
# Alignment of Protein Coding gene 
# ---------------------------------
echo "Alignment of protein coding gene: $GENE"

"${UTILS}"/translatorX.pl \
-i "${PATH1}"/"${GENE}"_multifasta_nodup.fasta \
-o "${PATH1}"/"${GENE}_${ALIGNMENT2}" \
-c "${CODE}" \
-p "${ALIGNMENT1}"

mkdir -p "${FOLDER}"/Alignments/"${ALIGNMENT2}"

mv "${PATH1}"/"${GENE}_${ALIGNMENT2}".nt_ali.fasta "${FOLDER}"/Alignments/"${ALIGNMENT2}"/"${GENE}_${ALIGNMENT2}_nt_align.fasta"
mv "${PATH1}"/"${GENE}_${ALIGNMENT2}".aa_ali.fasta "${FOLDER}"/Alignments/"${ALIGNMENT2}"/"${GENE}_${ALIGNMENT2}_prot_align.fasta"
mv "${PATH1}"/"${GENE}_${ALIGNMENT2}".html "${FOLDER}"/Alignments/"${ALIGNMENT2}"/"${GENE}_${ALIGNMENT2}_translatorX.html"
mv "${PATH1}"/"${GENE}_${ALIGNMENT2}".nt*_ali.fasta "${FOLDER}"/Alignments/"${ALIGNMENT2}"/
rm "${PATH1}"/"${GENE}_${ALIGNMENT2}".aa_based_codon_coloured.html
rm "${PATH1}"/"${GENE}_${ALIGNMENT2}".aaseqs*
# ---------------------------------

# ---------------------------------
# Convert alignment in codon data 
# ---------------------------------
echo "Running PAL2NAL for gene: $GENE"

cd "${FOLDER}"/Alignments/"${ALIGNMENT2}"

mkdir -p pal2nal_checks

"${UTILS}"/pal2nal.pl "${GENE}_${ALIGNMENT2}_prot_align.fasta" "${GENE}_${ALIGNMENT2}_nt_align.fasta" -output fasta > pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal.fasta.tmp"

seqkit seq -w0 pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal.fasta.tmp" > pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal.fasta"

rm pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal.fasta.tmp"
# ---------------------------------

# ---------------------------------
# Check for premature STOP
# ---------------------------------
echo "Check for the presence of premature STOP in alignment for gene: $GENE"

perl "${UTILS}"/premature_stops.pl -i pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal.fasta" -o pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_StopOK.fasta"
# ---------------------------------

# ---------------------------------
# Calculate Nonsynonymous position 
# ---------------------------------
echo "Calculate Nonsynonymous position in alignment of gene: $GENE"

mkdir -p KaKs

AMAS.py convert -d dna -f fasta -u phylip -i pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_StopOK.fasta" -c 48 --check-align

"${SOFT}"/Kaks_Calculator-main/KaKs_Calculator2-1.0/bin/Linux/KaKs_Calculator \
-i pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_stopOK.fasta-out.phy" \
-o KaKs/"${GENE}_${ALIGNMENT2}_pal2nal.fasta.kaks" \
-m YN

echo "To get Nonsynonymous Ka value run the following command: \n  
while read -r l; do sed -n '2p' "$l"_Mafft_pal2nal.fasta.kaks | awk -F '\t' -v fname="$l" '{ print fname " Ka: " $3}' ; done < Scripts/List_multifasta_Klebsiella.txt > NonSynonymous_perGeneAlignment.txt"
# ---------------------------------

# ---------------------------------
# Calculate Recombination of gene alignment
# ---------------------------------
echo "Estimate recombinant for the aligment gene: $GENE"

mkdir -p Recombination

Phi \
-r pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_stopOK.fasta-out.phy" \
-p 1000 \
-o \
-v > Recombination/"${GENE}_${ALIGNMENT2}_phipack_RecombinationTest.txt"

# gene has recombinant sites if pvalues of the three test P < 0.05
# ---------------------------------

# ---------------------------------
# Calculate nucleotide diversity (π) and Tajima's D statistic
# ---------------------------------
echo "Calculate nucleotide diversity (π) and Tajima's D statistic for gene : $GENE"

mkdir -p Dendropy

python "${UTILS}"/dendropy_popgenstat.py -i pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_StopOK.fasta" -o Dendropy/"${GENE}_${ALIGNMENT2}_NucleoDivTajima.txt"
# ---------------------------------

# ---------------------------------
# Make ML Gene tree 
# ---------------------------------
echo "Infereing maximum Likelihood gene tree for gene: $GENE"

mkdir -p GeneTree

cd GeneTree/

iqtree \
-s "${FOLDER}"/Alignments/"${ALIGNMENT2}"/pal2nal_checks/"${GENE}_${ALIGNMENT2}_pal2nal_stopOK.fasta" \
-m MFP \
--ufboot "${BOOTSTRAP}" \
-bnni \
--prefix "${GENE}" \
-T AUTO
# ---------------------------------

# ---------------------------------
mv "${FOLDER}"/Scripts/P03_MSA-QC-GeneTree."${SLURM_ARRAY_TASK_ID}".* "${FOLDER}/Scripts/Sbatch_log"
rm -f "${FOLDER}"/Scripts/P03_MSA-QC-GeneTree.*.*.out 
# ---------------------------------

