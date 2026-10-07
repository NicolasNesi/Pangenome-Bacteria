#!/usr/bin/env bash

# -----------------------------------------------
# Control file to submit Pangenome Prokaryote Pipeline
# -----------------------------------------------

# -----------------------------------------------
# Version v1.0
# -----------------------------------------------
# author: Nicolas Nesi
# University Caen Normandy, INSERM UMR 1311 DYNAMICURE
# Date: 12/11/2024
# -----------------------------------------------

# -----------------------------------------------
# READ INSTRUCTION BEFORE SUBMISSION
# -----------------------------------------------
# Configuration file for the submission of the pipeline for analysis of Illumina metatranscriptomics data

# Start point is a folder /dlocal/home/2019013/Data/$GROUP/$PROJECT with:

#   - subfolder Output_Illumina/ where you will copy your Illumina output folder of the $RUN

#   - subfolder Scripts/ with the Pangenome Pipeline scripts, List_samples_$RUN.txt
# 
# Sample name format: Start with [a-zA-Z], NO space, NO underscore! (e.g. S23-MTX-COV-NOSE)
#
# .
# └── GROUP/
#     └── PROJECT/
#         ├── Output_Illumina/
#         │   └── RUN/
#         │       └── Fastq/
#         │           ├── Sample1_R1.fastq.gz
#         │           ├── Sample1_R2.fastq.gz
#         │           ├── Sample2_R1.fastq.gz
#         │           └── Sample2_R2.fastq.gz
#         └── Scripts/
#             ├── Config_Pipeline_Pangenome_Prokaryote_Template.sh
#             ├── P01_DeNovoAssembly_Annotation_Template.sh
#             ├── P02_Orthologs_Template.sh
#             ├── P03_MSA_QC_GeneTree_Template.sh
#             └── List_samples_RUN.txt
#
# -----------------------------------------------

# -----------------------------------------------
# Submission of the control script from Scripts/ folder
# -----------------------------------------------
# Type:
# bash Config_Pipeline_Pangenome_Prokaryote_Template.sh
# -----------------------------------------------

# -----------------------------------------------#
# -----------------------------------------------#
# EDIT THE FOLLOWING LINES OF THE CONFIG SCRIPTS #
# -----------------------------------------------#
# -----------------------------------------------#

# -----------------------------------------------
# SLURM EMAIL
# -----------------------------------------------
# $EMAIL: Indicate your email address for the Slurm SBATCH line --mail-user= (e.g. nicolas.nesi@unicaen.fr)
EMAIL=""
# -----------------------------------------------

# -----------------------------------------------
# SLURM ARRAY JOBS
# -----------------------------------------------
# $ARRAY: numbers to be given to sbatch argument --array= corresponding to the number of lines in your List_samples_$RUN.txt file
    # A full 96-plate: 1-96
    # Only run some samples: 1,12,35-40
ARRAY=""
# -----------------------------------------------

# -----------------------------------------------
# DIRECTORIES and INPUT FILES PATH
# -----------------------------------------------
# $GROUP: path to the folder in /dlocal/home/201913/Data/ - Options available: 
    # BacterioCaen
    # Colab 
    # Rouen 
    # ViroCaen
    # Labeo
GROUP=""
# $PROJECT: name of the Project main directory in /dlocal/home/201913/Data/$GROUP (e.g. DYNAMICH, Microbiosthme_alice, DynamiK)
PROJECT=""
# $RUN: name of Run, corresponding to the Illumina folder created (e.g. )
RUN=""
# $PATH-SCRIPTS: path to scripts and List_samples_$RUN.txt files
PATHSCRIPTS="/dlocal/home/2019013/Data/${GROUP}/${PROJECT}/Scripts"
# -----------------------------------------------

# -----------------------------------------------
# CLEANING READS
# -----------------------------------------------
# $QSCORE: value for the --qualified_quality_phred argument of Fastp (e.g. 30)
QSCORE=""
# -----------------------------------------------

# -----------------------------------------------
# GENOME ASSEMBLY
# -----------------------------------------------
# $GENUS: genus name of your sample for Bakta annotation (e.g. Klebsiella)
GENUS=""
# -----------------------------------------------

# -----------------------------------------------
# PANGENOME PANAROO
# -----------------------------------------------
# $MODE: choose the --clean-mode for analysis - Option available:
  # strict (present in at least 5% of genome)
  # moderate (present in at least 1% of genome)
  # sensitive (does not delete any genes)
MODE=""
# $ALIGN: choose outputs multiple sequence alignments (MSA) of core genes or all genes - Options available:
  # core
  # pan
ALIGN=""
# $ALIGNER: choose the aligner you want to generate MSA - Option available:
  # clustal
  # mafft (default)
  # prank
  # none (= unaligned sequences use this option as script P03 will do codon MSA)
ALIGNER=""
# $CORETHRES: the frequency of a gene in the your sample required to classify it as 'core' and aligned (e.g. presnt in at least 98% of samples: 0.98)
CORETHRES=""
# -----------------------------------------------

# -----------------------------------------------
# CODON ALIGNMENT TRANSLATORX
# -----------------------------------------------
# $ALIGNMENT1: Aligner to use with translatorX.pl - option available:
   # C  (Clustalw)
   # F  (maFft)
   # M  (muscle)
   # P  (Prank !! currently doesn't work)
ALIGNMENT1=""
# $ALIGNMENT2: same as above - indicate full name of the aligner - option available:
   # Clustalw
   # Mafft
   # Muscle
   # Prank
ALIGNMENT2=""
# $CODE: genetic code to sue for codon translation e.g. 1 (standard)
CODE=""  
# -----------------------------------------------

# -----------------------------------------------
# GENE TREE IQTREE
# -----------------------------------------------
# $BOOTSTRAP: number of replicate for ulfrafast bootstrap (>=1000)
BOOTSTRAP=""
# -----------------------------------------------

# -----------------------------------------------
# SELECT WHAT YOU WANT TO RUN
# -----------------------------------------------
# $ANALYSES : indicate if run all the Pipeline
	# all
	# step_1 (DeNovo Assembly and Annotation)
	# step_2 (Pangenome Orthologs)
	# step_3 (MSA and GeneTree)
ANALYSES=""
# -----------------------------------------------

# -----------------------------------------------
# Submission of the script Pipeline_Metatranscriptomics_Illumina_Template.sh
# -----------------------------------------------

# Ensure required variables are set
if [ -z "$EMAIL" ] || [ -z "$ARRAY" ] || [ -z "$PATHSCRIPTS" ] || [ -z "$GROUP" ] || [ -z "$PROJECT" ] || [ -z "$RUN" ] || \
   [ -z "$QSCORE" ] || [ -z "$GENUS" ] || [ -z "$MODE" ] || [ -z "$ALIGNER" ] || [ -z "$CORETHRES" ] || \
   [ -z "$ALIGNMENT1" ] || [ -z "$ALIGNMENT2" ] || [ -z "$CODE" ] || [ -z "$BOOTSTRAP" ] || [ -z "$ANALYSES" ]; then
  echo "One or more required variables are not set."
  exit 1
fi

# indicate your email address 
if [[ -n "$EMAIL" ]]; then
    if ! grep -q "#SBATCH --mail-user=$EMAIL" job_script.sh; then
        sed -i "s/#SBATCH --mail-user=mymail/#SBATCH --mail-user=$EMAIL/" "${PATHSCRIPTS}"/P*_Template.sh
        echo "Replaced 'mymail' with $EMAIL in Template scripts"
    else
        echo "Email $EMAIL is already present in Template scripts"
    fi
else
    echo "Warning: \$EMAIL is not set. Please provide an email address."
fi

set -euo pipefail

# dos2unix List_samples_$RUN.txt
dos2unix "${PATHSCRIPTS}"/List_samples_"${RUN}".txt
 
 # indicate your email address 
if [[ -n "$EMAIL" ]]; then
    if ! grep -q "#SBATCH --mail-user=$EMAIL" job_script.sh; then
        sed -i "s/#SBATCH --mail-user=mymail/#SBATCH --mail-user=$EMAIL/" "${PATHSCRIPTS}"/P*_Template.sh
        echo "Replaced 'mymail' with $EMAIL in P*_Template.sh scripts"
    else
        echo "Email $EMAIL is already present in P*_Template.sh scripts"
    fi
else
    echo "Warning: \$EMAIL is not set. Please provide an email address."
fi

if [[ "$ANALYSES" == "all" ]]; then
   # Step 1: Submission of the script P01_DeNovoAssembly_Annotation_Template.sh
   echo "STEP 1 - DeNovo assembly and annotation"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_DeNovoAssembly_Annotation_Template.sh \
   $GROUP $PROJECT $RUN $QSCORE $GENUS)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

   # Step 2: Submission of the script P02_Orthologs_Template.sh
   echo "STEP 2 - Get Pangenome and orthologs using Panaroo"
   jid_step2=$(sbatch --parsable --dependency=afterok:${jid_step1} --array=$ARRAY $PATHSCRIPTS/P02_Orthologs_Template.sh \
   $GROUP $PROJECT $RUN $MODE $ALIGN $ALIGNER $CORETHRES)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"
   
   # Step 3: Submission of the script P03_MSA_QC_GeneTree_Template.sh
   echo "STEP 3 - Creation of MSA for CDSs, check of premature STOP codon, recombination, calculate KaKs, nucleotide diversity and Tajima's D statistic and create gene tree"
   jid_step3=$(sbatch --parsable --dependency=afterok:${jid_step1} $PATHSCRIPTS/P03_MSA_QC_GeneTree_Template.sh \
   $GROUP $PROJECT $RUN $ALIGNMENT1 $ALIGNMENT2 $CODE $BOOTSTRAP)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

elif [[ "$ANALYSES" == "step_1" ]]; then
   # Step 1: Submission of the script P01_DeNovoAssembly_Annotation_Template.sh
   echo "STEP 1 - DeNovo assembly and annotation"
   jid_step1=$(sbatch --parsable $PATHSCRIPTS/P01_DeNovoAssembly_Annotation_Template.sh \
   $GROUP $PROJECT $RUN $QSCORE $GENUS)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step1" ]; then
       echo "Failed to submit Step 1 job."
       exit 1
   fi
   echo "Step 1 job ID: ${jid_step1}"

elif [[ "$ANALYSES" == "step_2" ]]; then
   # Step 2: Submission of the script P02_Orthologs_Template.sh
   echo "STEP 2 - Get Pangenome and orthologs using Panaroo"
   jid_step2=$(sbatch --parsable --dependency=afterok:${jid_step1} --array=$ARRAY $PATHSCRIPTS/P02_Orthologs_Template.sh \
   $GROUP $PROJECT $RUN $MODE $ALIGN $ALIGNER $CORETHRES)

   # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step2" ]; then
       echo "Failed to submit Step 2 job."
       exit 1
   fi
   echo "Step 2 job ID: ${jid_step2}"

elif [[ "$ANALYSES" == "step_3" ]]; then
   # Step 3: Submission of the script P03_MSA_QC_GeneTree_Template.sh
   echo "STEP 3 - Creation of MSA for CDSs, check of premature STOP codon, recombination, calculate KaKs, nucleotide diversity and Tajima's D statistic and create gene tree"
   jid_step3=$(sbatch --parsable --dependency=afterok:${jid_step1} $PATHSCRIPTS/P03_MSA_QC_GeneTree_Template.sh \
   $GROUP $PROJECT $RUN $ALIGNMENT1 $ALIGNMENT2 $CODE $BOOTSTRAP)
   
    # Check if the job submission was successful
   if [ $? -ne 0 ] || [ -z "$jid_step3" ]; then
       echo "Failed to submit Step 3 job."
       exit 1
   fi
   echo "Step 3 job ID: ${jid_step3}"

else
  echo "Invalid value for ANALYSES variable. Please set it to 'all', 'step_1', 'step_2' or 'step_3'."
  exit 1
fi
# -----------------------------------------------
