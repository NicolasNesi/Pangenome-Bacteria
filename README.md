# Pangenome analyses for Bacteria

Pipeline for Pangenome analyses of Bacteria Whole Genome Amplification with Illumina reads in HPC Austral (group 2019013)

## Create List_samples_RUN.txt ##
To create ```List_samples_RUN.txt```, a .txt file is needed with only one sample name per line. Name must start with [a-zA-Z] have NO space and NO underscore:
```
S23-MTX-COV-NOSE
S59-MTX-HRV-LUNG
```
## Configuration file ##
The configuration file is needed for the Pangenome pipeline to function properly, where the following settings can be modified and customized to run the entire pipeline in a single command. Preferences of analyses can be modified to your projects specifications.

## Directory structure ##
The beginning directory structure should be:
```
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

*/ denotes directories
```

## Tools used
- Reads cleanning and QC [**fastp**]
- *DeNovo* genome assembly [**unicycler**]
- Genome assembly QC [**BandageNG** **Bowtie2** **checkm** **BUSCO**]
- Genome annotation [**bakta**]
- Average Nucleotide Identity [***average_nucleotide_identity.py***]
- Pangenome: get core, cloud and singleton orthologous CDSs [**Panaroo**]
- Alignment CDS [***translatorX.pl*** **pal2nal**]
- Recombination test [**PhiPack**]
- QC alignments, premature STOP codon, nucleotide diversity, Tajima's D statistic and KaKs [***premature_stop.pl*** **KaKs_Calculator**  **Dendropy**]
- Extraction intergenic regions [**Piggy**] SOON

 ## Outputs can be used for phylogenomics analyses (e.g. concatenation, MSC, coreSNP tree inferences), comparative genomics analyses (presence/absence genes, positive selection test) or GWAS test. ##

## Pipeline workflow ##

![Figure_Pipeline_Pangenome_Comparativegenomics](https://github.com/user-attachments/assets/1857b1e0-5df9-4dee-bdf5-abd867eaee1c)

