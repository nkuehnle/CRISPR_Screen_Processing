#!/bin/bash

# Note, ### indicates aspects that you may need to customize
# Create a folder called Raw_FASTQ
# Place your raw FASTQ files inside. You can rename the folder for your new project
# Activate any relevant Conda environments
### I created a single crisprenv running MAGeCK 0.5.9.4, cutadapt 3.1, bowtie 1.3, and samtools 1.9
### conda create -n crisprenv -c conda-forge -c bioconda -c default mageck=0.5.9.4 cutadapt=3.1 bowtie=1.3.0 samtools=1.9 bzip2=1.0.8 six
### Note, the six package is only needed for drugz (provides python 2/3 compatibility support)
# If Samtools passes errors about an SSL library, try running conda install -c conda-forge -c bioconda samtools=1.9 bzip2=1.0.8 again inside your new environment.
### Works as of December 2020
# Navigate via cd to the project directory and then use the sh command to run this script

################################# Variables to set ##########################################
### For ease, you may want to change the defaults.
CPU_CORES=1 # Set default cores to 1
PROJ_NAME="${PWD##*/}" # Set default name to the name of the current directory
LIBRARY=Brunello # Set default library to Brunello
ADAPTER="g cgaaacaccg"
NO_ALIGN="rc"
TRIM_LEN=20
MIN_LEN=20
NORM_METHOD="median"
###

while getopts p:n:l:a:x:t:m:c flag # Names of flags
do
    case "${flag}" in
        p) CPU_CORES=${OPTARG};; # CPU_CORES Flag -p
        n) PROJ_NAME=${OPTARG};; # PROJ_NAME Flag -n
        l) LIBRARY=${OPTARG};; # LIBRARY Flag -l
        a) ADAPTER=${OPTARG};; # Adapter sequence (g to indicate 5', a to indicate 3')
        x) NO_ALIGN=${OPTARG};; # Alignment direction to exclude
        t) TRIM_LEN=${OPTARG};; # Length to trim reads to after adapter removal
        m) MIN_LEN=${OPTARG};; # Minimum length of output trimmed reads
        c) NORM_METHOD=${OPTARG};; #Normalization method for MAGeCK
    esac
done

############################# Creating Directory/File Structure #############################
# Do not need to change anything here.
cd "$(dirname "$0")" || exit
mkdir Trim_FASTQ
mkdir Bowtie
mkdir MAGeCK

########################################## Trimming #########################################
# This will iterate through each raw FASTQ file and trim w/ cutadapt using the settings below.
### -g cut 5' Adaptor
### -l cut to length (from  the 3' end if positive, 5' if negative)
### -m minimum length -M maximum length filtering
# --cores = threads/cores
### --discard-untrimmed will get rid of reads without any adapters
# -o output file location and/or name
# Lastly give the name of the file we want to act on

for FNAME in Raw_FASTQ/*.fastq.gz
do
    SAMPLE_NAME="$(printf "$FNAME" | sed -e "s/^Raw_FASTQ\///" | sed -e "s/.fastq.gz$//")" #This will remove Raw_FASTQ/ and .fastq.gz from each file
    cutadapt -"$ADAPTER" --length "$TRIM_LEN" --discard-untrimmed -m "$MIN_LEN" --cores="$CPU_CORES" -o ./Trim_FASTQ/"$SAMPLE_NAME"_trim.fastq.gz "$FNAME" ### Tweak the arguments as neccessary, i.e. adaptors etc
done

################################## Creating Bowtie Indexes ###################################
# This will create indexes for the alignment.
# --threads is the number of threads/cores to use

awk '{print ">"$1"\n"$2}' Libraries/"$LIBRARY".txt > Bowtie/"$LIBRARY".fa # awk command to convert tab-delimited library file to fasta format
bowtie-build -f --threads "$CPU_CORES" Bowtie/"$LIBRARY".fa Bowtie/bowtie1_ind_"$LIBRARY" #Generate bowtie indexes

##################################### Bowtie Alignment #######################################
# This will iterate through each trimmed FASTQ file and align to the brunello library using Bowtie 1.3 and the following settings.
# -p is cores/threads for bowtie
# -S indicates we will be piping the output to samtools (no SAM intermediary file is made, instead we create a BAM file
# -x index naming convention/location
# -q FASTQ input file name
### -n is the number of mismatches in the see to allow (overall misamtches if seed length = the length of the sgRNA
### -l seed length (using 20, meaning the entire sgRNA is used as the seed
### -n + -l allows us to allow 1 mismatch in a 20 bp initial match or 1 mismatch overall
### --norc prevents us from looking for reverse compliment alignments
### -k ouput only the first alignment
### --best along with -k 1 ensures that the first alignment is the best, aka only the first/best alignment is output
# -@ threads/cores for samtools

for FNAME in Trim_FASTQ/*.fastq.gz
do
    SAMPLE_NAME="$(printf "$FNAME" | sed -e "s/^Trim_FASTQ\///" | sed -e "s/_trim.fastq.gz$//")" # This will remove Raw_FASTQ/ and __trim.fastq.gz from each file
    bowtie -p "$CPU_CORES" -S -x Bowtie/bowtie1_ind_"$LIBRARY" -q "$FNAME" -n 1 -l 20 -no"$NO_ALIGN" -k 1 --best | samtools view -@ "$CPU_CORES" -bS - > Bowtie/"$SAMPLE_NAME".bam ### Tweak the arguments as desired
done
################################## MAGeCK Count ###################################
# This will create the count table with one column for each initial fastq.gz file
### You may or may not want to look at MAGeCK's documentation to enable additional options.

cd Bowtie || exit

SAMPLE_LABELS=$(ls *.bam | xargs echo | sed 's/ /,/g' | sed 's/.bam//g') #Create a list of comma-separated labels from the bam file names (w/o .bam)
BAM_LIST=$(ls *.bam | xargs echo) #Create a list of space-separated .bam files
mageck count -l ../Libraries/Brunello.txt -n ../MAGeCK/"$PROJ_NAME" --sample-label "$SAMPLE_LABELS" --fastq $BAM_LIST ### Tweak the arguments as desired

cd ../MAGeCK

mkdir Counts
mkdir Tests
mkdir Counts/Logs
mkdir Counts/R_Output
mkdir Counts/Other
mkdir Tests/Logs
mkdir Tests/R_Output
mkdir Tests/Other
mkdir Tests/sgRNA_Results
mkdir Tests/Gene_Results

mv *.log ./Counts/Logs
mv *.R ./Counts/R_Output
mv *.Rmd ./Counts/R_Output
mv *.Rnw ./Counts/R_Output
mv *.count_normalized.txt ./Counts/Other
mv *.countsummary.txt ./Counts/Other

while IFS=$'\t' read -r Output_Name Treatment Control
do
    mageck test -k "$PROJ_NAME".count.txt -t "$Treatment" -c "$Control" -n "$Output_Name" --norm-method "$NORM_METHOD" --control-sgrna ../Libraries/Controls/"$LIBRARY"_controls.txt
done < ../MAGeCK_Tests_Table.txt

mv *.log ./Tests/Logs
mv *.sgrna_summary.txt ./Tests/sgRNA_Results
mv *.gene_summary.txt ./Tests/Gene_Results
mv *.R ./Tests/R_Output
mv *.Rmd ./Tests/R_Output
mv *.Rnw ./Tests/R_Output
