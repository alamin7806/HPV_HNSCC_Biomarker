#!/usr/bin/bash


# File: RNA seq analysis
# downloading SRA file (GSE211322)


var=$(tail -n +2 SraRunTable_GSE211322.txt | cut -d ',' -f 1)

for i in ${var}
        do
        if [[ -f ${i}.fastq.gz ]]
        then
                echo "the ${i} is present in the directory"
        else
                echo "downloading sra entry: ${i}"
                # Downloading Sra entry from GEO dataset
                fastq-dump --gzip --defline-qual '+' ${i}
                echo "completed downloading  ${i}"
        fi

        done

# Pre alignment QC of raw data
fastqc *.fastq.gz
multiqc .


# Downloading human reference transcript file

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz

# unzip reference file
gunzip gencode.v47.transcripts.fa.gz

# indexing salmon file

salmon index -t gencode.v47.transcripts.fa -p 12 -i salmon_index_1 --gencode

## Run Salmon
var=$(tail -n +2 SraRunTable_GSE211322.txt | cut -d ',' -f 1)
for i in ${var}
        do
                SRX=${i}.fastq.gz
                echo ${SRX}

salmon quant -i salmon_index_1 -l A -r ${SRX}  -o ${i}
        done











