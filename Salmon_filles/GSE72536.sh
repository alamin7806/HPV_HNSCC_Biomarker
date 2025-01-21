#!/usr/bin/bash


# File: RNA seq analysis
# downloading SRA file (GSE72536)
# Activate conda env
conda activate msa

var=$(tail -n +2 SraRunTable_GSE72536.txt | cut -d ',' -f 1)

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


## run salmon

tail -n +2 SraRunTable_536.txt | cut -d ',' -f 20 | sort | uniq | head -n 23 > geo_accesion.txt
var=$(cat geo_accesion.txt)
for i in ${var}

        do
                sra=$(grep ${i} SraRunTable_536.txt | cut -d ',' -f 1)
                sra=$(echo ${sra} | sed 's/ /.fastq.gz /g')
                SRX=${sra}.fastq.gz

salmon quant -i salmon_index_1 -l A -r ${SRX}  -o ${i}
        done





























