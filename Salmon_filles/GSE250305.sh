#!/usr/bin/bash


# File: RNA seq analysis
# downloading SRA file (GSE250305)


# activate the right conda env
conda activate msa

# get input data (uniprot ids)
SRA=$(tail -n +2 SraRunTable_GSE250305.txt | cut -d ',' -f 1)

# This is a loop for downloading the data

for i in ${SRA}
    do
        prefetch ${i}
        if [ -f ${i}.fastq.gz ]
            then
                echo "${i} already finished"
        else
                echo "(o) Convert SRA entry: ${i}"
                # downloading SRA entry
                fastq-dump --gzip --defline-qual '+' ${i}/${i}.sra
                echo "(o) Done convert  ${i}"
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

var=$(tail -n +2 SraRunTable_GSE250305.txt | cut -d ',' -f 1)
for i in ${var}
        do
                SRX=${i}.fastq.gz
                echo ${SRX}

salmon quant -i salmon_index_1 -l A -r ${SRX}  -o ${i}
        done


# pull transcript and gene id from the reference gene

grep -P -o "ENST\d{11}" gencode.v47.transcripts.fa > enst.txt
grep -P -o "ENSG\d{11}" gencode.v47.transcripts.fa > ensg.txt
paste -d ',' enst.txt ensg.txt > gene_id.txt





