library(tidyverse)
library(VennDiagram)
library(gplots)
library(pheatmap)
library(RSkittleBrewer)
library(ggVennDiagram)
library(purrr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(biomaRt)

setwd("C:/Users/hp/Desktop/HPV/GSE70462/RNAseq_data_analysis/GSE70463")
source("GSE70462.R")
setwd("C:/Users/hp/Desktop/HPV/GSE72536/RNASeq_data/GSE72536_rdata")
source("GSE72536.R")
setwd("C:/Users/hp/Desktop/HPV/GSE250305/RNASEQ_GSE250305")
source("GSE250305.R")
setwd("C:/Users/hp/Desktop/HPV/GSE211322/RNA_data_analysi/RnA_seq322/")
source("GSE211322.R")

ensg462 = annotated_462df$ensgene
ensg536 = annotated_536df$ensgene
ensg305 = annotated_305df$ensgene
ensg322 = annotated_322df$ensgene
## intersect genomic data

common_HPV_ensg = Reduce(intersect, list(ensg305, ensg462, ensg536, ensg322))


anno305_diff = anno305_df3$ensgene
anno462_diff = anno462_df3$ensgene
anno536_diff = anno536_df3$ensgene
anno322_diff = anno322_df3$ensgene

## Common significant genes

common_gene_diff = Reduce(intersect,list(anno305_diff,anno462_diff,
                                         anno536_diff, anno322_diff))

list_DE_gene = list(GSE250305= anno305_diff, 
                    GSE70462= anno462_diff, 
                    GSE72536= anno536_diff,
                    GSE211322= anno322_diff)

##  venn diagram

v = ggVennDiagram(list_DE_gene) +
  ggtitle("Venn diagram of overlapping DEG genes")

v
GO_results <- enrichGO(gene = common_gene_diff,
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", 
                       ont = "BP")

dotplot(GO_results, title = "Overlapping Genes")

entgene_DE = getBM(attributes = c("entrezgene_id"),
                   filters =c("ensembl_gene_id"),
                   values = common_gene_diff,
                   mart = ensemble111)

entgene_DE = as.character(entgene_DE$entrezgene_id)



ego_DE = enrichGO(gene = entgene_DE,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  readable = TRUE)

ego_DE
view(summary(ego_DE))
barplot(ego_DE, title = "Overlapping Genes")
dotplot(ego_DE, title = "Overlapping Genes")

## MF
ego_DE_MF = enrichGO(gene = entgene_DE,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  readable = TRUE)

ego_DE_MF
view(summary(ego_DE_MF))
barplot(ego_DE_MF, title = "Overlapping Genes Molecular function")
dotplot(ego_DE_MF, title = "Overlapping Genes Molecular function")


## individual gene expression

par(mfrow=c(1,4))
plotCounts(dds305, gene = "ENSG00000161860", 
           intgroup = "condition", main = "GSE250305") 

plotCounts(dds462, gene = "ENSG00000161860", 
           intgroup = "condition", main = "GSE70462")

plotCounts(dds536, gene = "ENSG00000161860", 
           intgroup = "condition", main = "GSE72536")

plotCounts(dds322, gene = "ENSG00000161860", 
           intgroup = "condition", main = "GSE211322")


## overlapping gene identification for PPI

ensemble111 = useEnsembl(biomart = "ensembl", version = 111)
view(listDatasets(ensemble111))
ensemble111 = useDataset("hsapiens_gene_ensembl", mart = ensemble111 )
view(listAttributes(ensemble111))
view(listFilters(ensemble111))
annotation_overlap = getBM(attributes = c("ensembl_gene_id", 
                                          "chromosome_name",
                                          "start_position", 
                                          "end_position",
                                          "strand",
                                          "gene_biotype",
                                          "external_gene_name",
                                          "description"),
                           filters =c("ensembl_gene_id"),
                           values = common_gene_diff,
                           mart = ensemble111)

write_tsv(annotation_overlap, "overlapping_genes.txt")

annotation_overlap$external_gene_name

write_tsv(common_gene_diff, "common_geneDiff.txt")



## finding overlapping gene in each dataset

anno_305_com = anno305_df3$ensgene %in% common_gene_diff
anno_305_com = anno305_df3[anno_305_com,]


anno_322_com = anno322_df3$ensgene %in% common_gene_diff
anno_322_com = anno322_df3[anno_322_com,]

anno_462_com = anno462_df3$ensgene %in% common_gene_diff
anno_462_com = anno462_df3[anno_462_com,]

anno_536_com = anno536_df3$ensgene %in% common_gene_diff
anno_536_com = anno536_df3[anno_536_com,]

write_tsv(anno_305_com, "anno_305_com.txt")
write_tsv(anno_322_com, "anno_322_com.txt")
write_tsv(anno_462_com, "anno_462_com.txt")
write_tsv(anno_536_com, "anno_536_com.txt")





