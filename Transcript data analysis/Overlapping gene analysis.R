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

ensg462 = annotated_462df$ensgene
ensg536 = annotated_536df$ensgene
ensg305 = annotated_305df$ensgene

## intersect genomic data

common_HPV_ensg = Reduce(intersect, list(ensg305, ensg462, ensg536))

ensg462_commonGene = annotated_462df[annotated_462df$ensgene
                                     %in% common_HPV_ensg, c("ensgene", "log2FoldChange")]

ensg536_commonGene = annotated_536df[annotated_536df$ensgene
                                     %in% common_HPV_ensg,c("ensgene", "log2FoldChange")]


ensg305_commonGene = annotated_305df[annotated_305df$ensgene
                                     %in% common_HPV_ensg,c("ensgene", "log2FoldChange")]

list_df = list(ensg305_commonGene, ensg462_commonGene,ensg536_commonGene)




common_folChange = purrr::reduce(list_df, left_join , by = "ensgene")
colnames(common_folChange)[c(2,3,4)] = c("GSE250305", "GSE70462","GSE72536")

correlate_gene = cor(common_folChange[,c("GSE250305", "GSE70462","GSE72536")], 
                     method = "pearson")


correlate_gene =  as.data.frame(correlate_gene)
rownames(correlate_gene) = colnames(correlate_gene)
rownames_to_column(correlate_gene, var = "GEO data set")

write_tsv(correlate_gene, "correlation of gene.txt")

anno305_diff = anno305_df3$ensgene
anno462_diff = anno462_df3$ensgene
anno536_diff = anno536_df3$ensgene

## Common significant genes

common_gene_diff = Reduce(intersect,list(anno305_diff,anno462_diff, anno536_diff))

list_DE_gene = list(GSE250305= anno305_diff, 
                    GSE70462= anno462_diff, 
                    GSE72536= anno536_diff)

##  venn diagram

v = ggVennDiagram(list_DE_gene) +
  ggtitle("Venn diagram of overlapping DEG genes")


GO_results <- enrichGO(gene = common_gene_diff,
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENSEMBL", 
                       ont = "BP")

entgene_DE = getBM(attributes = c("entrezgene_id"),
                   filters =c("ensembl_gene_id"),
                   values = common_gene_diff,
                   mart = ensemble111)

entgene_DE = as.character(entgene_DE$entrezgene_id)

entUni_DE = getBM(attributes = c("entrezgene_id"),
                  filters =c("ensembl_gene_id"),
                  values = common_HPV_ensg,
                  mart = ensemble111)
entUni_DE = as.character(entUni_DE$entrezgene_id)

ego_DE = enrichGO(gene = entgene_DE,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  universe = entUni_DE,
                  readable = TRUE)

ego_DE
view(summary(ego_DE))
barplot(ego_DE, title = "Overlapping Genes")
dotplot(ego_DE, title = "Overlapping Genes")




ekegg_DE = enrichKEGG(gene = entgene_DE,
                      universe = entUni_DE)

view(ekegg_DE)

## individual gene expression


plotCounts(dds305, gene = "ENSG00000077935", 
           intgroup = "condition", main = "GSE250305")

plotCounts(dds462, gene = "ENSG00000077935", 
           intgroup = "condition", main = "GSE70462")

plotCounts(dds536, gene = "ENSG00000077935", 
           intgroup = "condition", main = "GSE72536")


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

























