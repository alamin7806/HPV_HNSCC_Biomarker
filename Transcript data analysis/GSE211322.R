library(readr)
library(magrittr)
library(dplyr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(plotly)
library(tibble)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gplots)
library(AnnotationDbi)
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)

view(read_csv("SraRunTable_GSE211322.txt"))


table_322 = read_csv("SraRunTable_GSE211322.txt") %>%
  dplyr::select(Run, BioSample, Experiment, `Sample Name`, source_name, 
                genotype, tissue)
table_322 =  table_322[-c(1,7,17,18),]
view(table_322)


colnames(table_322)[4] = "Sample_name"

sample_file_322 = paste0(pull(table_322, Run), "/quant.sf")
sample_file_322
names(sample_file_322) = pull(table_322, Run)
sample_file_322
ref_gene = read_csv("gene_id.txt", col_names = c("enst", "ensg"))

## from transcript to gene level count
count_data_322 = tximport(files = sample_file_322,
                          type = "salmon",
                          tx2gene = ref_gene,
                          ignoreTxVersion = TRUE)



view(table_322)
table_322 = as.data.frame(table_322)
condition = factor(table_322$source_name)

repeat_HPV = c(rep("HPV_Negative", times =2),
               rep(c("HPV_Positive","HPV_Negative"), each =4),
               rep("HPV_Positive", times =10))
condition = factor(repeat_HPV)
table_322$condition = condition

## construct a DESeq DataSet from count data

deseq_data322 = DESeqDataSetFromTximport(txi = count_data_322,
                                         colData = table_322,
                                         design = ~condition)

## pre-filtering data

deseq_data322_norm = rowSums(counts(deseq_data322)) >=10
deseq_data322_norm = deseq_data322[deseq_data322_norm,]

## Differential expression analysis

dds322 = DESeq(deseq_data322_norm)
vst322 = varianceStabilizingTransformation(dds322)
plotPCA(vst322, intgroup = "condition") + ggtitle("GSE211322") + geom_point(size= 4)

dds322_res = results(dds322)
summary(dds322_res)


dds322_res_df = as.data.frame(dds322_res)
dds322_res_df = rownames_to_column(dds322_res_df, var = "ensgene")
dds322_f1 = dplyr::filter(dds322_res_df, complete.cases(dds322_res_df))
dds322_f2 = dplyr::filter(dds322_f1, padj <0.05)
dds322_f3 = dplyr::filter(dds322_f2, abs(log2FoldChange) > 1)
dds322_f1$test = dds322_f1$padj < 0.05 & abs(dds322_f1$log2FoldChange) > 1

g322 = ggplot(dds322_f1, aes(x= log2FoldChange, 
                             y=-log10(padj),
                             name=ensgene)) +
  geom_point(aes(colour = test), size =2, alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE211322")

g322
ggplotly(g322)

ensemble111 = useEnsembl(biomart = "ensembl", version = 111)
view(listDatasets(ensemble111))
ensemble111 = useDataset("hsapiens_gene_ensembl", mart = ensemble111 )
view(listAttributes(ensemble111))
view(listFilters(ensemble111))
annotation_322 = getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", 
                                      "end_position",
                                      "strand",
                                      "gene_biotype",
                                      "external_gene_name",
                                      "description"),
                       filters =c("ensembl_gene_id"),
                       values = dds322_f1$ensgene,
                       mart = ensemble111)
dim(annotation_322)
dim(dds322_f1)
annotated_322df = left_join(dds322_f1, annotation_322,
                            by = c("ensgene" = "ensembl_gene_id"))
dim(annotated_322df)

g322_annotated= ggplot(annotated_322df, aes(x= log2FoldChange, 
                                            y=-log10(padj),
                                            name= external_gene_name)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE211322")
g322_annotated

ggplotly(g322_annotated)
anno322_df2 = dplyr::filter(annotated_322df, padj <0.05)
anno322_df3 = dplyr::filter(anno322_df2, abs(log2FoldChange) > 1)


## GO analysis

entgene_322 = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = anno322_df3$ensgene,
                    mart = ensemble111)

entgene_322 = as.character(entgene_322$entrezgene_id)


ego_322 = enrichGO(gene = entgene_322,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP")

ego_322
view(summary(ego_322))
barplot(ego_322, title = "GSE211322")
dotplot(ego_322, title = "GSE211322")



ekg322 = enrichKEGG(gene = entgene_322)

barplot(ekg322, title = "GSE211322")
view(ekg322)
write_tsv(anno322_df3, "DE Genes_GSE211322.txt")
write_tsv(dds322_f3, "DE Genes_filtered_GSE211322.txt")





















































































































































