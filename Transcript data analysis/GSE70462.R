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
library(RSkittleBrewer)
library(gplots)
library(AnnotationDbi)




table_462 = read_csv("SraRunTable_GSE70462.txt") %>%
  dplyr::select(Run, dek_status, diagnosis, Experiment,
        `GEO_Accession (exp)`,hpv_status,`Sample Name`, source_name)
colnames(table_462)[c(5,7)] = c("GEO_Accession","Sample_name")
sample_file_462 = paste0(pull(table_462, Run), "/quant.sf")

names(sample_file_462) = pull(table_462, Run)

ref_gene = read_csv("gene_id.txt", col_names = c("enst", "ensg"))

## from transcript to gene level count
count_data_462 = tximport(files = sample_file_462,
                          type = "salmon",
                          tx2gene = ref_gene,
                          ignoreTxVersion = TRUE)

table_462 = as.data.frame(table_462)
condition = factor(rep(c("HPV_negative","HPV_positive"), each=2))

table_462$condition = condition

## construct a DESeq DataSet from count data
deseq_data462 = DESeqDataSetFromTximport(txi = count_data_462,
                                         colData = table_462,
                                         design = ~condition)
##pre-filtering data

deseq_data462_norm = rowSums(counts(deseq_data462)) >=10
deseq_data462_norm = deseq_data462[deseq_data462_norm,]

## Differential expression analysis

dds462 = DESeq(deseq_data462_norm)
vst462 = varianceStabilizingTransformation(dds462)
plotPCA(vst462, intgroup = "condition") + ggtitle("GSE70462")

dds462_res = results(dds462)
summary(dds462_res)

dds462_res_df = as.data.frame(dds462_res)
dds462_res_df = rownames_to_column(dds462_res_df, var = "ensgene")
dds462_f1 = dplyr::filter(dds462_res_df, complete.cases(dds462_res_df))
dds462_f2 = dplyr::filter(dds462_f1, padj <0.05)
dds462_f3 = dplyr::filter(dds462_f2, abs(log2FoldChange) > 1)
dds462_f1$test = dds462_f1$padj < 0.05 & abs(dds462_f1$log2FoldChange) > 1

## Volcano plot

g462 = ggplot(dds462_f1, aes(x= log2FoldChange, 
                                   y=-log10(padj),
                                   name=ensgene)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE70462")

g462
ggplotly(g462)
ensemble111 = useEnsembl(biomart = "ensembl", version = 111)
view(listDatasets(ensemble111))
ensemble111 = useDataset("hsapiens_gene_ensembl", mart = ensemble111 )
view(listAttributes(ensemble111))
view(listFilters(ensemble111))
annotation_462 = getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", 
                                      "end_position",
                                      "strand",
                                      "gene_biotype",
                                      "external_gene_name",
                                      "description"),
                       filters =c("ensembl_gene_id"),
                       values = dds462_f1$ensgene,
                       mart = ensemble111)
dim(annotation_462)
dim(dds462_f1)
annotated_462df = left_join(dds462_f1, annotation_462,
                            by = c("ensgene" = "ensembl_gene_id"))
dim(annotated_462df)

g462_annotated= ggplot(annotated_462df, aes(x= log2FoldChange, 
                                            y=-log10(padj),
                                            name= external_gene_name)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE70462")
g462_annotated

ggplotly(g462_annotated)
anno462_df2 = dplyr::filter(annotated_462df, padj <0.05)
anno462_df3 = dplyr::filter(anno462_df2, abs(log2FoldChange) > 1)


## GO analysis

entgene_462 = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = anno462_df3$ensgene,
                    mart = ensemble111)

entgene_462 = as.character(entgene_462$entrezgene_id)


ego_462 = enrichGO(gene = entgene_462,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP")


view(summary(ego_462))
barplot(ego_462, title = "GSE70462")
dotplot(ego_462, title = "GSE70462")

ekg462 = enrichKEGG(gene = entgene_462)

barplot(ekg462, title = "GSE70462")
view(ekg462)
write_tsv(anno462_df3, "DE Genes_GSE70462.txt")
write_tsv(dds462_f3, "DE Genes_filtered_GSE70462.txt")












