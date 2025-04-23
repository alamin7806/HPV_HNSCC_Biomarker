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


view(read_csv("SraRunTable_GSE250305.txt"))
table_305 = read_csv("SraRunTable_GSE250305.txt") %>%
  dplyr::select(Run, BioProject, BioSample, cell_line, cell_type, Experiment,
                genotype, `Library Name`,`Sample Name`, source_name,
                `SRA Study`)
view(table_305)

colnames(table_305)[c(8,9,11)] = c("Library_name","Sample_name", "SRA_Study")


sample_file_305 = paste0(pull(table_305, Run), "/quant.sf")
sample_file_305
names(sample_file_305) = pull(table_305, Sample_name)
sample_file_305
ref_gene = read_csv("gene_id.txt", col_names = c("enst", "ensg"))

## from transcript to gene level count
count_data_305 = tximport(files = sample_file_305,
                          type = "salmon",
                          tx2gene = ref_gene,
                          ignoreTxVersion = TRUE)
table_305 = as.data.frame(table_305)

condition = factor(c(rep("HPV_Positive",9), rep("HPV_Negative",12)))
condition
table_305$condition = condition
view(table_305)

## construct a DESeq DataSet from count data

deseq_data305 = DESeqDataSetFromTximport(txi = count_data_305,
                                         colData = table_305,
                                         design = ~condition)

## pre-filtering data

deseq_data305_norm = rowSums(counts(deseq_data305)) >=10
deseq_data305_norm = deseq_data305[deseq_data305_norm,]

## Differential expression analysis

dds305 = DESeq(deseq_data305_norm)
vst305 = varianceStabilizingTransformation(dds305)
plotPCA(vst305, intgroup = "condition") + ggtitle("GSE250305")
dds305_res = results(dds305)
summary(dds305_res)

dds305_res_df = as.data.frame(dds305_res)
dds305_res_df = rownames_to_column(dds305_res_df, var = "ensgene")
dds305_f1 = dplyr::filter(dds305_res_df, complete.cases(dds305_res_df))
dds305_f2 = dplyr::filter(dds305_f1, padj <0.05)
dds305_f3 = dplyr::filter(dds305_f2, abs(log2FoldChange) > 1)
dds305_f1$test = dds305_f1$padj < 0.05 & abs(dds305_f1$log2FoldChange) > 1
g305 = ggplot(dds305_f1, aes(x= log2FoldChange, 
                             y=-log10(padj),
                             name=ensgene)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE250305")

g305
ggplotly(g305)

ensemble111 = useEnsembl(biomart = "ensembl", version = 111)
view(listDatasets(ensemble111))
ensemble111 = useDataset("hsapiens_gene_ensembl", mart = ensemble111 )
view(listAttributes(ensemble111))
view(listFilters(ensemble111))
annotation_305 = getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", 
                                      "end_position",
                                      "strand",
                                      "gene_biotype",
                                      "external_gene_name",
                                      "description"),
                       filters =c("ensembl_gene_id"),
                       values = dds305_f1$ensgene,
                       mart = ensemble111)
dim(annotation_305)
dim(dds305_f1)
annotated_305df = left_join(dds305_f1, annotation_305,
                            by = c("ensgene" = "ensembl_gene_id"))
dim(annotated_305df)

g305_annotated= ggplot(annotated_305df, aes(x= log2FoldChange, 
                                            y=-log10(padj),
                                            name= external_gene_name)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE250305")

g305_annotated
ggplotly(g305_annotated)

anno305_df2 = dplyr::filter(annotated_305df, padj <0.05)
anno305_df3 = dplyr::filter(anno305_df2, abs(log2FoldChange) > 1)



## GO analysis

entgene_305 = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = anno305_df3$ensgene,
                    mart = ensemble111)

entgene_305 = as.character(entgene_305$entrezgene_id)



ego_305 = enrichGO(gene = entgene_305,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP")

ego_305
view(summary(ego_305))

barplot(ego_305, title = "GSE250305")
dotplot(ego_305, title = "GSE250305")


ekg305 = enrichKEGG(gene = entgene_305)

dotplot(ekg305, title = "GSE250305")

view(ekg305)
write_tsv(anno305_df3, "DE Genes_GSE250305.txt")
write_tsv(dds305_f3, "DE Genes_filtered_GSE250305.txt")

























































































































