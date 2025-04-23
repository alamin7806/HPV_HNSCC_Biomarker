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



view(read_csv("SraRunTable_536.txt"))
table_536 = read_csv("SraRunTable_536.txt") %>%
  dplyr::select(Run, BioProject, BioSample, Experiment,
                `GEO_Accession (exp)`,hpv_status,`Sample Name`, source_name,
                til_status,tissue, tumour_site, tumour_stage)
view(table_536)

colnames(table_536)[c(5,7)] = c("GEO_Accession","Sample_name")
table_536 = slice(table_536, seq(1,46, by =2))

sample_file_536 = paste0(pull(table_536, Sample_name), "/quant.sf")
sample_file_536
names(sample_file_536) = pull(table_536, Sample_name)
sample_file_536
ref_gene = read_csv("gene_id.txt", col_names = c("enst", "ensg"))

## from transcript to gene level count

count_data_536 = tximport(files = sample_file_536,
                          type = "salmon",
                          tx2gene = ref_gene,
                          ignoreTxVersion = TRUE)
table_536 = as.data.frame(table_536)

condition = factor(table_536$hpv_status)
condition
table_536$condition = condition


## construct a DESeq DataSet from count data
deseq_data536 = DESeqDataSetFromTximport(txi = count_data_536,
                                         colData = table_536,
                                         design = ~condition)

##pre-filtering data

deseq_data536_norm = rowSums(counts(deseq_data536)) >=10
deseq_data536_norm = deseq_data536[deseq_data536_norm,]

## Differential expression analysis

dds536 = DESeq(deseq_data536_norm)
vst536 = varianceStabilizingTransformation(dds536)
plotPCA(vst536, intgroup = "condition") + ggtitle("GSE772536")
dds536_res = results(dds536)
summary(dds536_res)


dds536_res_df = as.data.frame(dds536_res)
dds536_res_df = rownames_to_column(dds536_res_df, var = "ensgene")
dds536_f1 = dplyr::filter(dds536_res_df, complete.cases(dds536_res_df))
dds536_f2 = dplyr::filter(dds536_f1, padj <0.05)
dds536_f3 = dplyr::filter(dds536_f2, abs(log2FoldChange) > 1)
dds536_f1$test = dds536_f1$padj < 0.05 & 
  abs(dds536_f1$log2FoldChange) > 1
## Volcano plot
g536 = ggplot(dds536_f1, aes(x= log2FoldChange, 
                             y=-log10(padj),
                             name=ensgene)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE2536")
g536

ggplotly(g536)
ensemble111 = useEnsembl(biomart = "ensembl", version = 111)
view(listDatasets(ensemble111))
ensemble111 = useDataset("hsapiens_gene_ensembl", mart = ensemble111 )
view(listAttributes(ensemble111))
view(listFilters(ensemble111))
annotation_536 = getBM(attributes = c("ensembl_gene_id", 
                                      "chromosome_name",
                                      "start_position", 
                                      "end_position",
                                      "strand",
                                      "gene_biotype",
                                      "external_gene_name",
                                      "description"),
                       filters =c("ensembl_gene_id"),
                       values = dds536_f1$ensgene,
                       mart = ensemble111)
dim(annotation_536)
dim(dds536_f1)

annotated_536df = left_join(dds536_f1, annotation_536,
                            by = c("ensgene" = "ensembl_gene_id"))
dim(annotated_536df)

g536_annotated= ggplot(annotated_536df, aes(x= log2FoldChange, 
                                            y=-log10(padj),
                                            name= external_gene_name)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_vline(xintercept = -1, linetype = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Volcano Plot of GSE70462")
g536_annotated

ggplotly(g536_annotated)
anno536_df2 = dplyr::filter(annotated_536df, padj <0.05)
anno536_df3 = dplyr::filter(anno536_df2, abs(log2FoldChange) > 1)


## GO analysis

entgene_536 = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = anno536_df3$ensgene,
                    mart = ensemble111)

entgene_536 = as.character(entgene_536$entrezgene_id)


ego_536 = enrichGO(gene = entgene_536,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP")

ego_536

view(summary(ego_536))
barplot(ego_536, title = "GSE72536")
dotplot(ego_536, title = "GSE72536")


ekg536 = enrichKEGG(gene = entgene_536)

barplot(ekg536, title = "GSE72536")
view(ekg536)
write_tsv(anno536_df3, "DE Genes_GSE72536.txt")
write_tsv(dds536_f3, "DE Genes_filtered_GSE72536.txt")

















