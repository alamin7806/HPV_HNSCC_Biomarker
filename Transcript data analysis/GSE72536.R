library(readr)
library(magrittr)
library(dplyr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(plotly)
library(tibble)
library(biomaRt)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RSkittleBrewer)
library(gplots)
library(AnnotationDbi)
library(RColorBrewer)


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

## MA-plot
plotMA(dds536_res, main = "GSE72536")


dds536_res_df = as.data.frame(dds536_res)
dds536_res_df = rownames_to_column(dds536_res_df, var = "ensgene")
dds536_f1 = dplyr::filter(dds536_res_df, complete.cases(dds536_res_df))
dds536_f2 = dplyr::filter(dds536_f1, padj <0.05)
dds536_f3 = dplyr::filter(dds536_f2, abs(log2FoldChange) > 2)
dds536_f1$test = dds536_f1$padj < 0.05 & 
  abs(dds536_f1$log2FoldChange) > 2
## Volcano plot
g536 = ggplot(dds536_f1, aes(x= log2FoldChange, 
                             y=-log10(padj),
                             name=ensgene)) +
  geom_point(aes(colour = test), size =1.3, alpha = 0.5) +
  geom_vline(xintercept = 2, linetype = 3) +
  geom_vline(xintercept = -2, linetype = 3) +
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
anno536_df3 = dplyr::filter(anno536_df2, abs(log2FoldChange) > 2)

vst536_mat = assay(vst536)

xa = as.matrix(c(rep("HPV Negative", times = 4),
                 rep("HPV Positive", times = 10),
                 rep("HPV Negative", times = 9)))

## Significant Up Regulated Genes

up_genes536 = subset(anno536_df3, log2FoldChange > 0)
head(up_genes536)

top_up_536 <- head(up_genes536[order(up_genes536$log2FoldChange,
                                     decreasing = TRUE),],20)

top_exp536 <- vst536_mat[top_up_536$ensgene,]
pheatmap(top_exp536,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         labels_col = xa,
         col=brewer.pal(name="RdBu", n=11),
         main = "Significant upregulated Genes in GSE72536")




## Significant Down Regulated Genes

down_genes536 <- subset(anno536_df3,log2FoldChange < 0)
head(down_genes536)

top_down536 <- head(down_genes536[order(down_genes536$log2FoldChange, 
                                        decreasing = TRUE), ], 20)
top_down_exp536 <- vst536_mat[top_down536$ensgene,]
pheatmap(top_down_exp536,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         labels_col = xa,
         col=brewer.pal(name="RdBu", n=11),
         main = "Significant Down Regulated Genes in GSE72536")

## GO analysis

entgene_536 = getBM(attributes = c("entrezgene_id"),
                    filters =c("ensembl_gene_id"),
                    values = anno536_df3$ensgene,
                    mart = ensemble111)

entgene_536 = as.character(entgene_536$entrezgene_id)

entUni_536 = getBM(attributes = c("entrezgene_id"),
                   filters =c("ensembl_gene_id"),
                   values = annotated_536df$ensgene,
                   mart = ensemble111)
entUni_536 = as.character(entUni_536$entrezgene_id)

ego_536 = enrichGO(gene = entgene_536,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   universe = entUni_536,
                   readable = TRUE)

ego_536

view(summary(ego_536))
barplot(ego_536, title = "GSE72536")
dotplot(ego_536, title = "GSE72536")

foldchange_536 = anno536_df3$log2FoldChange
names(foldchange_536) = anno536_df3$external_gene_name  

cnetplot(ego_536,showCategory = 10, foldChange = foldchange_536) +
  ggtitle("GSE72536")

ekg536 = enrichKEGG(gene = entgene_536,
                    universe = entUni_536)
view(ekg536)
write_tsv(anno536_df3, "DE Genes_GSE72536.txt")
write_tsv(dds536_f3, "DE Genes_filtered_GSE72536.txt")

















