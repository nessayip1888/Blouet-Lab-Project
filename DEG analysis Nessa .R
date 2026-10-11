## Below is the code for my bioinformatic analysis using the hindmap dataset (Dowsett et al., 2021)
# Takes a pre-processed Seurat object, performs DEG analysis, then pathway analysis 
#Author: Nessa Yip 2026

#Initialise libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(viridis)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#Load file
filename <- file.choose()
obj <- readRDS(filename)

#tSNE dimensionality reduction 
DimPlot(obj,reduction = "tsne",group.by = "cell_type",label = TRUE)

#Subset cell type (choose from OPCs/NFOLs/OLs etc.)
Idents(obj)<-'cell_type'
subset_oligos<-subset(obj,subset=cell_type == 'OLs')
Idents(subset_oligos) <- "nutr.cond"
markers <- FindMarkers(
  subset_oligos,
  ident.1 = "fast",
  ident.2 = "adlib",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)
markers_filtered <- markers %>%
  rownames_to_column("ensembl_id") %>%
  filter (p_val < 0.05)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = markers_filtered$ensembl_id,
  mart = mart
)
markers_filtered <- markers_filtered %>%
  left_join(gene_map, by = c("ensembl_id" = "ensembl_gene_id"))
markers_filtered <- markers_filtered %>%
  mutate(gene_name = ifelse(mgi_symbol == "", ensembl_id, mgi_symbol))
up=markers_filtered%>%filter(avg_log2FC>=0.5& !is.na(ensembl_id))%>%arrange(desc(avg_log2FC))#set to whatever log2FC you want
down=markers_filtered%>%filter(avg_log2FC<=-0.5& !is.na(ensembl_id))%>%arrange((avg_log2FC))
write_csv(up,'Upregulated in fasted.csv')
write_csv(down,'Downregulated in fasted.csv')

#Put these genes into STRING for an ORA analysis
#Alternatively, do GSEA analysis 