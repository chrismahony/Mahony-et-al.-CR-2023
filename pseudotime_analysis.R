```{r}
library("Seurat")
library("stringr")
library("ggplot2")
library(patchwork)
library(dplyr)
options(bitmapType='cairo')
```

```{r}


load("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/analysis.RData")
rm(list=ls()[! ls() %in% c("integrated.all.all2")])

# Read in data
UMAP_coordinatesR <- read.csv("/rds/projects/2018/monteirr-lab/2019-12-19.SingleCellRNAseq/CM_cloupe_analysis/EGFPpa_CustomfilterDefault/projectionUMAPnew2.csv", stringsAsFactors = F, header = T, row.names = 1)

#edit rownames in excel
write.csv(UMAP_coordinatesR, "UMAP_coordinatesR.csv")
UMAP_coordinatesR <- read.csv("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/UMAP_coordinatesR.csv", stringsAsFactors = F, header = T, row.names = 1)

#filter to remove extra cells
UMAP_coordinates2 <- UMAP_coordinatesR[colnames(integrated.all.all2@assays$RNA),]


# conver to matrix
UMAP_coordinates_mat_new <- as(UMAP_coordinates2, "matrix")


# Create DimReducObject and add to object
integrated.all.all2[['UMAP']] <- CreateDimReducObject(embeddings = UMAP_coordinates_mat_new, key = "UMAP_", global = T, assay = "RNA")

DimPlot(integrated.all.all2, reduction = "UMAP", cols= c('Macs'='purple', 'Bcells'='#00BE67', 'Neutro'='pink', 'Ery-Prog'='#8494FF', 'Ery'='red', 'NKT cells'='#00B8E7', 'Kidney'='#7CAE00', 'cHSPCs'='#00BFC4', 'mHSPCs'='orange', 'Multi'='brown'), pt.size = 0.3)

```
```{r}
library(SeuratWrappers)
library(monocle3)
srt=integrated.all.all2
DefaultAssay(srt)<-'RNA'
srt.cds <- as.cell_data_set(srt)
srt.cds <- cluster_cells(cds = srt.cds, reduction_method = "UMAP")
srt.cds <- learn_graph(srt.cds, use_partition = TRUE)
rownames(srt.cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(srt.cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
plot_cells(srt.cds)
```
```{r}
srt.cds <- order_cells(srt.cds)
```
```{r}
srt.cds <- estimate_size_factors(srt.cds)
srt.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srt.cds[["RNA"]])
rowData(srt.cds)$gene_short_name <- row.names(rowData(srt.cds))
```

```{r}
plot_cells( cds = srt.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = F,
  label_roots = F,
  cell_size = 0.7
)

plot_cells( cds = srt.cds,
  color_cells_by = "ALL.split",
  show_trajectory_graph = F,
  label_roots = F,
  label_cell_groups = F,
  cell_size = 0.7
)
```
```{r}

WT_srt.cds <- srt.cds[,grepl("WT", colData(srt.cds)$orig.ident, ignore.case=TRUE)]
Mut_srt.cds <- srt.cds[,grepl("MUT", colData(srt.cds)$orig.ident, ignore.case=TRUE)]


WT_genes <- WT_srt.cds[row.names(subset(rowData(WT_srt.cds),
                 gene_short_name %in% c("alas2", "cebpa", "mpx", "lcp1"))),]

Mut_genes <- Mut_srt.cds[row.names(subset(rowData(Mut_srt.cds),
                 gene_short_name %in% c("alas2", "cebpa", "mpx", "lcp1"))),]


p1 <- plot_genes_in_pseudotime(WT_genes,
                         color_cells_by="pseudotime",  
                         min_expr=0.5)

p2 <- plot_genes_in_pseudotime(Mut_genes,
                         color_cells_by="pseudotime", 
                         min_expr=0.5)

p1 + p2

p1 <- p1 + scale_y_continuous(limits = c(0.5, 4))
p2 <- p2 + scale_y_continuous(limits = c(0.5, 4))
p1 + p2



```


