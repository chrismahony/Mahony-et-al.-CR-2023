---
title: "scRNAseqanalysis"
author: "ChrisMahony"
date: "07/05/2020"
output: html_document
---

```{r}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(EnhancedVolcano)
library(monocle)
options(bitmapType='cairo')
library(CellChat)
setwd("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM")
```




```{r}
wt_2data <- Read10X(data.dir = "/rds/projects/2018/monteirr-lab/2019-12-19.SingleCellRNAseq/WT/outs/filtered_feature_bc_matrix")
wt_2 <- CreateSeuratObject(counts = wt_2data, project = "WT", min.cells = 0, min.features = 0)
mito.features1 <- grep(pattern="^mt-", x=rownames(x=wt_2), value=T)
percent.mito1 <- Matrix::colSums(x = GetAssayData(object = wt_2, slot = "counts")[mito.features1,]) / Matrix::colSums(x = GetAssayData(object = wt_2, slot = "counts"))
wt_2[["percent.mito"]] <- percent.mito1
VlnPlot(object = wt_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
wt_2 <- subset(x = wt_2, subset = nFeature_RNA > 100 & nFeature_RNA <4000 & nCount_RNA > 100 & nCount_RNA < 30000 & percent.mito < 0.2)
wt_2 <- NormalizeData(object = wt_2, verbose = F)
wt_2 <- FindVariableFeatures(object = wt_2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.geneswt_2 <- rownames(wt_2)
wt_2 <- ScaleData(wt_2, features = all.geneswt_2)
```

Load Mutant samples

Then process:

```{r}
hom_2data <- Read10X(data.dir = "/rds/projects/2018/monteirr-lab/2019-12-19.SingleCellRNAseq/mutant/outs/filtered_feature_bc_matrix")
hom_2 <- CreateSeuratObject(counts = hom_2data, project = "MUT", min.cells = 0, min.features = 0)
mito.features2 <- grep(pattern="^mt-", x=rownames(x=hom_2), value=T)
percent.mito2 <- Matrix::colSums(x = GetAssayData(object = hom_2, slot = "counts")[mito.features2,]) / Matrix::colSums(x = GetAssayData(object = hom_2, slot = "counts"))
hom_2[["percent.mito"]] <- percent.mito2
VlnPlot(object = hom_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
hom_2 <- subset(x = hom_2, subset = nFeature_RNA > 100 & nFeature_RNA <4000 & nCount_RNA > 100 & nCount_RNA < 40000 & percent.mito < 0.2)
hom_2 <- NormalizeData(object = hom_2, verbose = F)
hom_2 <- FindVariableFeatures(object = hom_2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.geneshom_2 <- rownames(hom_2)
hom_2 <- ScaleData(hom_2, features = all.geneshom_2)
```


```{r}
reference.list <- c(wt_2,hom_2)
.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated.all.all <- IntegrateData(anchorset = .anchors, dims = 1:30)
DefaultAssay(object=integrated.all.all) <- "integrated"
integrated.all.all <- ScaleData(object = integrated.all.all, verbose=F)
integrated.all.all <- RunPCA(object = integrated.all.all, verbose=F)
ElbowPlot(object = integrated.all.all)
integrated.all.all <- FindNeighbors(object = integrated.all.all, dims = 1:15)
integrated.all.all <- FindClusters(object = integrated.all.all, dim.use= 1:15, resolution = 0.22)
integrated.all.all <- RunUMAP(object = integrated.all.all, reduction = "pca", dims = 1:15)
integrated.all.all <- RunTSNE(object = integrated.all.all, reduction = "pca", dims = 1:15)
p5 <- DimPlot(object = integrated.all.all, reduction = "umap", pt.size=0.5, label = T)
plot_grid(p5)
```


```{r}
DimPlot(object = integrated.all.all, reduction = "tsne", pt.size=0.5, label = T)

```

```{r}

ALLsplitIDs <- read.csv("/rds/projects/m/monteirr-gata2a-resub/figures/tableS7.csv", row.names = 1)
ALLsplitIDs$extra<-'1'
ALLsplitIDs <- ALLsplitIDs[colnames(integrated.all.all@assays$RNA),]
ALLsplitIDs<-as.data.frame(ALLsplitIDs)
length(rownames(ALLsplitIDs))
length(colnames(integrated.all.all))
integrated.all.all2<-integrated.all.all
integrated.all.all2 <- integrated.all.all2[,(colnames(integrated.all.all@assays$RNA) %in% row.names(ALLsplitIDs))]
ALLsplitIDs$extra<-NULL
integrated.all.all2<-AddMetaData(integrated.all.all2, ALLsplitIDs)
DimPlot(object = integrated.all.all2, reduction = "tsne", pt.size=0.5, label = T, group.by = "ALL.split")
table(integrated.all.all2$orig.ident)
```
```{r}

DimPlot(object = integrated.all.all2, reduction = "tsne", pt.size=0.5, label = T, group.by = "orig.ident")
```
```{r}
Idents(integrated.all.all2)<-'ALL.split'
levels(integrated.all.all2) <- c("cHSPCs", "mHSPCs", "Ery-Prog", "Ery", "Macs", "Neutro", "NKT cells", "Bcells", "Multi", "Kidney" )
integrated.all.all2 <- SCTransform(integrated.all.all2, verbose = FALSE)
DefaultAssay(integrated.all.all2)<-'SCT'
all.maerkerszf<-FindAllMarkers(integrated.all.all2, only.pos = T)
all.maerkerszf_top10 <- all.maerkerszf %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
DoHeatmap(integrated.all.all2, features = all.maerkerszf_top10$gene) + NoLegend()

all.maerkerszf_top20 <- all.maerkerszf %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC)

write.csv(all.maerkerszf_top20, "all.maerkerszf_top20.csv")

DoHeatmap(integrated.all.all2, features = all.maerkerszf_top10$gene, group.by = "ident", label = "False") + scale_fill_gradientn(colors = c("black", "yellow"))


```
```{r}
heatmapgns=c("si:dkey-261h17.1", "myb", "lmo2", "cebpa", "gata2a", "gata2b", "meis1b", "runx1", "stmn1a", "dnmt3ba", "angpt1", "gata1a", "alas2", "igfbp1a", "hemgn", "klf1", "c1qa", "mfap4", "mafba", "mpx", "lyz", "vamp8", "nkl.3", "il2rb", "pax5", "cd37", "cd79b", "sox8b", "tekt1", "spink2.2", "viml")

levels(integrated.all.all2) <- c("cHSPCs", "mHSPCs", "Ery-Prog", "Ery", "Macs", "Neutro", "NKT cells", "Bcells", "Multi", "Kidney" )

DefaultAssay(integrated.all.all2)<-'RNA'
integrated.all.all2<-NormalizeData(integrated.all.all2)
integrated.all.all2 <- ScaleData(integrated.all.all2, verbose = FALSE)
integrated.all.all2<-FindVariableFeatures(integrated.all.all2)
integrated.all.all2 <- RunPCA(integrated.all.all2, npcs = 30, verbose = FALSE)
DoHeatmap(integrated.all.all2, features = heatmapgns, group.by = "ident", label = "False") + scale_fill_gradientn(colors = c("black", "yellow"))

```
```{r}
DotPlot(integrated.all.all2, features = heatmapgns) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.7) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```


```{r}
levels(integrated.all.all2)
DimPlot(integrated.all.all2, label = T, cols= c('Macs'='purple', 'Bcells'='#00BE67', 'Neutro'='pink', 'Ery-Prog'='#8494FF', 'Ery'='red', 'NKT cells'='#00B8E7', 'Kidney'='#7CAE00', 'cHSPCs'='#00BFC4', 'mHSPCs'='orange', 'Multi'='brown'), pt.size = 0.5, reduction = "tsne" )


```
```{r}
DimPlot(integrated.all.all2, group.by = "orig.ident", reduction="tsne", pt.size = 0.5)
```
```{r}
#generate tSNE plots in loupe to keep the same colours
tsne_cord<-integrated.all.all2@reductions[["tsne"]]@cell.embeddings
write.csv(tsne_cord, "tsne_cord.csv")
rm(tsne_cord)


labels<-integrated.all.all2@meta.data
labels <- subset(labels, select = c(7))
write.csv(labels, "labels.csv")
rm(labels)

```


```{r}
DefaultAssay(integrated.all.all2)<-'RNA'

x <- VlnPlot(integrated.all.all2, features = c("si:dkey-261h17.1"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)


x <- VlnPlot(integrated.all.all2, features = c("gata2a"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)

x <- VlnPlot(integrated.all.all2, features = c("meis1b"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)


x <- VlnPlot(integrated.all.all2, features = c("spi1b"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)

x <- VlnPlot(integrated.all.all2, features = c("cebpa"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)

x <- VlnPlot(integrated.all.all2, features = c("hemgn"), ncol = 3, idents = c("cHSPCs", "mHSPCs", "Ery-Prog"), pt.size = 0) +  scale_x_discrete(limits = c("cHSPCs", "mHSPCs", "Ery-Prog"))
cols <- c("cHSPCs" = "#00BFC4", "mHSPCs" = "orange", "Ery-Prog" = "#8494FF")
x + scale_fill_manual(values = cols)

```



```{r}
integrated.all.all2$cluster_ID<-paste(integrated.all.all2$ALL.split, integrated.all.all2$orig.ident, sep = "_")
DefaultAssay(integrated.all.all2)<-'SCT'
Idents(integrated.all.all2) <- 'cluster_ID'
DE_mHSPCS_SCT <-FindMarkers(integrated.all.all2, ident.1 = "mHSPCs_MUT", ident.2 ="mHSPCs_WT")
DE_mHSPCS_SCT$gene<-rownames(DE_mHSPCS_SCT)
DE_Neutr <-FindMarkers(integrated.all.all2, ident.1 = "Neutro_MUT", ident.2 ="Neutro_WT")
DE_Neutr$gene<-rownames(DE_Neutr)
DE_Bcells <-FindMarkers(integrated.all.all2, ident.1 = "Bcells_MUT", ident.2 ="Bcells_WT")
DE_Bcells$gene<-rownames(DE_Bcells)
DE_macs <-FindMarkers(integrated.all.all2, ident.1 = "Macs_MUT", ident.2 ="Macs_WT")
DE_macs$gene<-rownames(DE_macs)
DE_ery <-FindMarkers(integrated.all.all2, ident.1 = "Ery_MUT", ident.2 ="Ery_WT")
DE_ery$gene<-rownames(DE_ery)
DE_NKTcells <-FindMarkers(integrated.all.all2, ident.1 = "NKT cells_MUT", ident.2 ="NKT cells_WT")
DE_NKTcells$gene<-rownames(DE_NKTcells)
DE_cHSPCs <-FindMarkers(integrated.all.all2, ident.1 = "cHSPCs_MUT", ident.2 ="cHSPCs_WT")
DE_cHSPCs$gene<-rownames(DE_cHSPCs)
DE_Ery_Prog <-FindMarkers(integrated.all.all2, ident.1 = "Ery-Prog_MUT", ident.2 ="Ery-Prog_WT" )
DE_Ery_Prog$gene<-rownames(DE_Ery_Prog)

DE_mHSPCS_SCT$cluster<-'mHSPCs'
DE_Neutr$cluster<-'Neutrophils'
DE_Bcells$cluster<-'Bcells'
DE_macs$cluster<-'Macrophages'
DE_ery$cluster<-'Erythrocytes'
DE_NKTcells$cluster<-'NK T cells'
DE_cHSPCs$cluster<-'cHSPCs'
DE_Ery_Prog$cluster<-'Erythrocyte progenitors'

all_DE<-rbind(DE_mHSPCS_SCT, DE_Neutr,  DE_Bcells, DE_macs, DE_ery, DE_NKTcells, DE_cHSPCs, DE_Ery_Prog)
write.csv(all_DE, "all_DE.csv")

```
```{r}
EnhancedVolcano(DE_mHSPCS_SCT,
    lab = DE_mHSPCS_SCT$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c('npm1a','mcm3','cebpa','lcp1', 'hells','hbba1.1','hbaa2', 'hemgn', 'ncl', 'lyz', 'spi1b'),
    title = 'DE_mHSPCS_SCT',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    max.overlaps = 300,
 
    )
```
```{r}
#DE_HSCgenes_from_paper <- read_excel("DE_HSCgenes_from_paper.xlsx") #DE otholougous gene from Wu et al., 2020, supplemental table 5.
DE_mHSPCS_SCT_up_from_paper <- DE_mHSPCS_SCT[(rownames(DE_mHSPCS_SCT) %in% DE_HSCgenes_from_paper$up_reg),]
DE_mHSPCS_SCT_up_from_paper <- DE_mHSPCS_SCT_up_from_paper[DE_mHSPCS_SCT_up_from_paper$avg_log2FC > 0.25, ]
DE_mHSPCS_SCT_up_from_paper <- DE_mHSPCS_SCT_up_from_paper[DE_mHSPCS_SCT_up_from_paper$p_val < 0.05, ]

DE_mHSPCS_SCT_down_from_paper <- DE_mHSPCS_SCT[(rownames(DE_mHSPCS_SCT) %in% DE_HSCgenes_from_paper$down_reg),]
DE_mHSPCS_SCT_down_from_paper <- DE_mHSPCS_SCT_down_from_paper[DE_mHSPCS_SCT_down_from_paper$avg_log2FC < -0.25, ]
DE_mHSPCS_SCT_down_from_paper <- DE_mHSPCS_SCT_down_from_paper[DE_mHSPCS_SCT_down_from_paper$p_val < 0.05, ]

all_from_paper<-rbind(DE_mHSPCS_SCT_up_from_paper, DE_mHSPCS_SCT_down_from_paper)



EnhancedVolcano(DE_mHSPCS_SCT,
    lab = DE_mHSPCS_SCT$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = rownames(all_from_paper),
    title = 'DE_mHSPCS_SCT',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    max.overlaps = 300,
 
    )
```
```{r}

DE_cHSPCs_SCT_up_from_paper <- DE_cHSPCs[(rownames(DE_cHSPCs) %in% DE_HSCgenes_from_paper$up_reg),]
DE_cHSPCs_SCT_up_from_paper <- DE_cHSPCs_SCT_up_from_paper[DE_cHSPCs_SCT_up_from_paper$avg_log2FC > 0.25, ]
DE_cHSPCs_SCT_up_from_paper <- DE_cHSPCs_SCT_up_from_paper[DE_cHSPCs_SCT_up_from_paper$p_val < 0.05, ]

DE_cHSPCs_SCT_down_from_paper <- DE_cHSPCs[(rownames(DE_cHSPCs) %in% DE_HSCgenes_from_paper$down_reg),]
DE_cHSPCs_SCT_down_from_paper <- DE_cHSPCs_SCT_down_from_paper[DE_cHSPCs_SCT_down_from_paper$avg_log2FC < -0.25, ]
DE_cHSPCs_SCT_down_from_paper <- DE_cHSPCs_SCT_down_from_paper[DE_cHSPCs_SCT_down_from_paper$p_val < 0.05, ]

all_from_paper_cHSPCs<-rbind(DE_cHSPCs_SCT_up_from_paper, DE_cHSPCs_SCT_down_from_paper)

EnhancedVolcano(DE_cHSPCs,
    lab = DE_cHSPCs$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = rownames(all_from_paper_cHSPCs),
    title = 'DE_cHSPCs',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    max.overlaps = 300,
 
    )

```


```{r}
EnhancedVolcano(DE_cHSPCs,
    lab = DE_cHSPCs$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("atp6v0e1"),
    title = 'DE_cHSPCs',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )



```
```{r}
#sense checking
mHSPCs_subset <- integrated.all.all2[,grepl("mHSPCs", integrated.all.all2$ALL.split, ignore.case=TRUE)]
Idents(mHSPCs_subset)<-'orig.ident'
levels(mHSPCs_subset)<-c('WT', 'MUT')
VlnPlot(mHSPCs_subset, features = c("cebpa", "lcp1", "npm1a"), pt.size = 0)
VlnPlot(mHSPCs_subset, features = c("cebpa", "hbba1.1", "hemgn"), pt.size = 0)
```
```{r}
#sense checking
cHSPCs_subset <- integrated.all.all2[,grepl("cHSPCs", integrated.all.all2$ALL.split, ignore.case=TRUE)]
Idents(cHSPCs_subset)<-'orig.ident'
levels(cHSPCs_subset)<-c('WT', 'MUT')
VlnPlot(cHSPCs_subset, features = c("cebpa", "lcp1", "npm1a"), pt.size = 0.0)
VlnPlot(cHSPCs_subset, features = c("cebpa", "hbba1.1", "hemgn"), pt.size = 0.0)
```
```{r}

EnhancedVolcano(DE_ery,
    lab = DE_ery$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("lyz", "hbba2", "hbaa2", "hbaa1", "lect2l", "alas2"),
    title = 'Ery',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )


```
```{r}
EnhancedVolcano(DE_Ery_Prog,
    lab = DE_Ery_Prog$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("si:ch211-5k11.8", "hbba1.1", "hbba2", "lyz", "hemgn", "hbaa1"),
    title = 'Ery_prog',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )


```
```{r}
EnhancedVolcano(DE_Bcells,
    lab = DE_Bcells$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("igic1s1", "igl4v8", "igl1c3", "hbba1", "hbaa1", "hbba2", "hbaa2", "lcp1"),
    title = 'DE_Bcells',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )
```
```{r}
EnhancedVolcano(DE_macs,
    lab = DE_macs$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("atp6v0e1"),
    title = 'DE_macs',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )
```
```{r}
EnhancedVolcano(DE_Neutr,
    lab = DE_Neutr$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("lcp1"),
    title = 'DE_Neutr',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )
```
```{r}
EnhancedVolcano(DE_NKTcells,
    lab = DE_NKTcells$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c("si:ch211-5k11.8"),
    title = 'DE_NKTcells',
    subtitle = "RNA, red=p<0.05 & FC > 0.25",
     pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.5,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5,
    
    )
```
```{r}

#data to write out for paper
#cell_embedding<-integrated.all.all2@active.ident
#write.csv(cell_embedding, "cell_embeddings_fig2.csv")
#rm(cell_embedding)


```





