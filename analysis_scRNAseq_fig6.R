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

wt_2data <- Read10X(data.dir = "/rds/projects/2018/monteirr-lab/2019-06-11.SingleCellRNAseq/WT/outs/filtered_feature_bc_matrix")
wt_2 <- CreateSeuratObject(counts = wt_2data, project = "WT", min.cells = 0, min.features = 0)
wt_2 <- NormalizeData(object = wt_2, verbose = F)
wt_2 <- FindVariableFeatures(object = wt_2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.geneswt_2 <- rownames(wt_2)
wt_2 <- ScaleData(wt_2, features = all.geneswt_2)

hom_2data <- Read10X(data.dir = "/rds/projects/2018/monteirr-lab/2019-06-11.SingleCellRNAseq/mutant/outs/filtered_feature_bc_matrix")
hom_2 <- CreateSeuratObject(counts = hom_2data, project = "MUT", min.cells = 0, min.features = 0)
hom_2 <- NormalizeData(object = hom_2, verbose = F)
hom_2 <- FindVariableFeatures(object = hom_2, selection.method = "vst", nfeatures = 2000, verbose=F)
all.geneshom_2 <- rownames(hom_2)
hom_2 <- ScaleData(hom_2, features = all.geneshom_2)

```
```{r}
reference.list <- c(wt_2,hom_2)
.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated.all.all <- IntegrateData(anchorset = .anchors, dims = 1:30)
```
```{r}

min(integrated.all.all$nCount_RNA)
VlnPlot(integrated.all.all, "nCount_RNA")

```


```{r}
tail(colnames(integrated.all.all))
tableS9 <- read_csv("/rds/projects/m/monteirr-gata2a-resub/figures/tableS9.csv")
barcodes_to_retain_CM=tableS9
ncol(integrated.all.all)
integrated.all.all <- subset(integrated.all.all, cells = barcodes_to_retain_CM$Barcode)
integrated.all.all <- SCTransform(integrated.all.all, verbose = FALSE)
DefaultAssay(integrated.all.all)<-'SCT'
Idents(integrated.all.all)<-'orig.ident'
DE_genes<-FindMarkers(integrated.all.all, ident.1 = "MUT", ident.2 = "WT")

EnhancedVolcano(DE_genes,
    lab = rownames(DE_genes),
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = c('npm1a','lyz','lect2l','hbba2', 'hbaa2'),
    title = 'DE_genes',
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
levels(integrated.all.all)<-c("MUT", "WT")
VlnPlot(integrated.all.all, features =  c('npm1a','lyz','lect2l','hbba2', 'hbaa2'), stack = T)+ geom_jitter()





```
```{r}
tSNE_fig6 <- read.csv("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/tSNE_fig6.csv", stringsAsFactors = F, header = T, row.names = 1)
tSNE_fig6 <- as(tSNE_fig6, "matrix")
integrated.all.all[['tSNE']] <- CreateDimReducObject(embeddings = tSNE_fig6, key = "tSNE_", global = T, assay = "RNA")
DimPlot(integrated.all.all, reduction = "tSNE")
```

