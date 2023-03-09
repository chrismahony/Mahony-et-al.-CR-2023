

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

data <- read.table("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/GSE146404_zebrafish_counts_change.3stage.merge.txt", quote="\"", comment.char="")
paper_data <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
paper_data <- NormalizeData(object = paper_data, verbose = F)
paper_data <- FindVariableFeatures(object = paper_data, selection.method = "vst", nfeatures = 2000, verbose=F)
all.genes <- rownames(paper_data)
paper_data <- ScaleData(paper_data, features = all.genes)
```
```{r}
GSM340206metadata <-read.csv("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/GSM340206metadata.csv", row.names=1)
paper_data<-AddMetaData(paper_data, GSM340206metadata)
```

```{r}
load("/rds/projects/m/monteirr-gata2a-resub/analysis_revision_CM/analysis.RData")
rm(list=ls()[! ls() %in% c("all.genes","data", "GSM340206metadata", "paper_data", "all.maerkerszf", "integrated.all.all2", "DE_mHSPCS_SCT")])
```

```{r}
all.maerkerszf_mHSPCs<-all.maerkerszf[all.maerkerszf$cluster=="mHSPCs",]
all.maerkerszf_mHSPCs <- all.maerkerszf_mHSPCs[all.maerkerszf_mHSPCs$p_val < 0.05, ]
Idents(paper_data)<-'celltype'
paper_data<- SCTransform(paper_data)
DefaultAssay(paper_data)<-'SCT'
allmarkers_paper<-FindAllMarkers(paper_data, only.pos = T)
```
```{r}
allmarkers_paper_HSPC1<-allmarkers_paper[allmarkers_paper$cluster=="HSPC1",]
allmarkers_paper_HSPC2<-allmarkers_paper[allmarkers_paper$cluster=="HSPC2",]
allmarkers_paper_HSPC3<-allmarkers_paper[allmarkers_paper$cluster=="HSPC3",]
allmarkers_paper_HSPC4<-allmarkers_paper[allmarkers_paper$cluster=="HSPC4",]

allmarkers_paper_HSPC1 <- allmarkers_paper_HSPC1[allmarkers_paper_HSPC1$p_val < 0.05, ]
allmarkers_paper_HSPC2 <- allmarkers_paper_HSPC2[allmarkers_paper_HSPC2$p_val < 0.05, ]
allmarkers_paper_HSPC3 <- allmarkers_paper_HSPC3[allmarkers_paper_HSPC3$p_val < 0.05, ]
allmarkers_paper_HSPC4 <- allmarkers_paper_HSPC4[allmarkers_paper_HSPC4$p_val < 0.05, ]

HSPC1_mHSPCoverlap <- all.maerkerszf_mHSPCs[all.maerkerszf_mHSPCs$gene %in% allmarkers_paper_HSPC1$gene,]
HSPC2_mHSPCoverlap <- all.maerkerszf_mHSPCs[all.maerkerszf_mHSPCs$gene %in% allmarkers_paper_HSPC2$gene,]
HSPC3_mHSPCoverlap <- all.maerkerszf_mHSPCs[all.maerkerszf_mHSPCs$gene %in% allmarkers_paper_HSPC3$gene,]
HSPC4_mHSPCoverlap <- all.maerkerszf_mHSPCs[all.maerkerszf_mHSPCs$gene %in% allmarkers_paper_HSPC4$gene,]
#goi<-c("tcp1", "npm1a", "pcna", "mcm3", "hells", "ncl")

HSPC3_mHSPCoverlap_DEmHSPCs <- DE_mHSPCS_SCT[DE_mHSPCS_SCT$gene %in% HSPC3_mHSPCoverlap$gene,]



write.csv(HSPC3_mHSPCoverlap_DEmHSPCs, "HSPC3_mHSPCoverlap_DEmHSPCs.csv")

EnhancedVolcano(HSPC3_mHSPCoverlap_DEmHSPCs,
    lab = HSPC3_mHSPCoverlap_DEmHSPCs$gene,
        x = 'avg_log2FC',
    y = 'p_val',
     selectLab = "mt-nd1",
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

