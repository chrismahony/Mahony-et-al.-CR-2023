```{r}


library(dplyr)
library(Signac)
library(Seurat)
library("AnnotationHub")
library(ggplot2)
library(patchwork)
#library(JASPAR2020)
#library(TFBSTools)
library(patchwork)
set.seed(1986)
options(bitmapType='cairo')
```


```{r}

# Load the peak/cell matrix
counts <- Read10X_h5(filename = "/rds/projects/m/mahonyc-gata2-mutated-mds-aml/scATACseqAugust2020_analysis/Chrystala_analysis/Files/Files/filtered_peak_bc_matrix.h5")
                    

# Load the metadata                     
metadata <- read.csv(
  file = "/rds/projects/m/mahonyc-gata2-mutated-mds-aml/scATACseqAugust2020_analysis/Chrystala_analysis/Files/Files/singlecell.csv",
  header = TRUE,
  row.names = 1)


# Create chromatin assay using the fragments file
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
   fragments = "/rds/projects/m/mahonyc-gata2-mutated-mds-aml/scATACseqAugust2020_analysis/Chrystala_analysis/Files/Files/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200)


# Create Seurat object with the chromatin assay
aggr <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata,
)


ah = AnnotationHub()
ah

ahDb <- query(ah, pattern = c("Danio rerio", "EnsDb", 99))
ahDb

# Retrieve Danio Rerio from query output
ahEdb <- ahDb[[1]]

# Retrieve all genes
gns <- genes(ahEdb)

# Retrieve genomic ranges
annotations <- GetGRangesFromEnsDb(ensdb = ahEdb)

# add the gene information to the object
genome(annotations) <- "danRer11"
Annotation(aggr) <- annotations



aggr
```



```{r}
aggr[['peaks']]
granges(aggr)
```

```{r}


save.image("aggr.RData")

load.image("aggr.RData")

```

```{r}
# compute nucleosome signal score per cell
aggr <- NucleosomeSignal(object = aggr)

# compute TSS enrichment score per cell  
aggr <- TSSEnrichment(object = aggr, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
aggr$pct_reads_in_peaks <- aggr$peak_region_fragments / aggr$passed_filters * 100
aggr$blacklist_ratio <- aggr$blacklist_region_fragments / aggr$peak_region_fragments

aggr$high.tss <- ifelse(aggr$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(aggr, group.by = 'high.tss') + NoLegend()


#save.image("aggr.RData")

```

```{r}
aggr$nucleosome_group <- ifelse(aggr$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = aggr, group.by = 'nucleosome_group')
```

```{r}
VlnPlot(
  object = aggr,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

```

```{r}
VlnPlot(
  object = aggr,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2,
  group.by = 'sample'
)

```

```{r}
#Based on metrics above remove cells that are outliers.
aggr <- subset(
  x = aggr,
  subset = peak_region_fragments < 6000 &
    pct_reads_in_peaks > 40 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 2 
)
aggr

save.image("aggr.RData")
```

```{r}
aggr <- RunTFIDF(aggr)
aggr <- FindTopFeatures(aggr, min.cutoff = 'q0')
aggr <- RunSVD(aggr)
DepthCor(aggr)
```



```{r}
DefaultAssay(aggr) <- 'peaks'

aggr <- RunUMAP(object = aggr, reduction = 'lsi', dims = 2:30)
aggr <- FindNeighbors(object = aggr, reduction = 'lsi', dims = 2:30)
aggr <- FindClusters(object = aggr, verbose = FALSE, algorithm = 3, resolution = 0.1)
DimPlot(object = aggr, label = TRUE) + NoLegend()

aggr <- RenameIdents(
  object = aggr,
  'neutr' = 'neutr',
  'WTery' = 'ery',
  'Mutery' = 'ery',
  'macs' = 'macs',
  'Bcells' = 'Bcells',
  'multi' = 'multi',
  'NKTcells' = 'NKTcells',
  'kidney' = 'kidney',
  'HSPCs' = 'HSPCs',
  'DN1' = 'DN1',
  'DN2' = 'DN2',
  'endo' = 'endo'
  )

aggr$cluster_named<-aggr@active.ident

DimPlot(object = aggr, label = TRUE,group.by = "cluster_named") + NoLegend()



write.table(closest_genes_DN2_closed_mutant, file = "table.txt")

```


```{r}

# Distribution of libraries accross the clusters
aggr$orig.ident<-aggr$library_ID
current.sample.ids <- c("MUT1","MUT2","MUT3", "WT", "WT1", "WT2", "WT3")
new.sample.ids <- c("MUT","MUT","WT", "WT", "WT", "WT", "MUT")
aggr@meta.data[["library_ID"]] <- plyr::mapvalues(x = aggr@meta.data[["library_ID"]], from = current.sample.ids, to = new.sample.ids)

aggr <- RunTSNE(object = aggr, reduction = 'lsi', dims = 2:30, spread= 2, min.dist= 0.8)
DimPlot(object = aggr, label = TRUE, reduction = "tsne")
DimPlot(object = aggr, group.by = "library_ID", label = TRUE, reduction = "tsne")

Idents(aggr)<-'cluster_named'
DimPlot(object = aggr, label = T, cols= c('macs'='purple', 'Bcells'='#00BE67', 'neutr'='pink', 'HSPCs'='#8494FF', 'ery'='red', 'NKTcells'='#00B8E7', 'endo'='grey', 'kidney'='#7CAE00', 'DN2'='#00BFC4', 'DN1'='grey', 'dHSPCs'='orange', 'multi'='brown'), pt.size = 0.5, reduction = "tsne" )

DimPlot(object = aggr, label = FALSE, group.by = 'library_ID', cols= c('WT'='red', 'MUT'='#0CB702'), pt.size = 2, reduction = "tsne" )

dev.off()


```

```{r}

#GeneACtivity had to be submitted through slurm for me to work!
gene.activities <- GeneActivity(aggr)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
aggr[['RNA']] <- CreateAssayObject(counts = gene.activities)
aggr <- NormalizeData(
  object = aggr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(aggr$nCount_RNA)
)

```


```{r}

DefaultAssay(aggr) <- 'RNA'

all.genes <- rownames(aggr)
aggr  <- ScaleData(aggr, features = all.genes)

aggr.markers <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)
top5RNA <- aggr.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(aggr, features = plotgns, group.by = "ident") + scale_fill_gradientn(colors = c("black", "yellow"))

#DN2=CHSPCs, HSPCs=ery-P
levels(aggr) <- c("DN2", "dHSPCs", "HSPCs", "ery", "macs", "neutr", "NKTcells", "Bcells", "multi", "kidney", "endo", "DN1" )
levels(aggr)
          
saveRDS(aggr, file = "aggr2.rds")

```



```{r}

DefaultAssay(aggr) <- 'RNA'

plotgns=c("si:dkey-261h17.1", "myb", "lmo2", "cebpa", "gata2a", "gata2b", "meis1b", "runx1", "stmn1a", "dnmt3ba", "angpt1", "fgfr2", "dnm3a", "gata1a", "alas2", "igfbp1a", "hemgn", "klf1", "c1qa", "c1qb", "mfap4", "mafba", "mpx", "lyz", "vamp8", "nkl.3", "il2rb", "pax5", "cd37", "cd79b", "sox8b", "tekt1", "spink2.2", "viml")

DotPlot(aggr, features = plotgns) + RotatedAxis()


```


```{r}

#check gata2a
DefaultAssay(aggr)<-'peaks'
Idents(aggr)<-'library_ID'
levels(aggr)<-c("WT" ,"MUT" )
CoveragePlot(
  object = aggr,
  region = "gata2a",
  annotation = TRUE,
  peaks = FALSE)


Idents(aggr)<-'library_ID'
CoveragePlot(
  object = aggr,
  region = "11-3849945-3854463",
  annotation = TRUE,
  peaks = FALSE)
  



```

```{r}
DefaultAssay(aggr) <- 'peaks'
Idents(aggr)<-'cluster_named'
#DN2=CHSPCs, HSPCs=ery-P
levels(aggr) <- c("DN2", "dHSPCs", "HSPCs", "ery", "macs", "neutr", "NKTcells", "Bcells", "multi", "kidney", "endo", "DN1" )
#extend up/downstream arguments added to optimise example plots
CoveragePlot(aggr, region = "mafba")
CoveragePlot(aggr, region = "mpx")
CoveragePlot(aggr, region = "pax5")
CoveragePlot(aggr, region = "hbba1")
CoveragePlot(aggr, region = "il2rb")


```






```{r}
DefaultAssay(aggr) <- 'peaks'

cols <- c("DN2" = "#00BFC4", "dHSPCs" = "orange", "HSPCs" = "#8494FF")

#spi1b
CoveragePlot(
  object = aggr,
  region = "7-32658315-32660197",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)

#gata2a
CoveragePlot(
  object = aggr,
  region = "11-3864838-3866177",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)

#meis1b
CoveragePlot(
  object = aggr,
  region = "13-5568491-5570566",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)

#cd34l
CoveragePlot(
  object = aggr,
  region = "23-32499561-32501198",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)

#cebpa
CoveragePlot(
  object = aggr,
  region = "7-38087323-38088681",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)

#hemgn
CoveragePlot(
  object = aggr,
  region = "1-26667655-26667991",
  idents = c("DN2", "dHSPCs", "HSPCs"),
  annotation = TRUE,
  peaks = FALSE,
  extend.downstream = 500
  
) + scale_fill_manual(values = cols)
```


```{r}
all.genes <- rownames(aggr)
aggr.scale  <- ScaleData(aggr, features = all.genes)


aggr.markers.peaks <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top3peaks <- aggr.markers.peaks %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

open_aggr <- rownames(aggr.markers.peaks[aggr.markers.peaks$avg_logFC > 0.5, ])
closest_genes_open_aggr <- ClosestFeature(aggr.scale, regions = open_aggr)


DoHeatmap(aggr.scale, features = top3peaks$gene) + NoLegend() + scale_fill_gradientn(colors = c("black", "yellow"))

DefaultAssay(aggr.scale) <- 'peaks'
DotPlot(aggr.scale, features = plotgns) + RotatedAxis()

```







```{r}
DefaultAssay(aggr) <- 'peaks'

aggr.markers.peaks <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5peaks <- aggr.markers.peaks %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

open_markers <- rownames(aggr.markers.peaks[aggr.markers.peaks$avg_logFC > 1,])
closest_genes_openmarkers <- ClosestFeature(aggr, regions = open_markers)
```

```{r}

#DN2=cHSPCs
Idents(aggr)<-"cluster_named"

DN2 <- aggr[,grepl("DN2", aggr@active.ident, ignore.case=TRUE)]


Idents(DN2) <- "library_ID"
levels(DN2) <- c("WT","MUT")

cov_plot <- CoveragePlot(
  object = DN2,
  idents = c("WT","MUT"),
  region = "cebpa",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 1500
      )

cov_plot

#npm1a
CoveragePlot(
  object = DN2,
  idents = c("WT","MUT"),
  region = "10-22033673-22036511",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 1500
      )
#lcp1
CoveragePlot(
  object = DN2,
  idents = c("WT","MUT"),
  region = "9-56262771-56274883",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500
      )

#hbba1
CoveragePlot(
  object = DN2,
  idents = c("WT","MUT"),
  region = "hbba1",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500
      )

#hemgn
CoveragePlot(
  object = DN2,
  idents = c("WT","MUT"),
  region = "1-26667655-26667991",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500,
  extend.upstream = 500
      )

markers.DN2 <- FindMarkers(DN2, ident.1 = "WT", ident.2 = "MUT", test.use = 'LR', latent.vars = 'nCount_peaks')

DN2_open_mutant <- rownames(markers.DN2[markers.DN2$avg_logFC < -0.2, ])
closest_genes_DN2_open_mutant <- ClosestFeature(DN2, regions = DN2_open_mutant)


DN2_closed_mutant <- rownames(markers.DN2[markers.DN2$avg_logFC > 0.2, ])
closest_genes_DN2_closed_mutant <- ClosestFeature(DN2, regions = DN2_closed_mutant)


write.table(closest_genes_DN2_closed_mutant, file = "DN2closedmut02.txt")
write.table(closest_genes_DN2_open_mutant, file = "DN2openmut02.txt")



```


```{r}

dHSPCs <- aggr[,grepl("dHSPCs", aggr@active.ident, ignore.case=TRUE)]


Idents(dHSPCs) <- "library_ID"
levels(dHSPCs) <- c("WT","MUT")

cov_plot <- CoveragePlot(
  object = dHSPCs,
  idents = c("WT","MUT"),
  region = "cebpa",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 1500
      )

cov_plot

#npm1a
CoveragePlot(
  object = dHSPCs,
  idents = c("WT","MUT"),
  region = "10-22033673-22036511",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 1500
      )
#lcp1
CoveragePlot(
  object = dHSPCs,
  idents = c("WT","MUT"),
  region = "9-56262771-56274883",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500
      )

#hbba1
CoveragePlot(
  object = dHSPCs,
  idents = c("WT","MUT"),
  region = "hbba1",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500
      )

#hemgn
CoveragePlot(
  object = dHSPCs,
  idents = c("WT","MUT"),
  region = "1-26667655-26667991",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downstream = 500,
  extend.upstream = 500
      )


```

```{r}

ery <- aggr[,grepl("ery", aggr@active.ident, ignore.case=TRUE)]

CoveragePlot(
  object = ery,
  idents = c("WT","MUT"),
  region = "8-21354564-21357062",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = -5000
      )

#nr4a1
CoveragePlot(
  object = ery,
  idents = c("WT","MUT"),
  region = "23-32160089-32165110",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = -50)

#slc25a37
CoveragePlot(
  object = ery,
  idents = c("WT","MUT"),
  region = "8-50189115-50191878",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = -50)



#generate peaks for ery cells to anlayse in Homer

markers.ERY <- FindMarkers(ery, ident.1 = "WT", ident.2 = "MUT", test.use = 'LR', latent.vars = 'nCount_peaks')

ERY_open_mutant <- rownames(markers.ERY[markers.ERY$avg_logFC < -0.2, ])
closest_genes_ERY_open_mutant <- ClosestFeature(ery, regions = ERY_open_mutant)


ERY_closed_mutant <- rownames(markers.ERY[markers.ERY$avg_logFC > 0.2, ])
closest_genes_ERY_closed_mutant <- ClosestFeature(ery, regions = ERY_closed_mutant)

closedmutERY <- data.frame(row.names= closest_genes_ERY_closed_mutant$query_region, group = rep("closed.mut", nrow(closest_genes_ERY_closed_mutant)))

openmutERY <- data.frame(row.names= closest_genes_ERY_open_mutant$query_region, group = rep("open.mut", nrow(closest_genes_ERY_open_mutant)))

write.table(closest_genes_ERY_closed_mutant, file = "ERYclosedmut.txt")
write.table(closest_genes_ERY_open_mutant, file = "ERYopenmut.txt")



```




```{r}


      

levels(aggr)
aggr <- RenameIdents(
  object = aggr,
  'DN2' = 'DN2',
  'dHSPCs' = 'dHSPCs',
  'HSPCs' = 'eryProg',
  'ery' = 'ery',
  'macs' = 'macs',
  'neutr' = 'neutr',
  'NKTcells' = 'NKTcells',
  'Bcells' = 'Bcells',
  'multi' = 'multi',
  'kidney' = 'kidney',
  'endo' = 'endo',
  'DN1' = 'DN1'
  )

eryProg <- aggr[,grepl("eryProg", aggr@active.ident, ignore.case=TRUE)]

```

####to be sorted####


```{r}
markers.DN2 <- FindMarkers(DN2, ident.1 = "WT", ident.2 = "MUT", test.use = 'LR', latent.vars = 'nCount_peaks')

DN2_open_mutant <- rownames(markers.DN2[markers.DN2$avg_logFC < -0.2, ])
closest_genes_DN2_open_mutant <- ClosestFeature(DN2, regions = DN2_open_mutant)


DN2_closed_mutant <- rownames(markers.DN2[markers.DN2$avg_logFC > 0.2, ])
closest_genes_DN2_closed_mutant <- ClosestFeature(DN2, regions = DN2_closed_mutant)


write.table(closest_genes_DN2_closed_mutant, file = "DN2closedmut02.txt")
write.table(closest_genes_DN2_open_mutant, file = "DN2openmut02.txt")




```

```{r}
levels(aggr)
aggr <- RenameIdents(
  object = aggr,
  'DN2' = 'DN2',
  'dHSPCs' = 'dHSPCs',
  'HSPCs' = 'eryProg',
  'ery' = 'ery',
  'macs' = 'macs',
  'neutr' = 'neutr',
  'NKTcells' = 'NKTcells',
  'Bcells' = 'Bcells',
  'multi' = 'multi',
  'kidney' = 'kidney',
  'endo' = 'endo',
  'DN1' = 'DN1'
  )

eryProg <- aggr[,grepl("eryProg", aggr@active.ident, ignore.case=TRUE)]


markers.eryProg <- FindMarkers(eryProg, ident.1 = "WT", ident.2 = "MUT", test.use = 'LR', latent.vars = 'nCount_peaks')

eryProg_open_mutant <- rownames(markers.eryProg[markers.eryProg$avg_logFC < -0.2, ])
closest_genes_eryProg_open_mutant <- ClosestFeature(eryProg, regions = eryProg_open_mutant)


eryProg_closed_mutant <- rownames(markers.eryProg[markers.eryProg$avg_logFC > 0.2, ])
closest_genes_eryProg_closed_mutant <- ClosestFeature(eryProg, regions = eryProg_closed_mutant)


write.table(closest_genes_eryProg_closed_mutant, file = "eryProgclosedmut02.txt")
write.table(closest_genes_eryProg_open_mutant, file = "eryProgopenmut02NEW.txt")



#to subset peaks and examine queary region directly
write.table(closest_genes_eryProg_open_mutant, file = "eryProgopenmut02NEW.txt")
library(dplyr)
eryProgopenmut02NEW$row_num <- seq.int(nrow(eryProgopenmut02NEW))
subset <- eryProgopenmut02NEW %>%
      filter(row_num %in% eryProg_OPEN_GATA1$PositionID)
#eryProg_OPEN_GATA1 is output from homer scanning open peaks for gata motif

#check peak region
cov_plot <- CoveragePlot(
  object = eryProg,
  idents = c("WT","MUT"),
  region = "gata1a",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = 10000
      )

cov_plot

CoverageBrowser(eryProg, region = "gata1a")




```
```{r}
#further homer analysis on ery prog populations
library(dplyr)
#1.ery prog
write.table(closest_genes_eryProg_open_mutant, file = "eryProgopenmut02NEW3.txt")

eryProgopenmut02NEW$row_num <- seq.int(nrow(eryProgopenmut02NEW))
eryPROGsubset_ongata1scl <- eryProgopenmut02NEW %>%
      filter(row_num %in% eryProg_OPEN_gata1scl$PositionID)
#eryProg_OPEN_GATA1 is output from homer scanning open peaks for gata motif

#check peak region

cov_plot <- CoveragePlot(
  object = eryProg,
  idents = c("WT","MUT"),
  region = "alas2",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = -5000
      )

cov_plot

CoverageBrowser(eryProg, region = "alas2")

CoverageBrowser(eryProg, region = "gata1a")

```

```{r}
#further homer analysis on ery prog populations
#1.ery 
write.table(closest_genes_ERY_open_mutant, file = "ERYopennmut.txt")

ERYopennmut$row_num <- seq.int(nrow(ERYopennmut))
ERYsubset_opengata1scl <- ERYopennmut %>%
      filter(row_num %in% ERY_OPEN_gata1scl$PositionID)
#ERY_OPEN_gata1scl is output from homer scanning open peaks for gata motif

ERYsubset_opengata1 <- ERYopennmut %>%
      filter(row_num %in% ERY_OPEN_GATA1$PositionID)
#ERY_OPEN_GATA1 is output from homer scanning open peaks for gata motif



#check peak region

cov_plot <- CoveragePlot(
  object = ery,
  idents = c("WT","MUT"),
  region = "alas2",
  annotation = TRUE,
  peaks = TRUE,
  tile = F,
  extend.downtream = -5000
      )

cov_plot

CoverageBrowser(ery, region = "alas2")

```




```{r}
dHSPCs <- aggr[,grepl("dHSPCs", aggr@active.ident, ignore.case=TRUE)]

current.sample.ids <- c("MUT1","MUT2","MUT3", "WT", "WT1", "WT2", "WT3")
new.sample.ids <- c("MUT","MUT","WT", "WT", "WT", "WT", "MUT")
dHSPCs@meta.data[["library_ID"]] <- plyr::mapvalues(x = dHSPCs@meta.data[["library_ID"]], from = current.sample.ids, to = new.sample.ids)

Idents(dHSPCs) <- "library_ID"
levels(dHSPCs) <- c("WT","MUT")

cov_plot <- CoveragePlot(
  object = ery,
  idents = c("WT","MUT"),
  region = "alas2",
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE,
  extend.upstream = 500
      )

cov_plot

Idents(dHSPCs)="sample"
Idents(dHSPCs)
table(Idents(kidney))

```



```{r}
library(JASPAR2020)
library(TFBSTools)
```

```{r}
Bcells <- aggr[,grepl("Bcells", aggr@active.ident, ignore.case=TRUE)]

current.sample.ids <- c("MUT1","MUT2","MUT3", "WT", "WT1", "WT2", "WT3")
new.sample.ids <- c("MUT","MUT","WT", "WT", "WT", "WT", "MUT")
Bcells@meta.data[["library_ID"]] <- plyr::mapvalues(x = Bcells@meta.data[["library_ID"]], from = current.sample.ids, to = new.sample.ids)

Idents(macs) <- "library_ID"
levels(macs) <- c("WT","MUT")


Idents(Bcells)="library_ID"
Idents(macs)
table(Idents(Bcells))


