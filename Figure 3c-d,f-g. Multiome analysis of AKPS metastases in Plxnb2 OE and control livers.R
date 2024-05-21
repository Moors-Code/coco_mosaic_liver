#This code performs the analysis of 4 multiomic single-nucleus RNA and ATAC sequencing datasets generated from livers of mice subjected to intraslenic AKPS injection
#following injection with AAV8-sgPlxnb2OE (2 mice) or with AAV8-sgNT (2 mice). Here's a summary of the key steps and functionalities:

# - Creation of a unified set of peaks: the code begins by reading in ATAC peak sets from different samples and merging them into a unified set of peaks using the GenomicRanges package.
# - Quantification of peaks in each dataset: the code quantifies the counts of these peaks in each dataset, filtering out bad quality peaks based on length criteria.
# - Integration of samples: The datasets are integrated using the FindIntegrationAnchors and IntegrateData functions in Seurat, allowing for joint analysis across different experimental conditions.
# - Annotation of cell types using the RNA space: the RNA assay of the integrated dataset is processed and cell types are annotated using the FindMarkers function.
# - Chromatin accessibility analysis in tumor cells: tumor cells are subsetted and the DNA accessibility data is processed.
# - Motif enrichment analysis: Motif enrichment is calculated to identify enriched transcription factor binding motifs in open chromatin regions.
# - Gene expression analysis after decontamination: The decontX package is used to remove ambient RNA, and differential gene expression analysis is performed on the decontaminated dataset.


######CREATE UNIFIED SET OF PEAKS BEFORE INTEGRATING SAMPLES######
#To create a unified set of peaks we can use functions from the GenomicRanges package. The reduce function from GenomicRanges will merge all intersecting peaks.

#install and load packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("GenomicRanges")

library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(magrittr)
library(dplyr)

# read in peak sets
peaks.OE <- read.table(
  file = "~/nuclei/OE_outs/atac_peaks.bed",
  col.names = c("chr", "start", "end"))

peaks.NT <- read.table(
  file = "~/nuclei/NT_outs/atac_peaks.bed",
  col.names = c("chr", "start", "end"))

peaks.OE.2 <- read.table(
  file = "~/nuclei/OE2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end"))

peaks.NT.2 <- read.table(
  file = "~/nuclei/NT2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end"))

# convert to genomic ranges
peaks.OE <- makeGRangesFromDataFrame(peaks.OE)
peaks.NT <- makeGRangesFromDataFrame(peaks.NT)
peaks.OE.2 <- makeGRangesFromDataFrame(peaks.OE.2)
peaks.NT.2 <- makeGRangesFromDataFrame(peaks.NT.2)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peaks.OE, peaks.NT, peaks.OE.2, peaks.NT.2))

# Filter out bad quality peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.OE <- read.table(
  file = "~/nuclei/OE_outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.NT <- read.table(
  file = "~/nuclei/NT_outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.OE.2 <- read.table(
  file = "~/nuclei/OE2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] 

md.NT.2 <- read.table(
  file = "~/nuclei/NT2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


# create fragment objects
frags.OE <- CreateFragmentObject(
  path = "~/nuclei/OE_outs/atac_fragments.tsv.gz",
  cells = rownames(md.OE),
  validate.fragments = F
)

frags.NT <- CreateFragmentObject(
  path = "~/nuclei/NT_outs/atac_fragments.tsv.gz",
  cells = rownames(md.NT),
  validate.fragments = F
)

frags.OE.2 <- CreateFragmentObject(
  path = "~/nuclei/OE2/outs/atac_fragments.tsv.gz",
  cells = rownames(md.OE.2),
  validate.fragments = F
)

frags.NT.2 <- CreateFragmentObject(
  path = "~/nuclei/NT2/outs/atac_fragments.tsv.gz",
  cells = rownames(md.NT.2),
  validate.fragments = F
)

#Quantify peaks in each dataset
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

OE.counts <- FeatureMatrix(
  fragments = frags.OE,
  features = combined.peaks,
  cells = rownames(md.OE)
)

NT.counts <- FeatureMatrix(
  fragments = frags.NT,
  features = combined.peaks,
  cells = rownames(md.NT)
)

OE.counts.2 <- FeatureMatrix(
  fragments = frags.OE.2,
  features = combined.peaks,
  cells = rownames(md.OE.2)
)

NT.counts.2 <- FeatureMatrix(
  fragments = frags.NT.2,
  features = combined.peaks,
  cells = rownames(md.NT.2)
)

#Create the objects
# get gene annotations 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA data
counts_OE <- Read10X_h5("~/nuclei/OE_outs/filtered_feature_bc_matrix.h5")
OE <- CreateSeuratObject(
  counts = counts_OE$`Gene Expression`,
  assay = "RNA", project = "OE" 
)

# create ATAC assay and add it to the object
OE[["ATAC"]] <- CreateChromatinAssay(
  counts = OE.counts[, colnames(counts_OE$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.OE,
  annotation = annotation
)
OE

#call peaks using MACS2
DefaultAssay(OE) <- "ATAC"
peaks_OE <- CallPeaks(OE, macs2.path="/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/nuclei"),
                      fragment.tempdir = tempdir("~/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_OE <- keepStandardChromosomes(peaks_OE, pruning.mode = "coarse")
peaks_OE <- subsetByOverlaps(x = peaks_OE, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(OE),
  features = peaks_OE,
  cells = colnames(OE)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
OE[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = frags.OE,
  annotation = annotation
)

#now for NT
counts_NT <- Read10X_h5("/nuclei/NT_outs/filtered_feature_bc_matrix.h5")

NT <- CreateSeuratObject(
  counts = counts_NT$`Gene Expression`,
  assay = "RNA", project = "NT" 
)

# create ATAC assay and add it to the object
NT[["ATAC"]] <- CreateChromatinAssay(
  counts = NT.counts[, colnames(counts_NT$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.NT,
  annotation = annotation
)
NT


#call peaks using MACS2
DefaultAssay(NT) <- "ATAC"
peaks_NT <- CallPeaks(NT, macs2.path="/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/nuclei"),
                      fragment.tempdir = tempdir("~/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_NT <- keepStandardChromosomes(peaks_NT, pruning.mode = "coarse")
peaks_NT <- subsetByOverlaps(x = peaks_NT, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts_NT <- FeatureMatrix(
  fragments = Fragments(NT),
  features = peaks_NT,
  cells = colnames(NT)
)

# create a NT assay using the MACS2 peak set and add it to the Seurat object
NT[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_NT,
  fragments = frags.NT,
  annotation = annotation
)


#now for OE.2
counts_OE.2 <- Read10X_h5("/nuclei/OE2/outs/filtered_feature_bc_matrix.h5")
OE.2 <- CreateSeuratObject(
  counts = counts_OE.2$`Gene Expression`,
  assay = "RNA", project = "OE.2" 
)

# create ATAC assay and add it to the object
OE.2[["ATAC"]] <- CreateChromatinAssay(
  counts = OE.counts.2[, colnames(counts_OE.2$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.OE.2,
  annotation = annotation
)
OE.2

#call peaks using MACS2
DefaultAssay(OE.2) <- "ATAC"
peaks_OE.2 <- CallPeaks(OE.2, macs2.path="/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/nuclei"),
                      fragment.tempdir = tempdir("~/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_OE.2 <- keepStandardChromosomes(peaks_OE.2, pruning.mode = "coarse")
peaks_OE.2 <- subsetByOverlaps(x = peaks_OE.2, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts.OE2 <- FeatureMatrix(
  fragments = Fragments(OE.2),
  features = peaks_OE.2,
  cells = colnames(OE.2)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
OE.2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts.OE2,
  fragments = frags.OE.2,
  annotation = annotation
)

#now for NT.2
counts_NT.2 <- Read10X_h5("/media/Coco/MOSAIC LIVER/Experiments/nuclei/NT2/outs/filtered_feature_bc_matrix.h5")

NT.2 <- CreateSeuratObject(
  counts = counts_NT.2$`Gene Expression`,
  assay = "RNA", project = "NT.2" 
)

# create ATAC assay and add it to the object
NT.2[["ATAC"]] <- CreateChromatinAssay(
  counts = NT.counts.2[, colnames(counts_NT.2$`Gene Expression`)],
  sep = c(":", "-"),
  fragments = frags.NT.2,
  annotation = annotation
)
NT.2


#call peaks using MACS2
DefaultAssay(NT.2) <- "ATAC"
peaks_NT.2 <- CallPeaks(NT.2, macs2.path="/miniconda3/envs/macs2-new/bin/macs2", outdir = tempdir("~/nuclei"),
                      fragment.tempdir = tempdir("~/nuclei"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_NT.2 <- keepStandardChromosomes(peaks_NT.2, pruning.mode = "coarse")
peaks_NT.2 <- subsetByOverlaps(x = peaks_NT.2, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts_NT.2 <- FeatureMatrix(
  fragments = Fragments(NT.2),
  features = peaks_NT.2,
  cells = colnames(NT.2)
)

# create a NT assay using the MACS2 peak set and add it to the Seurat object
NT.2[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_NT.2,
  fragments = frags.NT.2,
  annotation = annotation
)


#integrate samples
list <- list(OE,NT, OE.2, NT.2)
list.anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
integrated <- IntegrateData(anchorset = list.anchors, dims = 1:30)
saveRDS(integrated, file = "integrated.rds")

######## ANNOTATE CELL TYPES USING THE RNA SPACE ##########
DefaultAssay(integrated) <- "RNA"
integrated <- SCTransform(integrated)
integrated <- RunPCA(integrated)
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.3)
integrated <- RunUMAP(integrated, dims = 1:20, return.model=TRUE)
DimPlot(integrated, label = T)
DimPlot(integrated, label = T, group.by = "sample")
DimPlot(integrated)

#Annotating cell types
markers_integrated <- FindAllMarkers(object = integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_integrated %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
FeaturePlot(integrated, features = "Glul", order=T)
FeaturePlot(integrated, features = "Lgr5")
FeaturePlot(integrated, features = "Cyp2f2")
FeaturePlot(integrated, features = "Plxnb2", order=T)
FeaturePlot(integrated, features = "Epcam", order=T)
FeaturePlot(integrated, features = "Gpx2", order=T)
FeaturePlot(integrated, features = "Gpx2", order=T)

#rename clusters
current.cluster.ids <- c(0, 1, 2, 3, 4,5,6, 7,8,9,10, 11, 12, 13, 14,15,16)
new.cluster.ids <-  c("hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "endothelial cells",
                      "Kupffer cells", 
                      "hepatocytes",
                      "tumor cells", 
                      "hepatocytes",
                      "fibroblasts",
                      "hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "hepatocytes", 
                      "macrophages",
                      "lymphocytes")
integrated@meta.data$seurat_clusters <- plyr::mapvalues(x = integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
Idents(integrated) <- "seurat_clusters"
DimPlot(integrated, reduction = "umap", pt.size = .5, cols = rev(met.brewer("Hokusai1"))) + ggsave("Clusterdimplot.pdf", width =6,  height=5)

integrated$condition <- ifelse(test = integrated$sample %in% c("sgPlxnb2_1", "sgPlxnb2_2"), yes = "sgPlxnb2 OE", no = "sgNT")

####subset cell types
hepatocytes <- subset(integrated, idents = "hepatocytes")
tumor <- subset(integrated, idents = "tumor cells")

##### ATAC SEQ ANALYSIS IN TUMOR CELLS ############
#DNA accessibility data processing
DefaultAssay(tumor) <- "peaks"
tumor <- FindTopFeatures(tumor, min.cutoff = 5)
tumor <- RunTFIDF(tumor)
tumor <- RunSVD(tumor)
tumor <- RunUMAP(tumor, dims = 2:50, reduction = 'lsi')
DimPlot(tumor, group.by = 'condition', pt.size = 0.5)
DimPlot(tumor, cols=c("#8DC63F","#A7A9AC"), pt.size = 0.5) + ggsave("tumorATACdimplot.pdf", width =5,  height=3)
 
# build a joint neighbor graph using both assays
tumor <- FindMultiModalNeighbors(
  object = tumor,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
tumor <- RunUMAP(
  object = tumor,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(tumor, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend() 


# first compute the GC content for each peak
tumor <- RegionStats(tumor, genome = BSgenome.Mmusculus.UCSC.mm10)

# link peaks to genes
tumor <- LinkPeaks(
  object = tumor,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("Cdh1")
)

#We can visualize these links using the CoveragePlot() function, 
Idents(combined) <- "sample"
CoveragePlot(
  object = tumor,
  region = "Cdh1",
  features = "Cdh1",
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 10000
)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# add motif information
tumor <- AddMotifs(
  object = tumor,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = tumor,
  ident.1 = "sgPlxnb2 OE",
  ident.2 = "sgNT",
  only.pos = T,  #or set F, but why is F taking so much longer?
  test.use = 'LR', 
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# find peaks open 
open.peaks <- AccessiblePeaks(tumor, idents = c("sgPlxnb2 OE"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(tumor, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = tumor,
  features = top.da.peak,
  background=peaks.matched
)
View(enriched.motifs)

MotifPlot(
  object = tumor,
  motifs = head(rownames(enriched.motifs))
) + ggsave("Motifplot.pdf", width =8,  height=5)

############GENE EXPRESSION ANALYSIS##########
#here we used the package decontX to remove ambient RNA

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("celda")
library(celda)
vignette("decontX")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("eds")

install.packages("remotes")
remotes::install_github("compbiomed/singleCellTK")
library(singleCellTK)
library(Seurat)
library(magrittr)
library(dplyr)

#dexontX needs a single cell Experiment object 
integrated <- readRDS("~/NAS/Coco/MOSAIC LIVER/Experiments/nuclei/integrated.rds")
sce <- as.SingleCellExperiment(integrated)
sce.delta <- decontX(sce)

#check decontamination 
plotDecontXContamination(sce.delta)

#decontaminated count matrix 
head(decontXcounts(sce.delta))

#convert back to Seurat object
dec_integrated <- convertSCEToSeurat(sce.delta, copyDecontX = TRUE,copyColData = TRUE)
dec_integrated <- SCTransform(dec_integrated)
dec_integrated <- RunPCA(dec_integrated)
dec_integrated <- FindNeighbors(dec_integrated, dims = 1:20)
dec_integrated <- FindClusters(dec_integrated, resolution = 0.3)
dec_integrated <- RunUMAP(dec_integrated, dims = 1:20, return.model=TRUE)
DimPlot(dec_integrated, label = T)
DimPlot(dec_integrated, label = T, group.by = "sample")
FeaturePlot(dec_integrated, features = "decontX_contamination")

dec_integrated$contaminated <- ifelse(test = dec_integrated$decontX_contamination > 0.75, yes = "YES", no = "NO")
DimPlot(dec_integrated, label = T, group.by = "contaminated")
Idents(dec_integrated) <- "contaminated"
sub_dec_integrated <- subset(dec_integrated, idents = "NO")

DimPlot(sub_dec_integrated, group.by = "sample")
sub_dec_integrated <- SCTransform(sub_dec_integrated)
sub_dec_integrated <- RunPCA(sub_dec_integrated)
sub_dec_integrated <- FindNeighbors(sub_dec_integrated, dims = 1:20)
sub_dec_integrated <- FindClusters(sub_dec_integrated, resolution = 0.3)
sub_dec_integrated <- RunUMAP(sub_dec_integrated, dims = 1:20, return.model=TRUE)
DimPlot(sub_dec_integrated, label = T)
DimPlot(sub_dec_integrated, label = T, group.by = "sample")

#Annotating cell types
markers_sub_dec_integrated <- FindAllMarkers(object = sub_dec_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_sub_dec_integrated %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))
FeaturePlot(sub_dec_integrated, features = "Glul", order=T)
FeaturePlot(sub_dec_integrated, features = "Lgr5")
FeaturePlot(sub_dec_integrated, features = "Cyp2f2")
FeaturePlot(sub_dec_integrated, features = "Plxnb2", order=T)
FeaturePlot(sub_dec_integrated, features = "Epcam", order=T)
FeaturePlot(sub_dec_integrated, features = "Gpx2", order=T)
Idents(sub_dec_integrated) <- "seurat_clusters"
cluster_markers_sub_dec_integrated <- FindAllMarkers(object = sub_dec_integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
View(cluster_markers_sub_dec_integrated %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC))


#rename clusters
current.cluster.ids <- c(0, 1, 2, 3, 4,5,6, 7,8,9,10, 11, 12)
new.cluster.ids <-  c("hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "hepatocytes",
                      "endothelial cells",
                      "Kupffer cells", 
                      "hepatocytes",
                      "tumor cells", 
                      "hepatocytes",
                      "fibroblasts",
                      "hepatocytes",
                      "lymphocytes",
                      "macrophages"
                     )

sub_dec_integrated@meta.data$seurat_clusters <- plyr::mapvalues(x = sub_dec_integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sub_dec_integrated, reduction = "umap", pt.size = .3, cols = rev(met.brewer("Hokusai1"))) + ggsave("Clusterdimplot.pdf", width =6,  height=5)

sub_dec_integrated$condition <- ifelse(test = sub_dec_integrated$sample %in% c("sgPlxnb2_1", "sgPlxnb2_2"), yes = "sgPlxnb2 OE", no = "sgNT")
DimPlot(sub_dec_integrated, split.by = "condition")

#create color palette
library(MetBrewer)
pal <- met.brewer(name="Archambault", n=123, brew_type="continuous")

#plot cell type markers
library(ggplot2)
markers <- c("Cyp3a25", "Apoa2", "Pecam1", "Ptprb", "Clec4f", "Cd5l", "Lgr5", 
             "Epcam", "Col14a1", "Nrxn3", "Cd3e", "Cd4", "Cd8a","Cxcr5", "Spp1", "Gas6")
DotPlot(sub_dec_integrated, features = markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  theme(legend.position="right")  + scale_colour_gradientn(colours = pal)+
  labs(title = "cluster markers", y = "", x="")+
  ggsave("RNA_dotplot.pdf", width = 8, height = 3.5)


####subset cell types
hepatocytes <- subset(sub_dec_integrated, idents = "hepatocytes")
tumor <- subset(sub_dec_integrated, idents = "tumor cells")

#check Plxnb2 upregulation in hepatocytes
hepatocytes <- SCTransform(hepatocytes)
hepatocytes <- RunPCA(hepatocytes)
hepatocytes <- FindNeighbors(hepatocytes, dims = 1:20)
hepatocytes <- FindClusters(hepatocytes, resolution = 0.3)
hepatocytes <- RunUMAP(hepatocytes, dims = 1:20, return.model=TRUE)
saveRDS(hepatocytes, file="hepatocytes.rds")
DimPlot(hepatocytes, label = T, group.by = "sample")
VlnPlot(hepatocytes, features = "Plxnb2", log = T, pt.size = 0, group.by = "sample") +ylim(0.9, 3)
hepatocytes$condition <- ifelse(test = hepatocytes$sample %in% c("sgPlxnb2_1", "sgPlxnb2_2"), yes = "sgPlxnb2 OE", no = "sgNT")
DotPlot(hepatocytes, features = c("Plxnb2"), group.by = "condition")

Idents(hepatocytes) <-"condition"
VlnPlot(hepatocytes, features= "Plxnb2", group.by = "condition", pt.size = 0, adjust=3) +ylim(0, 3) +theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl", geom = "pointrange", color = "black")


#differential gene expression in tumor cells
tumor <- SCTransform(tumor)
tumor <- RunPCA(tumor)
tumor <- FindNeighbors(tumor, dims = 1:20)
tumor <- FindClusters(tumor, resolution = 0.3)
tumor <- RunUMAP(tumor, dims = 1:20, return.model=TRUE)
DimPlot(tumor, label = T, group.by = "condition")
tumor$condition <- ifelse(test = tumor$sample %in% c("sgPlxnb2_1", "sgPlxnb2_2"), yes = "sgPlxnb2 OE", no = "sgNT")
Idents(tumor) <-"condition"
tumor_OE_vs_NT_markers <- FindMarkers(object = tumor, ident.1 =  "sgPlxnb2 OE", ident.2 ="sgNT",
                                      only.pos = F, min.pct = 0.25, logfc.threshold = 0)
View(tumor_OE_vs_NT_markers)
tumor_OE_vs_NT_markers$genes <- rownames(tumor_OE_vs_NT_markers)
DotPlot(tumor, features = c("Cdh1", "Epcam", "Ctnnb1"), group.by = "condition")

#generate volcano plot
DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$genes <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), DEGs$genes, NA)
  
  plot <-  ggplot(data=DEGs, aes(x=avg_log2FC, y=-log10(p_val_adj), fill=factor(Significant), label = DEGs$genes) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=3, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10(p_val_adj)") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}
DEGs_volcano(tumor_OE_vs_NT_markers, 0.05, 0.75, "Plxnb2-OE vs NT", "grey",150, 1.7)+ggsave("Plxnb2OEvsNT_volcano.pdf", width =7,  height=5)


#KLF4 target genes
Klf4_Chea <- read.csv("Klf4_Chea.csv")
Klf4_targets <- list(c(str_to_title(Klf4_Chea$gene)))
Klf4_targets
tumor <-AddModuleScore(tumor, features= Klf4_targets,name = "Klf4_targets")
names(x = tumor[[]])
VlnPlot(tumor, features="Klf4_targets1", group.by = "condition", cols=met.brewer("Egypt", 2), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Klf4 target gene expression", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("KLf4_target_genes.pdf", width = 6, height = 6)
DotPlot(tumor, features= c("Klf4", "Klf4_targets1"), dot.scale = 10 ) +ggsave("KLf4_target_genes.pdf", width = 5, height = 2)

#HRC signature
coreHRCsig <- read.csv("coreHR.csv", sep="")
coreHRC <- list(c(str_to_title(coreHRCsig$gene)))
coreHRC
tumor <-AddModuleScore(tumor, features= coreHRC,name = "coreHRC")
names(x = tumor[[]])
DotPlot(tumor, features= c("coreHRC1"), dot.scale = 10 ) #+ggsave("coreHRC1.pdf", width = 5, height = 2)
VlnPlot(tumor, features="coreHRC1", group.by = "condition", cols=met.brewer("Egypt", 3), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "core HRC signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("coreHRC1.pdf", width = 6, height = 6) #p adj = 5.227036e-10


#iCMS3 signature
iCMS3sig <- read.csv("iCMS3.csv", sep="")
iCMS3 <- list(c(str_to_title(iCMS3sig$gene)))
iCMS3
tumor <-AddModuleScore(tumor, features= iCMS3,name = "iCMS3")
names(x = tumor[[]])
DotPlot(tumor, features= c("iCMS31"), dot.scale = 10 ) #+ggsave("coreHRC1.pdf", width = 5, height = 2)
VlnPlot(tumor, features="iCMS31", group.by = "condition", cols=met.brewer("Egypt", 3), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "iCMS3 signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("iCMS31.pdf", width = 6, height = 6) #p adj = 5.227036e-10


#gene set enrichment analysis
library(fgsea)
library(msigdbr)

#import datasets
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)

preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=avg_log2FC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  
  BP_x <- fgsea(pathways = BP, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

BP_tumor <- preranked_BP(tumor_OE_vs_NT_markers)
View(BP_tumor)

library(RColorBrewer) 
ggplot(BP_tumor %>% filter(abs(NES)>1 & pval<0.035) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=10)) + ggsave("tumor_PLxnb2OEvsNT_BP.pdf", width =20,  height=5)

