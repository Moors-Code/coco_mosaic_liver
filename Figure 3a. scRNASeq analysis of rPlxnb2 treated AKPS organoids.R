######## LOAD PACKAGES AND FUNCTIONS ##############
library(dplyr)
library(Seurat)
library(sctransform)
library(fgsea)
library(msigdbr)
library(rstatix)
library(stringr)
library("Nebulosa")
library(MetBrewer)
library(MoMAColors)

#######LOAD 10X FILES AND PREPROCESS SEURAT OBJECTS####
AKPS_rPlxnb2 <- Read10X(data.dir = "/rPlxnb2/filtered_feature_bc_matrix/")
AKPS_ctrl <- Read10X(data.dir = "/ctrl/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
AKPS_ctrl <- CreateSeuratObject(counts = AKPS_ctrl, project = "AKPS_ctrl", min.cells = 3, min.features = 200)
AKPS_rPlxnb2 <- CreateSeuratObject(counts = AKPS_rPlxnb2, project = "AKPS_rPlxnb2", min.cells = 3, min.features = 200)

#quality control
AKPS_ctrl[["percent.mt"]] <- PercentageFeatureSet(AKPS_ctrl, pattern = "^mt-")
AKPS_rPlxnb2[["percent.mt"]] <- PercentageFeatureSet(AKPS_rPlxnb2, pattern = "^mt-")

VlnPlot(AKPS_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AKPS_rPlxnb2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#merge both datasets
AKPS <- merge(AKPS_ctrl, AKPS_rPlxnb2, add.cell.ids = c("ctrl", "rPlxnb2"))
plot1 <- FeatureScatter(AKPS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AKPS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(AKPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#subset 
AKPS <- subset(AKPS, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & percent.mt < 25)

#normaliye and find variable features
AKPS <- NormalizeData(AKPS, normalization.method = "LogNormalize", scale.factor = 10000)
AKPS <- FindVariableFeatures(AKPS, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AKPS), 10)
all.genes <- rownames(AKPS)

#scale data, run PCA and cluster
AKPS <- ScaleData(AKPS, features = all.genes, vars.to.regress = "percent.mt")
AKPS <- RunPCA(AKPS, features = VariableFeatures(object = AKPS))
DimPlot(AKPS, reduction = "pca")
ElbowPlot(AKPS)
AKPS <- FindNeighbors(AKPS, dims = 1:10)
AKPS <- FindClusters(AKPS, resolution = 0.2)
AKPS <- RunUMAP(AKPS, dims = 1:10)
Idents(AKPS) <- "orig.ident"
DimPlot(AKPS)

#plot clusters
Idents(AKPS) <- "seurat_clusters"
DimPlot(AKPS, split.by = "orig.ident")

#find cluster markers
AKPS.markers <- FindAllMarkers(AKPS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(AKPS.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

#remove cluster 3 and 4 because of mito and contamination, respectively
AKPS_sub <- subset(AKPS, idents = c(0, 1, 2))

#reprocess subsetted data
all.genes <- rownames(AKPS_sub)
AKPS_sub <- ScaleData(AKPS_sub, features = all.genes, vars.to.regress = "percent.mt")
AKPS_sub <- RunPCA(AKPS_sub, features = VariableFeatures(object = AKPS_sub))
DimPlot(AKPS_sub, reduction = "pca")
ElbowPlot(AKPS_sub)
AKPS_sub <- FindNeighbors(AKPS_sub, dims = 1:10)
AKPS_sub <- FindClusters(AKPS_sub, resolution = 0.15)
AKPS_sub <- RunUMAP(AKPS_sub, dims = 1:10)
Idents(AKPS_sub) <- "seurat_clusters"
current.cluster.ids <- c(0, 1, 2)
DimPlot(AKPS_sub, cols=met.brewer("Isfahan2", 3)) + ggsave("Clusterdimplot.pdf", width =6,  height=5)
Idents(AKPS_sub) <- "orig.ident"

#plot
DimPlot(AKPS_sub, reduction = "umap", group.by = "orig.ident", cols=met.brewer("Egypt", 2)) + ggsave("Sampledimplot.pdf", width =6,  height=5)

####COMPOSITIONAL ANALYSIS####
# Calculate frequencies per cluster
numberofcells <- table(AKPS_sub$orig.ident, AKPS_sub$seurat_clusters)
totalcellssample <- colSums(numberofcells)
totalcellspercluster <- c(totalcellssample, sum(totalcellssample))

# Combine the tables
a <- cbind(numberofcells, totalcellssample)
b <- rbind(a, totalcellspercluster)

# Calculate percentages
percentage <- (b[1:2, ] / b[3, ]) * 100

# Create a dataframe with percentages
cluster_percentage <- as.data.frame(t(percentage))
rownames(cluster_percentage) <- c("cluster 0", "cluster 1", "cluster 2")
colnames(cluster_percentage) <- paste0("Sample ", 1:2)

# Plot
par(mar=c(6,8,2,14))
pdf(file="Clusterbreakdown.pdf")
barplot(as.matrix(cluster_percentage), horiz=TRUE,
        legend = TRUE, border=NA,
        args.legend=list(bty = "n", x=130, cex=.8),
        main = "Cluster breakdown per sample", 
        las = 1, 
        col= met.brewer("Isfahan2", 3))
dev.off()


######MARKERS and GO ANALYSIS######
########## GENE SET ENRICHMENT ANALYSIS #############
#import datasets
msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(BP)

#This function performs preranked gene set enrichment analysis (GSEA) using precomputed gene ranks for biological processes (BP).
#The function calculates rankings based on the log fold change (logFC), then uses these rankings to perform GSEA. It returns a dataframe containing 
#the results of preranked GSEA for biological processes, with modified pathway names for readability.
preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- x$genes%>% na.omit()
  head(ranks, 10)
  
  BP_x <- fgsea(pathways = BP, 
                stats = ranks%>% na.omit(),
                minSize=10,
                maxSize=500,
                nperm=1000000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

#calculate cluster markers and GSEA
Idents(AKPS_sub) <- "seurat_clusters"

#cluster 0
DEGs_C0 <- FindMarkers(AKPS_sub, ident.1 = 0, min.pct = 0.25, features = genes.use)
View(DEGs_C0)
BP_C0 <- preranked_BP(DEGs_C0)
View(BP_C0)

ggplot(BP_C0 %>% filter(padj<0.05) %>% head(n= 100), aes(reorder(tolower(pathway), NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(y="Normalized Enrichment Score", x= " ", title="Biological processes in Cluster 0") + 
  theme_classic()+
  theme(axis.text.y=element_text(size=10))# + ggsave("PLxnb2OEvsNT_Hallmarks.pdf", width =7,  height=5)

#cluster 1
DEGs_C1 <- FindMarkers(AKPS_sub, ident.1 = 1, min.pct = 0.25, features = genes.use)
View(DEGs_C1)
BP_C1 <- preranked_BP(DEGs_C1)
View(BP_C1)

ggplot(BP_C1 %>% filter(padj<0.05) %>% head(n= 100), aes(reorder(tolower(pathway), NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(y="Normalized Enrichment Score", x= " ", title="Biological processes in Cluster 1") + 
  theme_classic()+
  theme(axis.text.y=element_text(size=10))

#cluster 2
DEGs_C2 <- FindMarkers(AKPS_sub, ident.1 = 2, min.pct = 0.25, features = genes.use)
View(DEGs_C2)
BP_C2 <- preranked_BP(DEGs_C2)
View(BP_C2%>% filter(padj<0.05))


#######KLF4 AND HRC SIGNATURE##########
#get Klf4 target genes from here https://maayanlab.cloud/Harmonizome/gene_set/KLF4/CHEA+Transcription+Factor+Targets
Klf4_Chea <- read.csv("/media/Coco/MOSAIC LIVER/Experiments/rPlxnb2/Klf4_Chea.csv")
Klf4_targets <- list(c(str_to_title(Klf4_Chea$gene)))
Klf4_targets
AKPS_sub <-AddModuleScore(AKPS_sub, features= Klf4_targets,name = "Klf4_targets")
names(x = AKPS_sub[[]])
FeaturePlot(AKPS_sub, features="Klf4_targets1", split.by = "orig.ident")

DotPlot(AKPS_sub, features= c("Klf4", "Klf4_targets1"), dot.scale = 10 ) +ggsave("KLf4_target_genes.pdf", width = 5, height = 2)
Idents(AKPS_sub) <- "orig.ident"
AKPS_sub_Pl <- subset(AKPS_sub, idents = "AKPS_rPlxnb2")
AKPS_sub_ctrl <- subset(AKPS_sub, idents = "AKPS_ctrl")
p<-wilcox.test(AKPS_sub_Pl$Klf4_targets1, AKPS_sub_ctrl$Klf4_targets1, alternative = "two.sided") #p-value < 2.2e-16
adp <- p.adjust(p$p.value, method = "fdr", n = length(p))
adp


#get coreHRC signature from here https://www.nature.com/articles/s41586-022-05402-9#Sec43
coreHRCsig <- read.csv("/media/Coco/MOSAIC LIVER/Experiments/rPlxnb2/coreHR.csv", sep="")
coreHRC <- list(c(str_to_title(coreHRCsig$gene)))
coreHRC
AKPS_sub <-AddModuleScore(AKPS_sub, features= coreHRC,name = "coreHRC")
names(x = AKPS_sub[[]])
DotPlot(AKPS_sub, features= c("coreHRC1"), dot.scale = 10 ) +ggsave("coreHRC1.pdf", width = 5, height = 2)
VlnPlot(AKPS_sub, features="coreHRC1", group.by = "orig.ident", cols=met.brewer("Egypt", 3), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "core HRC signature", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  ggsave("coreHRC1.pdf", width = 6, height = 6) #p adj = 5.227036e-10
p<-wilcox.test(AKPS_sub_Pl$coreHRC1, AKPS_sub_ctrl$coreHRC1, alternative = "two.sided") #p-value < 2.2e-16
adp <- p.adjust(p$p.value, method = "fdr", n = length(p))
adp


########### EMT SCORE IN UNTREATED AKPS ORGANOIDS ##########
#check EMT and MET scores in untreated AKPS organoids
emt.score <- list(c(BP$GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION))
met.score <- list(c("Basp1", "Bmp4", "Cited1", "Ctnnb1", "Fzd7", "Gata3", "Gdnf", "Grem1", "Lif", 
                  "Pax2", "Pax8", "Pelo", "Sall1", "Six2", "Smo", "Stat1", "Tcf15", "Wnt4", "Wnt9b", 
                  "Wt1", "Klf4", "Ovol1", "Elf3", "Grhl2"))

#subset untreated AKPS
Idents(AKPS_sub) <- "orig.ident"
AKPS_sub_ctrl <- subset(AKPS_sub, idents = "AKPS_ctrl")

#process dataset
AKPS_sub_ctrl <- NormalizeData(AKPS_sub_ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
AKPS_sub_ctrl <- FindVariableFeatures(AKPS_sub_ctrl, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(AKPS_sub_ctrl)
AKPS_sub_ctrl <- ScaleData(AKPS_sub_ctrl, features = all.genes)
AKPS_sub_ctrl <- RunPCA(AKPS_sub_ctrl, features = VariableFeatures(object = AKPS_sub_ctrl))
DimPlot(AKPS_sub_ctrl, reduction = "pca")
AKPS_sub_ctrl <- FindNeighbors(AKPS_sub_ctrl, dims = 1:10)
AKPS_sub_ctrl <- FindClusters(AKPS_sub_ctrl, resolution = 0.1)
AKPS_sub_ctrl <- RunUMAP(AKPS_sub_ctrl, dims = 1:10)
Idents(AKPS_sub_ctrl) <- "seurat_clusters"
DimPlot(AKPS_sub_ctrl)

AKPS_sub_ctrl <-AddModuleScore(AKPS_sub_ctrl, features= emt.score,name = "EMT")
plot_density(AKPS_sub_ctrl, features = "EMT1")  + ggsave("AKPS_ctrl_EMT.pdf", width =5,  height=4)
AKPS_sub_ctrl <-AddModuleScore(AKPS_sub_ctrl, features= met.score,name = "MET")
plot_density(AKPS_sub_ctrl, features = "MET1") + ggsave("AKPS_ctrl_MET.pdf", width =5,  height=4) 

