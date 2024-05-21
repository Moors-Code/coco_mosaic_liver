#This code uses Seurat to analyze the Molecular Cartography dataset published in Handler et al, 2023.
#See here for detailed upstream processing and analysis of the dataset: https://github.com/Moors-Code/Fragment-sequencing/tree/main/Highly_multiplexed_FISH_analysis

library(Seurat)
library(ggplot2)
library(MetBrewer)

#load the dataset
Resolve_seurat_anno <- readRDS("Resolve_seurat_anno.rds")
DimPlot(Resolve_seurat_anno)
DimPlot(Resolve_seurat_anno, group.by = "samples")
Idents(Resolve_seurat_anno) <- "samples"

#split in metastatic livers and naive livers 
no_mets <- subset(Resolve_seurat_anno, idents = "no_mets")
mets <- subset(Resolve_seurat_anno, idents = "mets")

#subset hepatocytes
Idents(Resolve_seurat_anno) <- "annotation"
hepatocytes <- subset(Resolve_seurat_anno, idents = c("Hepatocytes_CV", "Hepatocytes_PV"))

#plot Plxnb2 expression in central vs portal hepatocytes (Extended Data Fig 5g)
VlnPlot(hepatocytes, features="Plxnb2", cols=met.brewer("Archambault", 2), pt.size = 0)+
  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + ylim(0.01,4)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(title = "", y = "Semaphorin expression", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
  ggsave("Plxnb2_PV_CV-8.337577e-164.pdf", width =5,  height=5)

dePlxnb2 <- FindMarkers(hepatocytes, ident.1= "Hepatocytes_PV", ident.2="Hepatocytes_CV", logfc.threshold=0)
View(dePlxnb2)


####plot Plxnb2 and class IV semaphorin expression in metastatic liver (Extended Data Fig 9a)
VlnPlot(mets, features = c("Sema4a", "Sema4c", "Sema4d", "Sema4g", "Plxnb2"), cols = c("yellow", "red", "green", "violet", "blue"), stack = TRUE, sort = TRUE,  flip = TRUE) + ylim(0.001,7) +
  theme(legend.position = "none") +ggsave("Semaexpression.pdf", width=6, height = 5)
