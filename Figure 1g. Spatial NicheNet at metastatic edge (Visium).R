### This code by Atefeh Lafzi is used to identify tumor-hepatocyte interactions at the metastatic leading edge

#####DATA LOADING AND PROCESSING WITH STUTILITY
#This code segment loads a Seurat object from a saved file, prepares an infotable for STUtility, creates a STUtility object, 
#and processes it to load spatial transcriptomics data. Then, it assigns predicted cell identities from the loaded Seurat object, 
#performs SCTransform normalization, and loads images into the Seurat object. Finally, it sets identity classes and overlays predicted cell identities on images.

library(magrittr)
library(dplyr)
library(STutility)

#load Visium data of human CRC hepatic metastasis, spots bulk data is deconvoluted with SPOTlight (Elodua-Bayes et al, 2021) using two published scRNAseq datasets (Lee et al, 2020 and Massalha et al, 2020)
#the analysis below was performed on two technical replicates of the same metastatic sample (both datasets availbale on Zenodo)
load("Met1rep1_seuratObj_V2.RData")

# Define infotable for STUtility
infoTable <- data.frame(samples=NA, spotfiles=NA, imgs=NA, json=NA, condition=NA)
infoTable <- rbind(infoTable, data.frame(samples="B1_FRB33CB70_Met/outs/filtered_feature_bc_matrix.h5", 
                                         spotfiles="B1_FRB33CB70_Met/outs/spatial/tissue_positions_list.csv", 
                                         imgs="B1_FRB33CB70_Met/outs/spatial/tissue_hires_image.png", 
                                         json="B1_FRB33CB70_Met/outs/spatial/scalefactors_json.json", condition="Met1"))
dim(infoTable)
infoTable <- infoTable[-1,]

STU <- InputFromTable(infotable = infoTable, 
                      min.gene.count = 0, 
                      min.gene.spots = 0,
                      min.spot.count = 0,
                      platform =  "Visium")

Met1STU <- STU
Met1.temp <- Met1
rownames(Met1.temp@meta.data) <- paste(rownames(Met1.temp@meta.data), "_1", sep = "")
Met1STU <- AddMetaData(Met1STU, Met1temp@meta.data[rownames(Met1STU@meta.data), "seurat_predicted.id"], "seurat_predicted.id")
Met1STU <- SCTransform(Met1STU)
Met1STU <- LoadImages(Met1STU, time.resolve = F, verbose = T)

Met1STU <- SetIdent(Met1STU, value = "seurat_predicted.id")
FeatureOverlay(Met1STU, features = "seurat_predicted.id", pt.size = 2, dark.theme = F)


######## IDENTIFY METASTASIS EDGE SPOTS ################
#This code segment identifies neighboring spots with respect to hepatocytes, overlays hepatocyte features on a spatial plot, 
#assigns labels to spots as either tumors or hepatocytes based on their interactions, adds these labels to the metadata, 
#and overlays the tumor or hepatocyte labels on the spatial plot for visualization.

# Identify neighboring spots with respect to hepatocytes
Met1STU <- RegionNeighbours(Met1STU, id = "Hepatocytes", verbose = TRUE)

# Overlay hepatocyte features on the spatial plot
FeatureOverlay(Met1STU, features = "nbs_Hepatocytes", cols = c("red", "lightgray"), pt.size = 2, dark.theme = F)

# Initialize vector to label spots as tumor or hepatocytes based on interactions
Hepa.nb.tumor <- rep(NA, nrow(Met1STU@meta.data))

# Label spots as tumor or hepatocytes based on conditions
Hepa.nb.tumor[which(Met1STU@meta.data$seurat_predicted.id == "Tumor_Carcinoma" & Met1STU@meta.data$nbs_Hepatocytes == "nbs_Hepatocytes")] <- "Tumor_nb"
Hepa.nb.tumor[which(Met1STU@meta.data$seurat_predicted.id == "Hepatocytes" & Met1STU@meta.data$nbs_Hepatocytes == "Hepatocytes")] <- "Hepatocytes"

# Add the labels to metadata
Met1STU <- AddMetaData(Met1STU, Hepa.nb.tumor, "nbs_Hepatocytes_Tumors")

# Overlay the tumor or hepatocyte labels on the spatial plot
FeatureOverlay(Met1STU, features = "nbs_Hepatocytes_Tumors", cols = c("red", "yellow"), pt.size = 2, dark.theme = F)

# Initialize vector to label spots as tumor-central or tumor-edge
Hepa.tumor.inner.edge <- rep(NA, nrow(Met1STU@meta.data))

# Label spots as tumor-central or tumor-edge based on conditions
for (i in 1:nrow(Met1STU@meta.data)){
  if (Met1STU@meta.data[i,"seurat_predicted.id"] == "Hepatocytes" & is.na(Met1STU@meta.data[i,"nbs_Hepatocytes_Tumors"])) {Hepa.tumor.inner.edge[i] <- "Hepatocyte_central"}
  else if (Met1STU@meta.data[i,"seurat_predicted.id"] == "Hepatocytes" & Met1STU@meta.data[i,"nbs_Hepatocytes_Tumors"] == "Hepatocytes") {Hepa.tumor.inner.edge[i] <- "Hepatocyte_TumorNeighbour"}
  if (Met1STU@meta.data[i,"seurat_predicted.id"] == "Tumor_Carcinoma" & is.na(Met1STU@meta.data[i,"nbs_Hepatocytes_Tumors"])) {Hepa.tumor.inner.edge[i] <- "Tumor_central"}
  else if (Met1STU@meta.data[i,"seurat_predicted.id"] == "Tumor_Carcinoma" & Met1STU@meta.data[i,"nbs_Hepatocytes_Tumors"] == "Tumor_nb") {Hepa.tumor.inner.edge[i] <- "Tumor_HepatoNeighbour"}
}
Met1STU <- AddMetaData(Met1STU, Hepa.tumor.inner.edge, "Hepa_tumor_inner_edge")

# Overlay tumor features on the spatial plot
FeatureOverlay(Met1STU, features = "Hepa_tumor_inner_edge", cols = c("green", "yellow", "red", "blue","grey"), pt.size = 2, dark.theme = F)


########DIFFERENTIAL GENE EXPRESSION################
#This section of the code performs differential gene expression analysis between hepatocyte-tumor neighbors and hepatocyte central spots, 
#as well as between tumor-hepatocyte neighbors and tumor central spots. It identifies markers for each group and visualizes them using heatmaps. 
#Finally, the results are saved for further analysis.

# Set identity class for further analysis
Met1STU <- SetIdent(Met1STU, value = "Hepa_tumor_inner_edge")

# Find markers for hepatocyte-tumor neighbor vs. hepatocyte central
hepato.edge.central.markers.Met1.rep <- FindMarkers(Met1STU, ident.1 = "Hepatocyte_TumorNeighbour", ident.2 = "Hepatocyte_central")

# Find markers for tumor-hepatocyte neighbor vs. tumor central
Tumor.edge.central.markers.Met1.rep <- FindMarkers(Met1STU, ident.1 = "Tumor_HepatoNeighbour", ident.2 = "Tumor_central")

# Save the results
save(hepato.edge.central.markers.Met1.rep, file = "DEGs_HepatoNeighbouringTumor_VS_HepatoCentral_Met1rep_V2.RData")
save(Tumor.edge.central.markers.Met1.rep, file = "DEGs_TumorNeighbouringHepato_VS_TumorCentral_Met1rep_V2.RData")

# Add gene names as a column for further analysis
hepato.edge.central.markers.Met1.rep$gene <- rownames(hepato.edge.central.markers.Met1.rep)

# Subset the Seurat object to include only the relevant cell subset
hepato.edge.central.Met1STU.subset <- SubsetSTData(Met1STU, expression = Hepa_tumor_inner_edge %in% c("Hepatocyte_TumorNeighbour", "Hepatocyte_central"))

# Sort markers by absolute average log2 fold change and select top markers
hepato.edge.central.Met1sorted.marks <- hepato.edge.central.markers.Met1.rep %>% top_n(n = 40, wt = abs(avg_log2FC))
hepato.edge.central.Met1sorted.marks <- hepato.edge.central.Met1sorted.marks[order(hepato.edge.central.Met1sorted.marks$avg_log2FC, decreasing = TRUE), ]

# Generate heatmap for the subset with selected markers
DoHeatmap(hepato.edge.central.Met1STU.subset, features = hepato.edge.central.Met1sorted.marks$gene, group.colors = c("red", "lightgray"), disp.min = -2, disp.max = 2) 

# Add gene names as a column for further analysis
Tumor.edge.central.markers.Met1.rep$gene <- rownames(Tumor.edge.central.markers.Met1.rep)

# Subset the Seurat object to include only the relevant cell subset
Tumor.edge.central.Met1STU.subset <- SubsetSTData(Met1STU, expression = Hepa_tumor_inner_edge %in% c("Tumor_HepatoNeighbour", "Tumor_central"))

# Sort markers by absolute average log2 fold change and select top markers
Tumor.edge.central.Met1sorted.marks <- Tumor.edge.central.markers.Met1.rep %>% top_n(n = 40, wt = abs(avg_log2FC))
Tumor.edge.central.Met1sorted.marks <- Tumor.edge.central.Met1sorted.marks[order(Tumor.edge.central.Met1sorted.marks$avg_log2FC, decreasing = TRUE), ]

# Generate heatmap for the subset with selected markers
DoHeatmap(Tumor.edge.central.Met1STU.subset, features = Tumor.edge.central.Met1sorted.marks$gene, group.colors = c("red", "lightgray"), disp.min = -2, disp.max = 2) 

#Volcano plot (Extended Data Fig 2f)

DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$Gene <- ifelse((DEGs$p_val_adj < p_treshold & abs(DEGs$avg_log2FC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=avg_log2FC, y=-log10(p_val_adj), fill=factor(Significant), label = DEGs$Gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=3, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}

DEGs_volcano(Tumor.edge.central.markers, 0.05, 0.5, "tumor edge vs central", "grey",20,1) +ggsave("DEGsmet1_tumor.pdf", width = 5, height = 5)



######### NICHENET ON NEIGHBORING SPOTS ###########
#This code segment performs gene expression analysis, predicts ligand activities, identifies active 
#ligand-target links, and prepares visualization data for ligand-target interactions. Finally, it creates a 
#heatmap plot for the ligand-target network and customizes its appearance.

library(nichenetr)
library(sparseMatrixStats)
Met1nichenet <- Met1STU
Met1nichenet@meta.data %>% head()

# Set cell identities to "Hepa_tumor_inner_edge"
Idents(Met1nichenet) <- "Hepa_tumor_inner_edge"

# Run PCA and UMAP on the Seurat object
Met1nichenet <- RunPCA(Met1nichenet, verbose = FALSE)  %>% RunUMAP(dims = 1:30)
DimPlot(Met1nichenet, reduction = "umap")

# Load ligand-target matrix from file
ligand_target_matrix = readRDS("ligand_target_matrix.rds")

# Display dimensions and a subset of the ligand-target matrix
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
dim(ligand_target_matrix)

# Transpose and convert RNA counts matrix to a dense matrix
expression = t(as.matrix(Met1nichenet@assays$RNA@counts))
sample_info = Met1nichenet@meta.data

# Extract IDs for Hepatocytes and Tumor spots based on cell identity
Hepato_ids = rownames(sample_info[which(sample_info$Hepa_tumor_inner_edge == "Hepatocyte_TumorNeighbour"),])
Tumor_ids = rownames(sample_info[which(sample_info$Hepa_tumor_inner_edge == "Tumor_HepatoNeighbour"),])

# Extract expressed genes for sender (Hepatocytes) and receiver (Tumor) cells
expressed_genes_sender = colnames(expression[Hepato_ids,])[which(colSums2(expression[Hepato_ids,]) > length(Hepato_ids)*0.1)]
expressed_genes_receiver = colnames(expression[Tumor_ids,])[which(colSums2(expression[Tumor_ids,]) > length(Tumor_ids)*0.1)]
length(expressed_genes_sender)
length(expressed_genes_receiver)

# Define gene set of interest and background expressed genes
geneset_oi = rownames(Tumor.edge.central.markers)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
length(geneset_oi)
length(background_expressed_genes)
length(intersect(geneset_oi, background_expressed_genes))

# Extract potential ligands and receptors based on lr_network
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

# Filter lr_network to include only expressed ligand-receptor pairs
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
save(lr_network_expressed, file= "Hepato_Tumor_all_expressed_ligandReceptorNetwork_NoFiltering.RData")


#Ligand activity analysis (Extended Data Fig 2g): calculate the ligand activity of each ligand to assess how well each 
#Hepatocyte-derived LR can predict the tumor edge gene set compared to the background of expressed genes 
# Extract potential ligands and calculate ligand activities
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson)

# Identify best upstream ligands based on Pearson correlation
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

# Extract active ligand-target links and prepare for visualization
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
head(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
head(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

# Save ligand activities, active ligand-target links, and visualization data
save(ligand_activities, file="Hepato_Tumor_active_ligand_activities_weighted_all_V2.RData")
save(active_ligand_target_links_df, file="Hepato_Tumor_active_ligand_target_links_df_V2.RData")
save(vis_ligand_target, file="Hepato_Tumor_vis_ligand_target_V2.RData")

# Create a heatmap plot for ligand-target network (Extended Data Fig 2g)
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Hepatocytes-ligands","Tumor_Tumor genes in Metastatic spots", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
                                                  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))
p_ligand_target_network +theme(text = element_text(size=15, face = "bold"),axis.text.x = element_text(size=15,face="bold"))


###########CALCULATE INTERSECTION WITH SCREENING HITS (Figure 1g) ###########
#This code segment intersects the top scoring hits of the screen with ligands having predicted activity,
#extracts corresponding rows from the ligand-receptor network, and generates a chord plot displaying the connections between ligands and receptors.

library(MAGeCKFlute)

#load screening results
SPH123_paired.gene_summary <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files/SPH123_paired.gene_summary.txt")
gdata = ReadRRA(SPH123_paired.gene_summary)
View(gdata)

#calculate top and bottom deciles (top-scoring hits)
quantile(gdata$Score, probs = seq(0, 1, 1/10))
nrow(gdata)
topdecile <- which(gdata$Score > 0.18)
bottomdecile <- which(gdata$Score < -0.19)
topscoringhits <- gdata[c(topdecile, bottomdecile),]
View(topscoringhits)
write.csv(topscoringhits, file="topscoringhits.csv")

#intersect top scoring hits of screen with LRs with predicted activity 
ligand_activities #109
overlap_active <- intersect(toupper(topscoringhits$id), ligand_activities$test_ligand)
overlap_active
[1] "PLXNB2"  "NECTIN2" "PSAP"    "APOE"    "MANF"    "NENF"    "FGA"     "PSEN1"   "APP"     "VTN"    
[11] "LTB"     "SAA1"    "HAMP"    "FAT1"    "HSPG2"   "GDF15"   "A2M"     "EFNA1"   "COL18A1" "PCDH1"  
[21] "BMP4"    "LIF"

# Extract rows from the ligand-receptor network containing overlapping active ligands
overlap_links_rows<-which(lr_network_expressed$from %in% overlap_active)
overlap_links <- lr_network_expressed[overlap_links_rows,]
View(overlap_links)
circos_links = overlap_links

#make chord plot 
library(circlize)
circos.par(gap.degree = gaps)
pdf(file="circos_plot.pdf",  width=4, height=4)
chordDiagram(circos_links, annotationTrack = c("name","grid"))
dev.off()


