#This script identifies potential interactions between screen hits and differentially expressed between liver mets and primary tumors by:

#Data Preprocessing: It loads scRNA-seq datasets from Che et al and Wang et al. Then, it preprocesses the data using Seurat, 
#including steps like filtering cells based on quality control metrics, normalization, dimensionality reduction, clustering, and UMAP visualization.

#Differential Expression Analysis: After identifying epithelial cells, it performs DEA between liver metastases and primary CRC tumors for both datasets 
#using the FindMarkers function in Seurat. It saves the results to CSV files.

#Querying External Databases: It queries several external databases (CellPhoneDB, CellTalkDB, NicheNet) to find interactors of differentially expressed genes (DEGs). 
#It then processes and merges the results from these databases.

#Intersecting Results: It intersects the results obtained from querying the external databases with top scoring hits of the screen to identify potential interactions.


####### CALCULATE DEREGULATED GENES BETWEEN LIVER METS AND MATCHED PRIMARY CRC TUMORS IN TWO INDEPENDENT DATASETS ######
library(dplyr)
library(Seurat)
library(sctransform)
library(tibble)
library(rstatix)
library(stringr)
library(MetBrewer)


###CHE ET AL####
#Load the Che et al dataset
Che <- Read10X(data.dir = "Published_datasets/Che_liver_mets")
Che <- CreateSeuratObject(counts = Che, project = "Che", min.cells = 3, min.features = 200,  meta.data = "matrix.mtx.gz")

#preprocess Seurat object
Che[["percent.mt"]] <- PercentageFeatureSet(Che, pattern = "^MT-") 
VlnPlot(Che, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Che, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Che, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(Che, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Che <- subset(Che, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 50)
Che <- NormalizeData(Che, normalization.method = "LogNormalize", scale.factor = 10000)
Che <- FindVariableFeatures(Che, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Che), 10)
Che <- ScaleData(Che, features = all.genes, vars.to.regress = "percent.mt")
Che <- RunPCA(Che, features = VariableFeatures(object = Che))
DimPlot(Che, reduction = "pca")
ElbowPlot(Che)
Che <- FindNeighbors(Che, dims = 1:20)
Che <- FindClusters(Che, resolution = 0.3)
Che <- RunUMAP(Che, dims = 1:20)

#identify and subset epithelial cells, then process again 
Idents(Che) <- "seurat_clusters"
DimPlot(Che, group.by = "orig.ident")
Idents(Che) <- "orig.ident"
DimPlot(Che, label=T)
FeaturePlot(Che, features=c("EPCAM", "CDH1"), order=T)
epithelial <- subset(Che, idents =c(8,10,11))
epithelial <- NormalizeData(epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
epithelial <- FindVariableFeatures(epithelial, selection.method = "vst", nfeatures = 2000)
epithelial <- ScaleData(epithelial, features = all.genes, vars.to.regress = "percent.mt")
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
DimPlot(epithelial, reduction = "pca")
ElbowPlot(epithelial)
epithelial <- FindNeighbors(epithelial, dims = 1:20)
epithelial <- FindClusters(epithelial, resolution = 0.3)
epithelial <- RunUMAP(epithelial, dims = 1:20)
DimPlot(epithelial)

#Identify sample origin (patient) and location (CRC primary vs liver mets)
#Extract barcodes from row names
barcodes<- rownames(epithelial@meta.data)
head(barcodes)

# Split barcodes to extract sample information (patient and location)
sample <- str_split(barcodes, pattern="_", simplify =T)[,2:3]
unique(sample)
epithelial$patient<- str_split(rownames(epithelial@meta.data), pattern="_", simplify =T)[,2]
DimPlot(epithelial, group.by="patient")

epithelial$location<- str_split(rownames(epithelial@meta.data), pattern="_", simplify =T)[,3]
DimPlot(epithelial, group.by="location")

# Perform differential expression analysis between liver mets (LM) and CRC primary (CRC)
Idents(epithelial) <- "location"
mets.vs.primary <- FindMarkers(epithelial, ident.1 = "LM", ident.2 = "CRC")
write.csv(mets.vs.primary, file = "mets.vs.primary.Che.csv")


###WANG ET AL####
#Load the Wang et al dataset

counts <- read.delim("/Published_datasets/Wang_liver_mets/GSM7058755_non_immune_counts.txt.gz",row.names=1, header=T)
head(counts[1:5,1:5])
GSM7058755_non_immune_meta <- read.delim("/media/Coco/MOSAIC LIVER/Published_datasets/Wang_liver_mets/GSM7058755_non_immune_meta.txt.gz", row.names=1)

#Transpose meta data as Seurat expects meta data to have cell names as rows and meta data values as columns
Wang <- CreateSeuratObject(counts = counts, project = "Wang", meta.data = data.frame(t(GSM7058755_non_immune_meta))) 
Wang <- AddMetaData(object = Wang, metadata = GSM7058755_non_immune_meta$organ, col.name = "organ")
Wang$organ

#Processing
Wang[["percent.mt"]] <- PercentageFeatureSet(Wang, pattern = "^MT-") 
VlnPlot(Wang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Wang, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Wang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(Wang, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Wang <- subset(Wang, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 50)
Wang <- NormalizeData(Wang, normalization.method = "LogNormalize", scale.factor = 10000)
Wang <- FindVariableFeatures(Wang, selection.method = "vst", nfeatures = 2000)
Wang <- ScaleData(Wang, features = all.genes, vars.to.regress = "percent.mt")
Wang <- RunPCA(Wang, features = VariableFeatures(object = Wang))
DimPlot(Wang, reduction = "pca")
ElbowPlot(Wang)
Wang <- FindNeighbors(Wang, dims = 1:20)
Wang <- FindClusters(Wang, resolution = 0.3)
Wang <- RunUMAP(Wang, dims = 1:20)

#identify and subset epithelial cells, then process again 
Idents(Wang) <- "seurat_clusters"
FeaturePlot(Wang, features=c("EPCAM", "CDH1"), order=T)
epithelial <- subset(Wang, idents =c(0, 2, 5,7,8))
epithelial <- NormalizeData(epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
epithelial <- FindVariableFeatures(epithelial, selection.method = "vst", nfeatures = 2000)
epithelial <- ScaleData(epithelial,  features = all.genes, vars.to.regress = "percent.mt")
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
DimPlot(epithelial, reduction = "pca")
ElbowPlot(epithelial)
epithelial <- FindNeighbors(epithelial, dims = 1:20)
epithelial <- FindClusters(epithelial, resolution = 0.3)
epithelial <- RunUMAP(epithelial, dims = 1:20)
DimPlot(epithelial, group.by = "organ")

# Perform differential expression analysis between liver mets (LCT) and CRC primary (CCT)
Idents(epithelial) <- "organ"
LCT_vs_CCT <- FindMarkers(epithelial, ident.1 = "LCT", ident.2 = "CCT")
View(LCT_vs_CCT)
write.csv(LCT_vs_CCT, "livermetsvsprimary.Wang.csv")


##########FOR BOTH THE CHE AND WANG DATASET, QUERY CELLTALK BD, CELLPHONE DB AND NICHENET DATABASES FOR INTERACTORS (code from Elena Guido Vinzoni######
#below is the code used for the the Che dataset, the Wang dataset was processed the same way 
library(openxlsx)
library(tidyverse)
library(dplyr)
library(stringr)

#load list of degs
che <- read.csv("mets.vs.primary.Che.csv")

# Filter for the DEGs between liver mets and primary
posDEG <- che %>% filter(`avg_log2FC` > 0)
negDEG <- che %>% filter(`avg_log2FC` < 0)

# Save
write.xlsx(posDEG, "posDEGs_liver.xlsx")
write.xlsx(negDEG, "negDEGs_liver.xlsx")

setwd("/media/Elena/Just an idea/Che_scRNAseq//data")
posDEG <- read.xlsx("posDEGs_liver.xlsx")
negDEG <- read.xlsx("negDEGs_liver.xlsx")


####CELLPHONE DB
# make txt files to query CPDB
pos_names <- NULL
neg_names <- NULL
pos_names$names <- posDEG$X
neg_names$names <- negDEG$X
writeLines(pos_names$names, "CPDB_query_pos.txt")
writeLines(neg_names$names, "CPDB_query_neg.txt")

# each is a bit long 
indices <- seq_along(pos_names[[1]])
group_factor <- cut(indices, breaks=4, labels=FALSE)
split_lists <- split(pos_names[[1]], group_factor)

for (i in seq_along(split_lists)) {
  file_name <- sprintf("CPDB_query_pos_%d.txt", i)
  writeLines(split_lists[[i]], con=file_name)
}

#now for the negative DEGs
indices <- seq_along(neg_names[[1]])
group_factor <- cut(indices, breaks=6, labels=FALSE)
split_lists <- split(neg_names[[1]], group_factor)

for (i in seq_along(split_lists)) {
  file_name <- sprintf("CPDB_query_neg_%d.txt", i)
  writeLines(split_lists[[i]], con=file_name)
}

# Now query the genes on https://www.cellphonedb.org/
# The queried results are 4 files for pos and 6 for neg
# merge the files together
pos_1 <- read.csv("cpdb_search_results_pos_1.03_04_2024.10_53_54.csv")
pos_2 <- read.csv("cpdb_search_results_pos_2.03_04_2024.10_54_50.csv")
pos_3 <- read.csv("cpdb_search_results_pos_3.03_04_2024.10_55_16.csv")
pos_4 <- read.csv("cpdb_search_results_pos_4.03_04_2024.10_55_51.csv")
neg_1 <- read.csv("cpdb_search_results_neg_1.03_04_2024.10_56_26.csv")
neg_2 <- read.csv("cpdb_search_results_neg_2.03_04_2024.10_56_49.csv")
neg_3 <- read.csv("cpdb_search_results_neg_3.03_04_2024.10_57_14.csv")
neg_4 <- read.csv("cpdb_search_results_neg_4.03_04_2024.10_57_41.csv")
neg_5 <- read.csv("cpdb_search_results_neg_5.03_04_2024.10_58_03.csv")
neg_6 <- read.csv("cpdb_search_results_neg_6.03_04_2024.10_58_30.csv")

list_of_pos <- list(pos_1, pos_2, pos_3, pos_4) 
pos_result <- do.call(rbind, list_of_pos)

list_of_neg <- list(neg_1, neg_2, neg_3, neg_4, neg_5, neg_6) 
neg_result <- do.call(rbind, list_of_neg)

pos_result$row_number <- seq_len(nrow(pos_result))
neg_result$row_number <- seq_len(nrow(neg_result))

# Fill empty gene names
# Filter rows where Gene name B is empty
filtered_data <- pos_result %>% filter(is.na(`Gene.name.B`) | `Gene.name.B` == "")
# Extract values between 'complex:' and '_complex' in Partner B column
filtered_data$Partner_B_extracted <- str_extract(filtered_data$`Partner.B`, "(?<=complex:)(.+?)(?=_complex$|_receptor$|$)")
# Replace empty spaces in Gene name B column with extracted values
filtered_data$`Gene.name.B` <- filtered_data$Partner_B_extracted
filtered_pos_A <- filtered_data 

# Filter rows where Gene name B is empty
empty_rows <- which(is.na(pos_result$`Gene.name.B`) | pos_result$`Gene.name.B` == "")
# Replace empty Gene name B cells with filtered genes using row indices
pos_result$`Gene.name.B`[empty_rows] <- filtered_pos_A$`Gene.name.B`
pos_result_filled <- pos_result 

# Filter rows where Gene name A is empty
filtered_data <- pos_result %>% filter(is.na(`Gene.name.A`) | `Gene.name.A` == "")
# Extract values between 'complex:' and '_complex' in Partner B column
filtered_data$Partner_A_extracted <- str_extract(filtered_data$`Partner.A`, "(?<=complex:)(.+?)(?=_complex$|_receptor$|$)")
# Replace empty spaces in Gene name B column with extracted values
filtered_data$`Gene.name.A` <- filtered_data$Partner_A_extracted
filtered_pos_B <- filtered_data 

# Filter rows where Gene name A is empty
empty_rows <- which(is.na(pos_result$`Gene.name.A`) | pos_result$`Gene.name.A` == "")
# Replace empty Gene name A cells with filtered genes using row indices
pos_result$`Gene.name.A`[empty_rows] <- filtered_pos_B$`Gene.name.A`
pos_result_filled <- pos_result


# Same for the negatively regulated genes
# Filter rows where Gene name B is empty
filtered_data <- neg_result %>% filter(is.na(`Gene.name.B`) | `Gene.name.B` == "")
# Extract values between 'complex:' and '_complex' in Partner B column
filtered_data$Partner_B_extracted <- str_extract(filtered_data$`Partner.B`, "(?<=complex:)(.+?)(?=_complex$|_receptor$|$)")
# Replace empty spaces in Gene name B column with extracted values
filtered_data$`Gene.name.B` <- filtered_data$Partner_B_extracted
filtered_neg_A <- filtered_data 

# Filter rows where Gene name B is empty
empty_rows <- which(is.na(neg_result$`Gene.name.B`) | neg_result$`Gene.name.B` == "")
# Replace empty Gene name B cells with filtered genes using row indices
neg_result$`Gene.name.B`[empty_rows] <- filtered_neg_A$`Gene.name.B`
neg_result_filled <- neg_result 

# Filter rows where Gene name A is empty
filtered_data <- neg_result %>% filter(is.na(`Gene.name.A`) | `Gene.name.A` == "")
# Extract values between 'complex:' and '_complex' in Partner B column
filtered_data$Partner_A_extracted <- str_extract(filtered_data$`Partner.A`, "(?<=complex:)(.+?)(?=_complex$|_receptor$|$)")
# Replace empty spaces in Gene name B column with extracted values
filtered_data$`Gene.name.A` <- filtered_data$Partner_A_extracted
filtered_neg_B <- filtered_data 

# Filter rows where Gene name A is empty
empty_rows <- which(is.na(neg_result$`Gene.name.A`) | neg_result$`Gene.name.A` == "")
# Replace empty Gene name A cells with filtered genes using row indices
neg_result$`Gene.name.A`[empty_rows] <- filtered_neg_B$`Gene.name.A`
neg_result_filled <- neg_result

#####
# Filter rows with targets in Gene name A col and select Gene name B and Partner B
pos_DEGs_as_geneA <- pos_result_filled %>%
  rowwise() %>% 
  filter(any(sapply(pos_names$names, function(x) str_detect(`Gene.name.A`, x) | str_detect(`Partner.A`, x)))) %>%
  ungroup()

pos_DEGs_as_geneB <- pos_result_filled %>%
  rowwise() %>% 
  filter(any(sapply(pos_names$names, function(x) str_detect(`Gene.name.B`, x) | str_detect(`Partner.B`, x)))) %>%
  ungroup()

neg_DEGs_as_geneA <- neg_result_filled %>%
  rowwise() %>% 
  filter(any(sapply(neg_names$names, function(x) str_detect(`Gene.name.A`, x) | str_detect(`Partner.A`, x)))) %>%
  ungroup()

neg_DEGs_as_geneB <- neg_result_filled %>%
  rowwise() %>% 
  filter(any(sapply(neg_names$names, function(x) str_detect(`Gene.name.B`, x) | str_detect(`Partner.B`, x)))) %>%
  ungroup()

#Save the lists
write.xlsx(pos_DEGs_as_geneA, "CPDB_posDEGs_as_geneA.xlsx")
write.xlsx(pos_DEGs_as_geneB, "CPDB_posDEGs_as_geneB.xlsx")
write.xlsx(neg_DEGs_as_geneA, "CPDB_negDEGs_as_geneA.xlsx")
write.xlsx(neg_DEGs_as_geneB, "CPDB_negDEGs_as_geneB.xlsx")

#merge them
#we make a new dataframe with a column of the cancer posDEGs and the other are the interactors

pos_A_select <- pos_DEGs_as_geneA %>%
  select(Gene.name.A, Gene.name.B) %>%
  rename(pos_DEGs = Gene.name.A, interactor = Gene.name.B)
pos_B_select <- pos_DEGs_as_geneB %>%
  select(Gene.name.B, Gene.name.A) %>%
  rename(pos_DEGs = Gene.name.B, interactor = Gene.name.A)

merged_posDEGs <- rbind(pos_A_select, pos_B_select)

write.xlsx(merged_posDEGs, "CPDB_posDEGs.xlsx")


neg_A_select <- neg_DEGs_as_geneA %>%
  select(Gene.name.A, Gene.name.B) %>%
  rename(neg_DEGs = Gene.name.A, interactor = Gene.name.B)
neg_B_select <- neg_DEGs_as_geneB %>%
  select(Gene.name.B, Gene.name.A) %>%
  rename(neg_DEGs = Gene.name.B, interactor = Gene.name.A)

merged_negDEGs <- rbind(neg_A_select, neg_B_select)

write.xlsx(merged_negDEGs, "CPDB_negDEGs.xlsx")

# Now we can intersect 
#first we just see which genes that pop up are part of the top deciles:
# top deciles separated into proximal and distal
hits <- read.csv("topscoringhits_human.csv")
phits = hits %>%
  select(screen_proximal = proximal) %>% filter(!is.na(screen_proximal), screen_proximal != "")
dhits = hits %>%
  select(screen_distal = distal) %>% filter(!is.na(screen_distal), screen_distal != "")

# Intersect pos DEGs and proximal hits 
intersect_pos_prox <- merged_posDEGs %>%
  filter(`interactor` %in% phits$screen_proximal) %>%
  select(`interactor`, `pos_DEGs`)
# Intersect pos DEGs and distal hits
intersect_pos_dist <- merged_posDEGs %>%
  filter(`interactor` %in% dhits$screen_distal) %>%
  select(`interactor`, `pos_DEGs`)

# Intersect neg DEGs and proximal hits
intersect_neg_prox <- merged_negDEGs %>%
  filter(`interactor` %in% phits$screen_proximal) %>%
  select(`interactor`, `neg_DEGs`)
# Intersect neg DEGs and distal hits
intersect_neg_dist <- merged_negDEGs %>%
  filter(`interactor` %in% dhits$screen_distal) %>%
  select(`interactor`, `neg_DEGs`)


# Save
write.xlsx(intersect_pos_prox, "CPDB_intersect_pos_prox.xlsx")
write.xlsx(intersect_pos_dist, "CPDB_intersect_pos_dist.xlsx")
write.xlsx(intersect_neg_prox, "CPDB_intersect_neg_prox.xlsx")
write.xlsx(intersect_neg_dist, "CPDB_intersect_neg_dist.xlsx")

####CELLTALK DB
CellTalkDBInteractions <- read.delim("human_lr_pair_CellTalkDB.txt", header = TRUE, sep = "\t")

#get the list of positive and negative DEGs
posDEG <- read.xlsx("posDEGs_liver.xlsx")
negDEG <- read.xlsx("negDEGs_liver.xlsx")
pos_names <- NULL
neg_names <- NULL
pos_names$names <- posDEG$X
neg_names$names <- negDEG$X

#get the list of pre available LRs of the screen hits (top decile)
hits_CTDB_prox_lig <- read.xlsx("CellTalkDB_interactors_top60prox.xlsx")
hits_CTDB_dist_lig <- read.xlsx("CellTalkDB_interactors_top60dist.xlsx")

file_path <- "CellTalkDB_interactors_top60prox.xlsx"
hits_CTDB_prox_rec <- read.xlsx(file_path, sheet = 2)
file_path <- "CellTalkDB_interactors_top60dist.xlsx"
hits_CTDB_dist_rec <- read.xlsx(file_path, sheet = 2)

hits_CTDB_prox <- hits_CTDB_prox_lig %>%
  select(receptor_gene_symbol, ligand_list) %>%
  rename(interactor = receptor_gene_symbol, prox_hit = ligand_list)
hits_CTDB_prox_R <- hits_CTDB_prox_rec %>%
  select(ligand_gene_symbol, receptor_list) %>%
  rename(interactor = ligand_gene_symbol, prox_hit = receptor_list)

hits_CTDB_prox <- rbind(hits_CTDB_prox, hits_CTDB_prox_R)

write.xlsx(hits_CTDB_prox, "CTDB_prox_merged.xlsx")

hits_CTDB_dist <- hits_CTDB_dist_lig %>%
  select(receptor_gene_symbol, ligand_list) %>%
  rename(interactor = receptor_gene_symbol, dist_hit = ligand_list)
hits_CTDB_dist_R <- hits_CTDB_dist_rec %>%
  select(ligand_gene_symbol, receptor_list) %>%
  rename(interactor = ligand_gene_symbol, dist_hit = receptor_list)

hits_CTDB_dist <- rbind(hits_CTDB_dist, hits_CTDB_dist_R)

write.xlsx(hits_CTDB_dist, "CTDB_dist_merged.xlsx")

# Intersect proximal hits and pos DEGs
intersect_prox_pos <- hits_CTDB_prox %>%
  filter(`interactor` %in% pos_names$names) %>%
  select(`interactor`, `prox_hit`)
# same result as the intersection the other way around TREM2-APOE
intersect_prox_pos
# Intersect distal hits and pos DEGs 
intersect_dist_pos <- hits_CTDB_dist %>%
  filter(`interactor` %in% pos_names$names) %>%
  select(`interactor`, `dist_hit`)
# similar
intersect_dist_pos
# Intersect proximal hits and neg DEGs 
intersect_prox_neg <- hits_CTDB_prox %>%
  filter(`interactor` %in% neg_names$names) %>%
  select(`interactor`, `prox_hit`)
# similar
intersect_prox_neg
# Intersect distal hits and neg DEGs
intersect_dist_neg <- hits_CTDB_dist %>%
  filter(`interactor` %in% neg_names$names) %>%
  select(`interactor`, `dist_hit`)
# similar
intersect_dist_neg

#Save
write.xlsx(intersect_prox_pos, "CTDB_intersect_prox_pos.xlsx")
write.xlsx(intersect_dist_pos, "CTDB_intersect_dist_pos.xlsx")
write.xlsx(intersect_prox_neg, "CTDB_intersect_prox_neg.xlsx")
write.xlsx(intersect_dist_neg, "CTDB_intersect_dist_neg.xlsx")

##########NICHE NET
posDEG <- read.xlsx("posDEGs_liver.xlsx")
negDEG <- read.xlsx("negDEGs_liver.xlsx")
pos_names <- NULL
neg_names <- NULL
pos_names$names <- posDEG$X
neg_names$names <- negDEG$X
#get the list of pre available LRs of the screen hits (top decile)


hits_NN_prox_lig <- read.xlsx("NN_interactors_top60prox.xlsx")
hits_NN_dist_lig <- read.xlsx("NN_interactors_top60dist.xlsx")
file_path <- "NN_interactors_top60prox.xlsx"
hits_NN_prox_rec <- read.xlsx(file_path, sheet = 2)
file_path <- "NN_interactors_top60dist.xlsx"
hits_NN_dist_rec <- read.xlsx(file_path, sheet = 2)

hits_NN_prox <- hits_NN_prox_lig %>%
  select(to, ligand_list) %>%
  rename(interactor = to, prox_hit = ligand_list)
hits_NN_prox_R <- hits_NN_prox_rec %>%
  select(from, receptor_list) %>%
  rename(interactor = from, prox_hit = receptor_list)

hits_NN_prox <- rbind(hits_NN_prox, hits_NN_prox_R)

write.xlsx(hits_NN_prox, "NN_prox_merged.xlsx")

hits_NN_dist <- hits_NN_dist_lig %>%
  select(to, ligand_list) %>%
  rename(interactor = to, dist_hit = ligand_list)
hits_NN_dist_R <- hits_NN_dist_rec %>%
  select(from, receptor_list) %>%
  rename(interactor = from, dist_hit = receptor_list)

hits_NN_dist <- rbind(hits_NN_dist, hits_NN_dist_R)

write.xlsx(hits_NN_dist, "NN_dist_merged.xlsx")

# Intersect proximal hits and pos DEGs
intersect_prox_pos <- hits_NN_prox %>%
  filter(`interactor` %in% pos_names$names) %>%
  select(`interactor`, `prox_hit`)
# same result as the intersection the other way around TREM2-APOE
intersect_prox_pos
# Intersect distal hits and pos DEGs 
intersect_dist_pos <- hits_NN_dist %>%
  filter(`interactor` %in% pos_names$names) %>%
  select(`interactor`, `dist_hit`)
# similar
intersect_dist_pos
# Intersect proximal hits and neg DEGs 
intersect_prox_neg <- hits_NN_prox %>%
  filter(`interactor` %in% neg_names$names) %>%
  select(`interactor`, `prox_hit`)
# similar
intersect_prox_neg
# Intersect distal hits and neg DEGs
intersect_dist_neg <- hits_NN_dist %>%
  filter(`interactor` %in% neg_names$names) %>%
  select(`interactor`, `dist_hit`)
# similar
intersect_dist_neg

#Save
write.xlsx(intersect_prox_pos, "NN_intersect_prox_pos.xlsx")
write.xlsx(intersect_dist_pos, "NN_intersect_dist_pos.xlsx")
write.xlsx(intersect_prox_neg, "NN_intersect_prox_neg.xlsx")
write.xlsx(intersect_dist_neg, "NN_intersect_dist_neg.xlsx")



######### MERGE RESULTS FROM ALL THREE DATABASES########
#start with pos-prox
NN_pos_prox <- read.xlsx("NN_intersect_pos_prox.xlsx")
NN_prox_pos <- read.xlsx("NN_intersect_prox_pos.xlsx")
CTDB_pos_prox <- read.xlsx("CTDB_intersect_pos_prox.xlsx")
CTDB_prox_pos <- read.xlsx("CTDB_intersect_prox_pos.xlsx")
CPDB_pos_prox <- read.xlsx("CPDB_intersect_pos_prox.xlsx")
# CPDB_prox_pos <- read.xlsx("CPDB_intersect_prox_pos.xlsx") does not exist (yet?)

# First we merge the ones that have the same structure:
all_pos_prox <- rbind(NN_pos_prox,CTDB_pos_prox )
all_pos_prox <- rbind(all_pos_prox,CPDB_pos_prox)

all_prox_pos <- rbind(NN_prox_pos,CTDB_prox_pos)
all_prox_pos$prox_hit <- str_split(all_prox_pos$prox_hit, pattern = ",\\s*")
all_prox_pos <- all_prox_pos %>% unnest(prox_hit)

unique_all_pos_prox <- all_pos_prox %>% distinct()
unique_all_prox_pos <- all_prox_pos %>% distinct()

#merge
unique_all_pos_prox <- unique_all_pos_prox %>%
  select(interactor, pos_DEGs) %>%
  rename(prox_hit = interactor, pos_DEGs = pos_DEGs)
unique_all_prox_pos <- unique_all_prox_pos %>%
  select(interactor, prox_hit) %>%
  rename(pos_DEGs = interactor, prox_hit = prox_hit)

all_pp <- rbind(unique_all_pos_prox, unique_all_prox_pos)
all_pp <- all_pp %>% distinct()

#save
write.xlsx(all_pp, "all_pp.xlsx")


#second with with pos-dist
NN_pos_dist <- read.xlsx("NN_intersect_pos_dist.xlsx")
NN_dist_pos <- read.xlsx("NN_intersect_dist_pos.xlsx")
CTDB_pos_dist <- read.xlsx("CTDB_intersect_pos_dist.xlsx")
CTDB_dist_pos <- read.xlsx("CTDB_intersect_dist_pos.xlsx")
CPDB_pos_dist <- read.xlsx("CPDB_intersect_pos_dist.xlsx")
# CPDB_dist_pos <- read.xlsx("CPDB_intersect_dist_pos.xlsx") does not exist (yet?)

# First we merge the ones that have the same structure:
all_pos_dist <- rbind(NN_pos_dist,CTDB_pos_dist )
all_pos_dist <- rbind(all_pos_dist,CPDB_pos_dist)

all_dist_pos <- rbind(NN_dist_pos,CTDB_dist_pos)
all_dist_pos$dist_hit <- str_split(all_dist_pos$dist_hit, pattern = ",\\s*")
all_dist_pos <- all_dist_pos %>% unnest(dist_hit)

unique_all_pos_dist <- all_pos_dist %>% distinct()
unique_all_dist_pos <- all_dist_pos %>% distinct()

#merge
unique_all_pos_dist <- unique_all_pos_dist %>%
  select(interactor, pos_DEGs) %>%
  rename(dist_hit = interactor, pos_DEGs = pos_DEGs)
unique_all_dist_pos <- unique_all_dist_pos %>%
  select(interactor, dist_hit) %>%
  rename(pos_DEGs = interactor, dist_hit = dist_hit)

all_pd <- rbind(unique_all_pos_dist, unique_all_dist_pos)
all_pd <- all_pd %>% distinct()
#save
write.xlsx(all_pd, "all_pd.xlsx")


#third with with neg-dist
NN_neg_dist <- read.xlsx("NN_intersect_neg_dist.xlsx")
NN_dist_neg <- read.xlsx("NN_intersect_dist_neg.xlsx")
CTDB_neg_dist <- read.xlsx("CTDB_intersect_neg_dist.xlsx")
CTDB_dist_neg <- read.xlsx("CTDB_intersect_dist_neg.xlsx")
CPDB_neg_dist <- read.xlsx("CPDB_intersect_neg_dist.xlsx")
# CPDB_dist_neg <- read.xlsx("CPDB_intersect_dist_neg.xlsx") does not exist (yet?)

# First we merge the ones that have the same structure:
all_neg_dist <- rbind(NN_neg_dist,CTDB_neg_dist )
all_neg_dist <- rbind(all_neg_dist,CPDB_neg_dist)

all_dist_neg <- rbind(NN_dist_neg,CTDB_dist_neg)
all_dist_neg$dist_hit <- str_split(all_dist_neg$dist_hit, pattern = ",\\s*")
all_dist_neg <- all_dist_neg %>% unnest(dist_hit)

unique_all_neg_dist <- all_neg_dist %>% distinct()
unique_all_dist_neg <- all_dist_neg %>% distinct()

#merge
unique_all_neg_dist <- unique_all_neg_dist %>%
  select(interactor, neg_DEGs) %>%
  rename(dist_hit = interactor, neg_DEGs = neg_DEGs)
unique_all_dist_neg <- unique_all_dist_neg %>%
  select(interactor, dist_hit) %>%
  rename(neg_DEGs = interactor, dist_hit = dist_hit)

all_nd <- rbind(unique_all_neg_dist, unique_all_dist_neg)
all_nd <- all_nd %>% distinct()

#save
write.xlsx(all_nd, "all_nd.xlsx")

                    
#finally with with neg-prox
NN_neg_prox <- read.xlsx("NN_intersect_neg_prox.xlsx")
NN_prox_neg <- read.xlsx("NN_intersect_prox_neg.xlsx")
CTDB_neg_prox <- read.xlsx("CTDB_intersect_neg_prox.xlsx")
CTDB_prox_neg <- read.xlsx("CTDB_intersect_prox_neg.xlsx")
CPDB_neg_prox <- read.xlsx("CPDB_intersect_neg_prox.xlsx")

# First we merge the ones that have the same structure:
all_neg_prox <- rbind(NN_neg_prox,CTDB_neg_prox )
all_neg_prox <- rbind(all_neg_prox,CPDB_neg_prox)

all_prox_neg <- rbind(NN_prox_neg,CTDB_prox_neg)
all_prox_neg$prox_hit <- str_split(all_prox_neg$prox_hit, pattern = ",\\s*")
all_prox_neg <- all_prox_neg %>% unnest(prox_hit)

unique_all_neg_prox <- all_neg_prox %>% distinct()
unique_all_prox_neg <- all_prox_neg %>% distinct()

#merge
unique_all_neg_prox <- unique_all_neg_prox %>%
  select(interactor, neg_DEGs) %>%
  rename(prox_hit = interactor, neg_DEGs = neg_DEGs)
unique_all_prox_neg <- unique_all_prox_neg %>%
  select(interactor, prox_hit) %>%
  rename(neg_DEGs = interactor, prox_hit = prox_hit)

all_np <- rbind(unique_all_neg_prox, unique_all_prox_neg)
all_np <- all_np %>% distinct()
#save
write.xlsx(all_np, "all_np.xlsx")


#use circlize package to make chord plots in Extended Data Fig2h (one for proximal hits and all deregulated genes, one for distal hits)



                    
