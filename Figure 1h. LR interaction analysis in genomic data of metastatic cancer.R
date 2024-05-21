#Code by Elena Guido Vinzoni to identify LR interactions between screening hits and frequenctly mutated LRs in hepatic metastases
#The code below exemplifies the workflow using the CellTalk database: 

#1. extracting interactors of screening hits quering the CellTalk database (the NicheNet and CellPhone DB databases were queried in the same way)
#The code focuses first on the proximal hits, identifying ligands and receptors within these hits and their interactors.
# Interactors and their frequency are saved in Excel files.
# The same process is repeated for distal top 60 hits.

#2. intersecting with mutational data from Nguyen et al (s4bliver_unique)
# Following the extraction of LR interactions, the code intersects the results with mutational data from Nguyen et al.
# Mutated genes are extracted and compared with LR interactors.
# Exact matches are found, and the results are saved in Excel files, separately for proximal and distal LR interactions.

# Finally, for each proximal and distal LR interaction, the code identifies gene amplifications and deletions.
# The results are saved in Excel files for further analysis and visualization.

#results from the analyses with Celltalk DB, Cellphone DB and NicheNet were then manually merged and plotted in R using the ggalluvial package to generate Fig 1h

#######EXTRACTING INTERACTORS OF SCREENING HITS FROM CELLTALK DB##########
library(openxlsx)
library(tidyverse)

#finding genes from coco's hits in CellTalkDB
top_60_hits <- read.csv("topscoringhits_human.csv", header = T)
CellTalkDBInteractions <- read.delim("human_lr_pair_CellTalkDB.txt", header = TRUE, sep = "\t")
library(tidyverse)

#PROXIMAL TOP60
# find ligands and receptors within hits amongst the top proximal hits
ligands = CellTalkDBInteractions %>% pull(ligand_gene_symbol) %>% unique()
hits_ligands = intersect(ligands,top_60_hits$proximal)

receptors = CellTalkDBInteractions %>% pull(receptor_gene_symbol) %>% unique()
hits_receptors = intersect(receptors,top_60_hits$proximal)

# we then got the ligands amongst hits and checked for their interactors
lr_network_hits_ligands = CellTalkDBInteractions %>% dplyr::filter(ligand_gene_symbol %in% hits_ligands)
#  then got the receptors amongst hits and checked for their interactors
lr_network_hits_receptors = CellTalkDBInteractions %>% dplyr::filter(receptor_gene_symbol %in% hits_receptors)

# save the files 
write.xlsx(lr_network_hits_ligands, "CellTalkDB_top60_prox_ligands.xlsx")
write.xlsx(lr_network_hits_receptors, "CellTalkDB_top60_prox_receptors.xlsx")

# look for unique interactors and their frquency 
interactors = c(lr_network_hits_ligands$receptor_gene_symbol, lr_network_hits_receptors$ligand_gene_symbol)
freq_interactors = as.data.frame(table(interactors))

# save the files 
write.xlsx(freq_interactors, "CellTalkDB_top60_prox_interactors.xlsx")

# for each interactor which genes it interacts with
library(dplyr)
hits_as_ligands <- lr_network_hits_ligands %>%
  group_by(receptor_gene_symbol) %>%
  summarise(ligand_list = paste(ligand_gene_symbol, collapse = ", ")) %>%
  ungroup()

hits_as_receptors <- lr_network_hits_receptors %>%
  group_by(ligand_gene_symbol) %>%
  summarise(receptor_list = paste(receptor_gene_symbol, collapse = ", ")) %>%
  ungroup()

wb <- createWorkbook()
# Add data frames as sheets to the Excel workbook
addWorksheet(wb, "coco_top60prox_as_ligands")
writeData(wb, sheet = 1, hits_as_ligands, startCol = 1, startRow = 1)

addWorksheet(wb, "coco_top60prox_as_receptors")
writeData(wb, sheet = 2, hits_as_receptors, startCol = 1, startRow = 1)

# Save the Excel workbook to a file
saveWorkbook(wb, "CellTalkDB_interactors_top60prox.xlsx", overwrite = TRUE)


#DISTAL TOP60
# find ligands and receptors withing Coco's hits amongst the top proximal hits
ligands = CellTalkDBInteractions %>% pull(ligand_gene_symbol) %>% unique()
hits_ligands = intersect(ligands,top_60_hits$distal)

receptors = CellTalkDBInteractions %>% pull(receptor_gene_symbol) %>% unique()
hits_receptors = intersect(receptors,top_60_hits$distal)

# we then got the ligands that are in Cocos hits and checked for their interactors
lr_network_hits_ligands = CellTalkDBInteractions %>% dplyr::filter(ligand_gene_symbol %in% hits_ligands)
#  then got the receptors that are in Cocos hits and checked for their interactors
lr_network_hits_receptors = CellTalkDBInteractions %>% dplyr::filter(receptor_gene_symbol %in% hits_receptors)

# save the files 
write.xlsx(lr_network_hits_ligands, "CellTalkDB_top60_dist_ligands.xlsx")
write.xlsx(lr_network_hits_receptors, "CellTalkDB_top60_dist_receptors.xlsx")

# look for unique interactors and their frquency 
interactors = c(lr_network_hits_ligands$receptor_gene_symbol, lr_network_hits_receptors$ligand_gene_symbol)
freq_interactors = as.data.frame(table(interactors))

# save the files 
write.xlsx(freq_interactors, "CellTalkDB_top60_dist_interactors.xlsx")

# for each interactor which genes it interacts with

library(dplyr)
hits_as_ligands <- lr_network_hits_ligands %>%
  group_by(receptor_gene_symbol) %>%
  summarise(ligand_list = paste(ligand_gene_symbol, collapse = ", ")) %>%
  ungroup()

hits_as_receptors <- lr_network_hits_receptors %>%
  group_by(ligand_gene_symbol) %>%
  summarise(receptor_list = paste(receptor_gene_symbol, collapse = ", ")) %>%
  ungroup()

wb <- createWorkbook()
# Add data frames as sheets to the Excel workbook
addWorksheet(wb, "coco_top60dist_as_ligands")
writeData(wb, sheet = 1, hits_as_ligands, startCol = 1, startRow = 1)

addWorksheet(wb, "coco_top60dist_as_receptors")
writeData(wb, sheet = 2, hits_as_receptors, startCol = 1, startRow = 1)

# Save the Excel workbook to a file
saveWorkbook(wb, "CellTalkDB_interactors_top60dist.xlsx", overwrite = TRUE)



#######INTERSECTION WITH MUTATIONAL DATA########
#find intersect between screen-targets-interactors and mutational data from Nguyen et al (s4bliver_unique)
mutations <- read.xlsx("Nguyen_sb4_liver_unique_1.xlsx")

#Distal
mutated_genes <- gsub("^(.*)_.*$", "\\1", mutations$s4bliver_unique)
ctdb_interactors_dist <- read.xlsx("CellTalkDB_interactors_top60dist.xlsx", sheet=3)

# Find exact matches with mutated_genes
exact_matches <- intersect(mutated_genes, ctdb_interactors_dist$interactor)

# Filter matching genes based on exact matches
matching_genes <- as.data.frame(ctdb_interactors_dist[ctdb_interactors_dist$interactor %in% exact_matches, ])
colnames(matching_genes) <- c("exact_matches","coco_hit_dist")

#Save
write.xlsx(matching_genes, "top60_dist_LRs_in_s4b.xlsx")


#Proximal
ctdb_interactors_prox <- read.xlsx("CellTalkDB_interactors_top60prox.xlsx", sheet=3)
mutated_genes <- gsub("^(.*)_.*$", "\\1", mutations$s4bliver_unique)

# Find exact matches with mutated_genes
exact_matches <- intersect(mutated_genes, ctdb_interactors_prox$interactor)

# Filter matching genes based on exact matches
matching_genes <- as.data.frame(ctdb_interactors_prox[ctdb_interactors_prox$interactor %in% exact_matches, ])
colnames(matching_genes) <- c("exact_matches","coco_hit_prox")

#Save
write.xlsx(matching_genes, "top60_prox_LRs_in_s4b.xlsx")



#---------------------------------------------------
#Proximal
# Create an empty vector to store matched entries
matched_entries <- c()
# Loop through each gene in matching_genes$exact_matches
for (gene in matching_genes$exact_matches) {
  # Use grep to find exact matches (ignoring the description part after "_")
  matches <- grep(paste0("^", gene, "_"), mutations$s4bliver_unique, value = TRUE)
  # Append the matched entries to the vector
  matched_entries <- c(matched_entries, matches)
}
# Create a new dataframe with matched entries
matched_df <- as.data.frame(mutations[mutations$s4bliver_unique %in% matched_entries, ])
colnames(matched_df) <- "gene_name"



######CHECK WHICH GENES WERE AMPLIFIED AND WHICH DELETED#######
#now that we got all the entries in s4b that correspond to LRs of the targets we can divide them in amp and dels
# Use grep to filter entries ending with "_Amplification"
amplification_entries <- grep("^(.*?)_Amplification", matched_df$gene_name, value = TRUE)
# Extract gene names (remove "_Amplification" part)
gene_names <- sub("^(.*)_Amplification$", "\\1", amplification_entries)
# Create a new dataframe amp_df with extracted gene names
amp_df <- data.frame(gene_name = gene_names)
colnames(amp_df) <- "Amplifications"

# Use grep to filter entries ending with "_Deletion"
deletion_entries <- grep("^(.*?)_Deletion", matched_df$gene_name, value = TRUE)
# Extract gene names (remove "_Deletion" part)
gene_names <- sub("^(.*)_Deletion$", "\\1", deletion_entries)
# Create a new dataframe amp_df with extracted gene names
del_df <- data.frame(gene_name = gene_names)
colnames(del_df) <- "Deletions"

# Saving: Create a new Excel workbook
wb <- createWorkbook()
# Add data frames to the workbook as separate sheets
addWorksheet(wb, "Sheet1")
writeData(wb, sheet = 1, x = matched_df, startCol = "A", startRow = 1)
writeData(wb, sheet = 1, x = matching_genes, startCol = "B", startRow = 1)
writeData(wb, sheet = 1, x = del_df, startCol = "E", startRow = 1)
writeData(wb, sheet = 1, x = amp_df, startCol = "D", startRow = 1)

# Save the workbook to an Excel file
saveWorkbook(wb, "matched_s4b_top60prox.xlsx", overwrite = TRUE)

# Distal
# Create an empty vector to store matched entries
matched_entries <- c()
# Loop through each gene in matching_genes$exact_matches
for (gene in matching_genes$exact_matches) {
  # Use grep to find exact matches (ignoring the description part after "_")
  matches <- grep(paste0("^", gene, "_"), mutations$s4bliver_unique, value = TRUE)
  # Append the matched entries to the vector
  matched_entries <- c(matched_entries, matches)
}
# Create a new dataframe with matched entries
matched_df <- as.data.frame(mutations[mutations$s4bliver_unique %in% matched_entries, ])
colnames(matched_df) <- "gene_name"

#now that we got all the entries in s4b that correspond to LRs of the targets we can divide them in amp and dels
# Use grep to filter entries ending with "_Amplification"
amplification_entries <- grep("^(.*?)_Amplification", matched_df$gene_name, value = TRUE)
# Extract gene names (remove "_Amplification" part)
gene_names <- sub("^(.*)_Amplification$", "\\1", amplification_entries)
# Create a new dataframe amp_df with extracted gene names
amp_df <- data.frame(gene_name = gene_names)
colnames(amp_df) <- "Amplifications"

# Use grep to filter entries ending with "_Deletion"
deletion_entries <- grep("^(.*?)_Deletion", matched_df$gene_name, value = TRUE)
# Extract gene names (remove "_Deletion" part)
gene_names <- sub("^(.*)_Deletion$", "\\1", deletion_entries)
# Create a new dataframe amp_df with extracted gene names
del_df <- data.frame(gene_name = gene_names)
colnames(del_df) <- "Deletions"

# Saving: Create a new Excel workbook
wb <- createWorkbook()
# Add data frames to the workbook as separate sheets
addWorksheet(wb, "Sheet1")
writeData(wb, sheet = 1, x = matched_df, startCol = "A", startRow = 1)
writeData(wb, sheet = 1, x = matching_genes, startCol = "B", startRow = 1)
writeData(wb, sheet = 1, x = del_df, startCol = "E", startRow = 1)
writeData(wb, sheet = 1, x = amp_df, startCol = "D", startRow = 1)

# Save the workbook to an Excel file
saveWorkbook(wb, "matched_s4b_top60dist.xlsx", overwrite = TRUE)
