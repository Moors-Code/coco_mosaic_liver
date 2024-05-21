#Code for the analysis of bulk RNA sequencing data of AKPS Sema4KD (quadruple sempahorin shRNA) vs emppty vector controls

#load packages
library(rlang)
library(devtools)
library(BiocManager)
library(biomaRt)
library("dplyr")
library(ggplot2)
library(RColorBrewer) 
library(ggraph)
library(wesanderson)
library(fgsea)
library(msigdbr)
library(ggrepel)
library(cowplot)
library(magrittr)
library(data.table)
library(stringr)
library(tidyr)
library("edgeR")

######### ANNOTATE AND CREATE COUNT MATRICES #########
#this function annotates the UMI count data with external gene names using Ensembl biomart for mouse species. The resulting dataframe contains gene expression data with annotated gene names.
Annotation_mouse <- function(zUMI_output, prefix) {
  ensembl<-useEnsembl(biomart="ensembl")
  list<-listDatasets(ensembl)
  mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version = 95)
  attributes<-listAttributes(mart)
  gene_ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  colnames(dataframe)<-paste(colnames(dataframe), prefix,sep = "_") 
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  join<-dataframe%>%
    left_join(dplyr::select(gene_ids,1:2))
  length(unique(join$external_gene_name))
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 #make all empty value to zero
  rownames(join)<-join$external_gene_name
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  return(join)
}


#read in dgecounts.rds file after zUMI demultiplexing, every pool contains multiple samples with individual barcodes introduced with the mcSCRB-seq protocol
#barcode to sample assignment
#barcode          sample               
# "ACAGTG_poolA"  mets in AMO25713          
#"AGTTCC_poolA"   AKPS sLPmCherry 3
#"ATCACG_poolA"   KPC shRNA EV 1
#"ATTCCT_poolA"   KPC parental
# "CACTCA_poolA"  AKPS shRNA Sema4s 2       
# "CTATAC_poolA"  KPC shRNA Sema4s 2
#"GCCAAT_poolA"   mets in AMO25715          
#"TACAGC_poolA"   KPC shRNA EV 2
#"TATAAT_poolA"   KPC shRNA EV 3
#"AGTCAA_poolB"   AKPS sLPmCherry 2
#"CCAACA_poolB"   AKPS shRNA EV 3           
#"CTAGCT_poolB"   KPC shRNA Sema4s 1
#"GTAGAG_poolB"   AKPS shRNA EV 1           
#"TAATCG_poolB"   AKPS shRNA Sema4s 1       
#"TCATTC_poolB"   mets in AMO25297          
#"TCGAAG_poolB"   AKPS shRNA Sema4s 3       
# "TGACCA_poolB"  mets in AMO25712          
# "TTAGGC_poolB"  mets in AMO25710          
# "ATCACG_poolC"  mets in AMO25294          
# "ATGAGC_poolC"  AKPS shRNA EV 2           
#"CAGATC_poolC"   mets in AMO25856          
#"CCGTCC_poolC"   mets derived org line 1
#"GACGAC_poolC"   mets derived org line 2
#"GAGTGG_poolC"   mets in AMO25934          
#"GTGAAA_poolC"   AKPS sLPmCherry 1

poolA <- readRDS("poolA/zUMIs_output/expressionbulk_orgs_round2.dgecounts.rds")
colnames(poolA$umicount$exon$all)
poolA_ann <- Annotation_mouse(poolA, "poolA")
colnames(poolA_ann) <- c("mets in AMO25713", "AKPS sLPmCherry 3","KPC parental", 
                         "AKPS shRNA Sema4s 2", "KPC shRNA EV 1", "KPC shRNA Sema4s 2", 
                         "mets in AMO25715", "KPC shRNA EV 2", "KPC shRNA EV 3")
poolA_ann$Gene <- rownames(poolA_ann)

poolB <- readRDS("poolB/zUMIs_output/expression/bulk_orgs_round2.dgecounts.rds")
colnames(poolB$umicount$exon$all)
poolB_ann <- Annotation_mouse(poolB, "poolB")
colnames(poolB_ann) <- c("AKPS sLPmCherry 2", 
                         "AKPS shRNA EV 3", 
                         "KPC shRNA Sema4s 1",
                         "AKPS shRNA EV 1", 
                         "AKPS shRNA Sema4s 1", 
                         "mets in AMO25297", 
                         "AKPS shRNA Sema4s 3", 
                         "mets in AMO25712", 
                         "mets in AMO25710")      
poolB_ann$Gene <- rownames(poolB_ann) 

poolC <- readRDS("poolC/zUMIs_output/expression/bulk_orgs_round2.dgecounts.rds")
colnames(poolC$umicount$exon$all)
poolC_ann <- Annotation_mouse(poolC, "poolC")
View(poolC_ann)
colnames(poolC_ann) <- c("mets in AMO25294" , "AKPS shRNA EV 2", "mets in AMO25856",
                                        "mets derived org line 1",  "mets derived org line 2",
                                        "mets in AMO25934" ,"AKPS sLPmCherry 1")
poolC_ann$Gene <- rownames(poolC_ann)


#create sample count matrix
poolAB <- merge(poolA_ann, poolB_ann,by="Gene",  stringsAsFactors = TRUE)
colnames(poolAB)
poolABC <- merge(poolAB, poolC_ann,by="Gene",  stringsAsFactors = TRUE)
colnames(poolABC)

all_samples_mat <- as.matrix(poolABC[2:26])
rownames(all_samples_mat)<- poolABC$Gene

########## DIFFERENTIAL GENE EXPRESSION TESTING WITH edgeR ###############
# Subset data to include specific columns containing data from AKPS Sema4KD and AKPS EV control
AKPS_Sema4_KD <- all_samples_mat[,c(5,11, 13,14, 16, 20)]
head(AKPS_Sema4_KD)

# Define group factors
group <- factor(c(2,1,1,2,2,1))

# Create a DGEList object
y <- DGEList(counts=AKPS_Sema4_KD,group=group)
y$samples

# Create a barplot of library sizes (the las argument rotates the axis names)
barplot(y$samples$lib.size,names=colnames(y), las=2, main = "Barplot of library sizes")

# Filter reads by counts and normalize the library
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# Create design matrix
design <- model.matrix(~group)

# Estimate dispersion
y <- estimateDisp(y,design)

# Plot MDS
plotMDS(y)

# Calculate log-transformed counts per million (logCPM)
logcpm <- cpm(y, log=TRUE)
View(logcpm)

#fit the model
design <- model.matrix(~group)
fit <- glmQLFit(y, design)

# Conduct differential expression analysis
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)

# Add gene names to the results table
qlf.2vs1$table$Gene <- rownames(qlf.2vs1$table)
View(qlf.2vs1$table)

# Export all genes
all_qlf.2vs1 = topTags(qlf.2vs1, n = Inf)
dim(qlf.2vs1)
head(qlf.2vs1$table)
qlf.2vs1$table$Gene <- rownames(qlf.2vs1$table)
View(qlf.2vs1$table)

# Extract corrected DEGs using False Discovery Rate (FDR) adjustment
corrected_DEGs <- na.omit(topTags(qlf_liver, n = Inf, adjust.method = "fdr", sort.by = "PValue", p.value = 1)$table)
View(corrected_DEGs)

#plot results

#this function generates a volcano plot to visualize differentially expressed genes (DEGs) based on specified thresholds for False Discovery Rate (FDR) and absolute log2 fold change. 
#It marks DEGs meeting the thresholds as 'Significant' and assigns gene names to them. The volcano plot is customized with colors, labels, axis limits, and title. The function returns the generated volcano plot.
DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$gene <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), rownames(DEGs), NA)
  
  plot <-  ggplot(data=DEGs, aes(x=logFC, y=-log10(FDR), fill=factor(Significant), label = DEGs$gene) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=4, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10(FDR)") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}

corrected_DEGs$gene <- rownames(corrected_DEGs)
DEGs_volcano(corrected_DEGs, 0.02, 3, "Sema-KD vs. control", "grey",3.5, 5) +ggsave("KD vs control AKPS.pdf", width =8,  height=8)


##############GENE SET ENRICHMENT ANALYSIS###############
#import dataset
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

BP_enrichment<- preranked_BP(qlf.2vs1$table)
View(BP_enrichment)

ggplot(BP_enrichment %>% filter(abs(NES)>2.5 & pval<0.005) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=10)) + ggsave("SemaKD_BP.pdf", width =15,  height=8)


