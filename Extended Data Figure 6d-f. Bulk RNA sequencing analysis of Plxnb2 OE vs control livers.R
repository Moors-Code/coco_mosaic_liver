#Code for the analysis of bulk RNA sequencing data of Plxnb2 OE vs control (sgNT) livers

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
library(org.Mm.eg.db)

#import bulk counts 
counts <- read.csv("bulk_counts.csv", row.names=1)
colnames(counts)
counts <-counts[-grep("__",rownames(counts)),]
colnames(counts)

#subset specific samples 

AAV <- counts[,1:8]
colnames(AAV)
colnames(AAV) <- c("liver_OE_1","liver_OE_2", "liver_OE_3", "liver_KO_1", "liver_NT_1", "liver_KO_2", "liver_KO_3", "liver_NT_2")


######### DIFFERENTIAL GENE EXPRESSION ANALYSIS WITH edgeR ##########
#Create DGEList object
group <- factor(c("liver_OE","liver_OE", "liver_OE", "liver_KO", "liver_NT", "liver_KO", "liver_KO", "liver_NT")) 
y <- DGEList(counts=AAV,group=group)
y$samples
barplot(y$samples$lib.size,names=colnames(y), las=2, main = "Barplot of library sizes")

# Filter reads by counts
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# Create a design matrix for the analysis
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y,design, robust=T)
y$common.dispersion
plotMDS(y, top= 500)

# Calculate log2 counts per million (logCPM) and view the results
logcpm <- cpm(y, log=TRUE)
View(logcpm)
write.table(logcpm, file="logcpm_OE_NT.txt")

# Fit a quasi-likelihood negative binomial generalized linear model (GLM) with robust estimation
fit <- glmQLFit(y, design, robust=T, adjust.method ="fdr")
head(fit$coefficients)
plotQLDisp(fit)

# Define the contrast for the differential expression analysis
con<-makeContrasts(liver_OE - liver_NT,levels=design)

# Perform quasi-likelihood F-test for the defined contrast
qlf_OE<-glmQLFTest(fit,contrast=con)

#annotate with external gene names
qlf_OE$table$genes <- mapIds(org.Mm.eg.db, 
                                keys=gsub("\\.[0-9]{1,2}","",rownames(qlf_OE)),
                                keytype = "ENSEMBL", column="SYMBOL", 
                                multiVals = "first")

View(qlf_OE$table)
tr<-glmTreat(fit,contrast=con,lfc=log2(1.5))
summary(decideTests(tr))

# Extract corrected DEGs using False Discovery Rate (FDR) adjustment
corrected_DEGs_OE <- topTags(qlf_OE, n = Inf, adjust.method = "fdr", sort.by = "PValue", p.value = 1)$table
View(corrected_DEGs_OE)


##generate volcano plot
#this function generates a volcano plot to visualize differentially expressed genes (DEGs) based on specified thresholds for False Discovery Rate (FDR) and absolute log2 fold change. 
#It marks DEGs meeting the thresholds as 'Significant' and assigns gene names to them. The volcano plot is customized with colors, labels, axis limits, and title. The function returns the generated volcano plot.
DEGs_volcano <- function(DEGs, p_treshold, FC_treshold, title, color, upperylim, xlim) {
  DEGs$Significant <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), "Significant", "Not Sig")
  DEGs$genes <- ifelse((DEGs$FDR < p_treshold & abs(DEGs$logFC) > FC_treshold ), DEGs$genes, NA)
  
  plot <-  ggplot(data=DEGs, aes(x=logFC, y=-log10(FDR), fill=factor(Significant), label = DEGs$genes) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point(shape = 21,color=color)+
    scale_fill_manual(values = c(color, "black")) +
    geom_text_repel(size=3, max.overlaps = Inf)+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10(FDR)") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}

DEGs_volcano(corrected_DEGs_OE, 0.5, 6, "Plxnb2-OE vs NT", "grey",1.7, 10)+ggsave("Plxnb2OEvsNT_volcano.pdf", width =7,  height=5)

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

BP_enrichment_OE <- preranked_BP(DEGs_OE)
View(BP_enrichment_OE)

ggplot(BP_enrichment_OE %>% filter(abs(NES)>1 & padj<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=padj)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=10)) + ggsave("PLxnb2OEvsNT_BP.pdf", width =7,  height=5)



########CELL TYPE ENRICHMENT WITH ENRICHR########
#save DEGs
write.csv(na.omit(corrected_DEGs_OE), file="DEGS_Plxnb2OEvsNT.csv") #this list of genes was used in EnrichR to estimate cell type composition

#import saved results
Enrichr_celltypedeconTM_downinOE <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/bulk_liver_and_lungs/Enrichr_celltypedeconTM_downinOE.txt")
head(Enrichr_celltypedeconTM_downinOE)
Enrichr_celltypedeconTM_upinOE <- read.delim("/media/Coco/MOSAIC LIVER/Experiments/bulk_liver_and_lungs/Enrichr_celltypedeconTM_upinOE.txt")
head(Enrichr_celltypedeconTM_upinOE)

# Select columns of interest from Enrichr_celltypedeconTM_downinOE dataframe
down <- Enrichr_celltypedeconTM_downinOE[, c(1,4,7)]

# Add a new column "group" with value "DOWN" to the down dataframe
down$group <- "DOWN"

# Invert the values of the "Odds.Ratio" column in the down dataframe
down$Odds.Ratio <- -down$Odds.Ratio

# Select columns of interest from Enrichr_celltypedeconTM_upinOE dataframe
up <- Enrichr_celltypedeconTM_upinOE[, c(1,4,7)]

# Add a new column "group" with value "UP" to the up dataframe
up$group <- "UP"

# Bind the down and up dataframes row-wise to create a new dataframe "plot"
plot <- rbind(down, up)

# Plotting the data using ggplot2 library
ggplot(data=plot, aes(x=plot$Odds.Ratio, y=-log10(plot$Adjusted.P.value), fill=factor(group)) ) + 
    theme_bw()+ 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank())+
    geom_point() +
    ylim(0,1.5)+ xlim(-4,4)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed",  color = "black", size=0.1)+
    geom_vline(xintercept = 0, linetype="dashed",  color = "black", size=0.1) +ggsave("Plxnb2OEvsNT_celltype.pdf", width =8,  height=6)
