#Code for the analysis of bulk RNA sequencing data of livers of mice bearing orthotopic CRC tumors or naive mice

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
counts <-counts[-grep("__",rownames(counts)),]
colnames(counts)

#subset specific samples 
liver <- counts[,9:20]
colnames(liver)
colnames(liver) <- c("liver_AKPS_1", "liver_AKPS_2", "liver_AKPS_3", "liver_AKPS_4", "liver_sham_1", "liver_sham_2")


######### DIFFERENTIAL GENE EXPRESSION ANALYSIS WITH edgeR ##########
#Create DGEList object
group <- factor(c("liver_AKPS","liver_AKPS", "liver_AKPS", "liver_AKPS", "liver_sham", "liver_sham")) 
y <- DGEList(counts=liver,group=group)
y$samples
barplot(y$samples$lib.size,names=colnames(y), las=2, main = "Barplot of library sizes")


# Filter reads by counts
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~0 + group)
y <- estimateDisp(y,design, robust=T)
y$common.dispersion
plotMDS(y, top= 5000)
logcpm <- cpm(y, log=TRUE)
View(logcpm)

design <- model.matrix(~0 +group)
colnames(design) <- levels(group)
design
fit <- glmQLFit(y, design, robust=T)
head(fit$coefficients)
plotQLDisp(fit)

con<-makeContrasts(liver_AKPS - liver_sham,levels=design)
qlf_liver<-glmQLFTest(fit,contrast=con)

#annotate with external gene names
qlf_liver$table$genes <- mapIds(org.Mm.eg.db, 
                               keys=gsub("\\.[0-9]{1,2}","",rownames(qlf_liver)),
                               keytype = "ENSEMBL", column="SYMBOL", 
                               multiVals = "first")
View(qlf_liver$table)
tr<-glmTreat(fit,contrast=con,lfc=log2(1.5))
summary(decideTests(tr))

# Extract corrected DEGs using False Discovery Rate (FDR) adjustment
corrected_DEGs <- na.omit(topTags(qlf_liver, n = Inf, adjust.method = "fdr", sort.by = "PValue", p.value = 1)$table)
View(corrected_DEGs)

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
    geom_text_repel(size=3, max.overlaps = Inf, segment.size = 0.2,)+
    geom_label_repel(label = ifelse(DEGs_liver$genelabels == TRUE, as.character(DEGs_liver$genes),""), 
                     size=4, segment.size = 0.2, max.overlaps= Inf, color="black", fill="white")+
    xlab("log2 fold change") + ylim(0,upperylim)+xlim(-xlim,xlim)+
    ylab("-log10 adjusted p-value") +  labs(title = title)+
    theme(legend.position = "none", text = element_text(size=12),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
  return(plot)
}

DEGs_liver<-corrected_DEGs

# Initialize a new column named "genelabels" in the DEGs_liver dataframe with empty strings
DEGs_liver$genelabels <- ""

# Populate the "genelabels" column with TRUE if the value in the "genes" column is "Plxnb2", otherwise FALSE
DEGs_liver$genelabels <- ifelse(DEGs_liver$genes == "Plxnb2", TRUE, FALSE)

#plot
DEGs_volcano(DEGs_liver, 0.05, 10, "AKPS vs sham liver", "grey",9, 16)# +ggsave("liver_AKPSvsSham.pdf", width =7,  height=5)

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

BP_enrichment <- preranked_BP(DEGs_liver)
View(BP_enrichment)

ggplot(BP_enrichment %>% filter(abs(NES)>1 & padj<0.05) %>% head(n= 100), aes(reorder(pathway, NES), NES)) +
  geom_point(aes(size=size, colour=pval)) +
  scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[1:11]) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=" ") + 
  theme_classic()+
  labs(y="")+   
  theme(axis.text.y=element_text(size=10)) + ggsave("liver_BP.pdf", width =10,  height=8)


