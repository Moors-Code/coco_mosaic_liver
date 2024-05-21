## ANALYSIS PIPELINE

#demultiplexing with Bcl2Fastq, trimming with cutadapt, alignment with Bowtie2, and counting and testing with MAGeCK were performed in the terminal
#all the other analysis steps like plotting and gene set enrichment analysis on R

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#1. DEMULTIPLEXING FASTQS WITH Bcl2Fastq (Illumina), each sample has a unique barcode in the P7 primer. 

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#2. TRIMMING with cutadapt > cutadapt loop

cutadapt -g CACCG -o sampleX_trimup.fastq ../demux_fastqs/sampleX_merged.fastq.gz
cutadapt -a GTTTT -o sampleX_trimmed.fastq sampleX_trimup.fastq

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#3. ALIGNMENT with Bowtie2

awk -F ',' '{print ">"$1"\n"$2}' library2.csv > library2.fa #build bowtie index 
bowtie2-build library2.fa bowtie2_ind_library
bowtie2 -x bowtie2_ind_library -U ../cutadapt_output/sampleX_trimmed.fastq --norc | samtools view -bS - > sampleX.bam

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#4. LIBRARY RETENTION: calculate correlation between plasmid prep and post injection library (Extended Data Fig 1)

mageck count -l library2.csv -n premets_SPH1 --sample-label "premets,SPH1" --fastq 3071.bam SPH1.bam --norm-method total

premets_SPH1.count_normalized <- read.delim("premets_SPH1.count_normalized.txt")
View(premets_SPH1.count_normalized)

library("ggpubr")
ggscatter(premets_SPH1.count_normalized, x = "preinj", y = "plasmid", 
          add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
          fill = "#83C441",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "sgRNA count hepatocytes", ylab = "sgRNA count plasmid")+
          ggsave("plasmid_preinj.pdf", width = 5, height = 5)

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#5. LIBRARY STATS (Extended Data Fig 1)

#count all Cre samples to get stats in countsummary
mageck count -l library2.csv -n SPH_sum --sample-label "distal,proximal" --fastq 3174distal.bam,3068distal.bam,3039distal.bam,3379distal.bam,3425distal.bam,3434distal.bam,3436proximal.bam 3174proximal.bam,3068proximal.bam,3039proximal.bam,3379proximal.bam,3425proximal.bam,3434proximal.bam,3436distal.bam  
SPH_sum.count_normalized <- read.delim("SPH_sum.count_normalized.txt", header=T)
View(SPH_sum.count_normalized)
SPH_sum.count_normalized$ratio <- SPH_sum.count_normalized$proximal/SPH_sum.count_normalized$distal

SPH_sum.countsummary <- read.delim("SPH_sum.countsummary.txt")
View(SPH_sum.countsummary)

#percent mapped reads
a<- ggplot(data=SPH_sum.countsummary, aes(x=File, y=Percentage, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+ ylim(0,1)+
  scale_x_discrete(guide = guide_axis(angle = 45))

#zerocounts
b <-ggplot(data=SPH_sum.countsummary, aes(x=File, y=Zerocounts, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+ ylim(0,50)+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+
  scale_x_discrete(guide = guide_axis(angle = 45))

#givi index (sgRNA evenness)
c <- ggplot(data=SPH_sum.countsummary, aes(x=File, y=GiniIndex, fill=Label)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+ ylim(0,1)+
  theme_classic()+ scale_fill_manual(values=c('#81C341','#D12026'))+
  scale_x_discrete(guide = guide_axis(angle = 45))

ggarrange(a,b,c, ncol = 1, nrow = 3) +ggsave("stats.pdf", width = 5, height = 11)


#check correlation between library batches (3 independent plasmid preps, SPH1 = Hi1, SPH2 = Hi2, SPH3 = Hi3)
cd /NAS/Coco/MOSAIC\ LIVER/Experiments/Screen1/Fastqs/MAGECK/Bam_files
mageck count -l library.csv -n SPH1 --sample-label "SPH1" --fastq A3.bam --norm-method total
mageck count -l library.csv -n SPH2 --sample-label "SPH2" --fastq SPH2.bam --norm-method total
mageck count -l library.csv -n SPH3 --sample-label "SPH3" --fastq SPH3.bam --norm-method total

#merge all counts
`SPH1.count` <- read.delim("SPH1.count.txt")
View(SPH1.count)
`SPH2.count` <- read.delim("SPH2.count.txt")
View(SPH2.count)
`SPH3.count` <- read.delim("SPH3.count.txt")
View(SPH3.count)

SPH1_2_3 <- list(`SPH1.count` ,`SPH2.count`, `SPH3.count`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
colnames(SPH1_2_3)
head(SPH1_2_3)
SPH1_2_3_counts <- SPH1_2_3[, c(2,3,4,7,10)]
colnames(SPH1_2_3_counts) <- c( "sgRNA" ,"Gene", "SPH1", "SPH2", "SPH3")

View(SPH1_2_3_counts)

#correlation
library("ggpubr")
a <- ggscatter(SPH1_2_3_counts, x = "SPH1", y = "SPH2", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI1_counts", ylab = "HI2_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)

b <- ggscatter(SPH1_2_3_counts, x = "SPH1", y = "SPH3", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI1_counts", ylab = "HI3_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)

c <- ggscatter(SPH1_2_3_counts, x = "SPH2", y = "SPH3", 
               add = "reg.line", conf.int = TRUE, conf.int.level = 0.95, color = "black",
               fill = "#83C441",
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "HI2_counts", ylab = "HI3_counts")#+ggsave("plasmid_preinj.pdf", width = 5, height = 5)
ggarrange(a,b,c,ncol = 3, nrow = 1) +ggsave("plasmid_SPH123.pdf", width = 10, height = 3)

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#6. COVERAGE (Extended Data Fig 1)

Coverageplot <- read.csv("Coverageplot.csv")
View(Coverageplot)

ggplot(data=Coverageplot, aes(x=Sample, y=Total, fill=library)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  theme_classic() + scale_fill_manual(values=c('#81C341','#818641', "#2F8641"))+
  scale_x_discrete(guide = guide_axis(angle = 45))  +ggsave("coverage.pdf", width = 4, height = 3)
 
#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#7. ANALYSIS OF EACH MOUSE INDIVIDUALLY, use mageck to count and to perform robust rank aggregation test

#SPH1
mageck count -l library2.csv -n 3068 --sample-label "distal,proximal" --fastq 3068distal.bam 3068proximal.bam --norm-method total
mageck test -k 3068.count_normalized.txt -t proximal -c distal -n 3068
`3068.gene_summary` <- read.delim("3068.gene_summary.txt")
View(`3068.gene_summary`)

mageck count -l library2.csv -n 3174 --sample-label "distal,proximal" --fastq 3174distal.bam 3174proximal.bam --norm-method total
mageck test -k 3174.count_normalized.txt -t proximal -c distal -n 3174
`3174.gene_summary` <- read.delim("3174.gene_summary.txt")
View(`3174.gene_summary`)

#nocre
mageck count -l library2.csv -n 3070 --sample-label "distal,proximal" --fastq 3070distal.bam 3070proximal.bam --norm-method total
mageck test -k 3070.count_normalized.txt -t proximal -c distal -n 3070
`3070.gene_summary` <- read.delim("3070.gene_summary.txt")
View(`3070.gene_summary`)


#SPH2
mageck count -l library2.csv -n 3039 --sample-label "distal,proximal" --fastq 3039distal.bam 3039proximal.bam --norm-method total
mageck test -k 3039.count_normalized.txt -t proximal -c distal -n 3039
`3039.gene_summary` <- read.delim("3039.gene_summary.txt")
View(`3039.gene_summary`)

mageck count -l library2.csv -n 3379 --sample-label "distal,proximal" --fastq 3379distal.bam 3379proximal.bam --norm-method total
mageck test -k 3379.count_normalized.txt -t proximal -c distal -n 3379
`3379.gene_summary` <- read.delim("3379.gene_summary.txt")
View(`3379.gene_summary`)

mageck count -l library2.csv -n 3425 --sample-label "distal,proximal" --fastq 3425distal.bam 3425proximal.bam --norm-method total
mageck test -k 3378.count_normalized.txt -t proximal -c distal -n 3378
`3378.gene_summary` <- read.delim("3378.gene_summary.txt")
View(`3378.gene_summary`)

#nocre
mageck count -l library2.csv -n 3378 --sample-label "distal,proximal" --fastq 3378distal.bam 3378proximal.bam --norm-method total
mageck test -k 3381.count_normalized.txt -t proximal -c distal -n 3381
`3381.gene_summary` <- read.delim("3381.gene_summary.txt")
View(`3381.gene_summary`)

mageck count -l library2.csv -n 3381 --sample-label "distal,proximal" --fastq 3381distal.bam 3381proximal.bam --norm-method total
mageck test -k 3381.count_normalized.txt -t proximal -c distal -n 3381
`3381.gene_summary` <- read.delim("3381.gene_summary.txt")
View(`3381.gene_summary`)


#SPH3
mageck count -l library2.csv -n 3434 --sample-label "distal,proximal" --fastq 3434distal.bam 3434proximal.bam --norm-method total
mageck test -k 3434.count_normalized.txt -t proximal -c distal -n 3434
`3434.gene_summary` <- read.delim("3434.gene_summary.txt")
View(`3434.gene_summary`)

mageck count -l library2.csv -n 3436 --sample-label "distal,proximal" --fastq 3436proximal.bam 3436distal.bam --norm-method total #proximal and distal were switched during lib prep
mageck test -k 3436.count_normalized.txt -t proximal -c distal -n 3436
`3436.gene_summary` <- read.delim("3436.gene_summary.txt")
View(`3436.gene_summary`)

#nocre
mageck count -l library2.csv -n 3428 --sample-label "distal,proximal" --fastq 3428distal.bam 3428proximal.bam --norm-method total
mageck test -k 3428.count_normalized.txt -t proximal -c distal -n 3428
`3428.gene_summary` <- read.delim("3428.gene_summary.txt")
View(`3428.gene_summary`)

mageck count -l library2.csv -n 3431 --sample-label "distal,proximal" --fastq 3431distal.bam 3431proximal.bam --norm-method total
mageck test -k 3431.count_normalized.txt -t proximal -c distal -n 3431
`3431.gene_summary` <- read.delim("3431.gene_summary.txt")
View(`3431.gene_summary`)

#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#8. PAIRED ANALYSIS OF EACH LIBRARY BATCH (Extended Data Fig 1)

#SPH1
`3068.count_normalized` <- read.delim("3068.count_normalized.txt")
`3174.count_normalized` <- read.delim("3174.count_normalized.txt")

SPH1 <- list(`3174.count_normalized` ,`3068.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
head(SPH1)
SPH1_counts <- SPH1[, c(2,3,4,5,8,9)]
colnames(SPH1_counts) <- c( "sgRNA" ,"Gene", "3174d","3174p","3068d","3068p")
head(SPH1_counts)
write.table(SPH1_counts, "SPH1_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH1_counts.txt -t 3174p,3068p -c 3174d,3068d -n SPH1_paired --paired
SPH1_paired.gene_summary <- read.delim("SPH1_paired.gene_summary.txt")
View(SPH1_paired.gene_summary)

library(MAGECKFlute)
gdata1 = ReadRRA(SPH1_paired.gene_summary)
View(gdata1)
gdata1$LogFDR = -log10(gdata1$FDR)
ScatterView(gdata1, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 20, y_cut=0.2, x_cut = 0.05, max)
VolcanoView(gdata1, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1, max.overlaps=Inf)+ 
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH1_paired.pdf", width = 5, height = 4)


#SPH2
`3039.count_normalized` <- read.delim("3039.count_normalized.txt")
`3379.count_normalized` <- read.delim("3379.count_normalized.txt")
`3425.count_normalized` <- read.delim("3425.count_normalized.txt") 

SPH2 <- list(`3039.count_normalized` ,`3379.count_normalized`, `3425.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
head(SPH2)
SPH2_counts <- SPH2[, c(2,3,4,5,8,9,12,13)]
colnames(SPH2_counts) <- c( "sgRNA" ,"Gene", "3039d","3039p","3379d","3379p", "3425d", "3425p" )
head(SPH2_counts)
write.table(SPH2_counts, "SPH2_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH2_counts.txt -t 3039p,3379p,3425p -c 3039d,3379d,3425d -n SPH2_paired --paired
SPH2_paired.gene_summary <- read.delim("SPH2_paired.gene_summary.txt")
View(SPH2_paired.gene_summary)

gdata2 = ReadRRA(SPH2_paired.gene_summary)
View(gdata2)
gdata2$LogFDR = -log10(gdata2$FDR)
ScatterView(gdata2, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 20, y_cut=0.2, x_cut = 0.05, max)
VolcanoView(gdata2, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1, max.overlaps=Inf)+ 
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH2_paired.pdf", width = 5, height = 4)


#SPH3
`3434.count_normalized` <- read.delim("3434.count_normalized.txt")
`3436.count_normalized` <- read.delim("3436.count_normalized.txt")

SPH3 <- list(`3434.count_normalized`,`3436.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  # NB: below we are IMPUTING the genes where they miss in some samples, with value 0:
  # many ways possible, I choose one simple enough but also ok for big data sets:
  # https://stackoverflow.com/questions/8161836/how-do-i-replace-na-values-with-zeros-in-an-r-dataframe
  mutate_all(~replace(., is.na(.), 0))
head(SPH3)
SPH3_counts <- SPH3[, c(2,3,4,5,8,9)]
colnames(SPH3_counts) <- c( "sgRNA" ,"Gene", "3434d","3434p","3436d","3436p")
head(SPH3_counts)
write.table(SPH3_counts, "SPH3_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH3_counts.txt -t 3434p,3436p -c 3434d,3436d -n SPH3_paired --paired
SPH3_paired.gene_summary <- read.delim("SPH3_paired.gene_summary.txt")
head(SPH3_paired.gene_summary)
View(SPH3_paired.gene_summary)

gdata3 = ReadRRA(SPH3_paired.gene_summary)
View(gdata3)
gdata3$LogFDR = -log10(gdata3$FDR)
ScatterView(gdata3, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata3, x = "Score", y = "FDR", Label = "id", x_cut = 0.05, y_cut = 0.965, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026')) +#+ ylim(0,3)+ xlim(-2,2) #+
  ggsave("/media/Coco/MOSAIC LIVER/Manuscript/Graphs/VolcanoSPH3_paired.pdf", width = 5, height = 4)


#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#8. PAIRED ANALYSIS OF ALL MICE AND BATCHES (Fig 1)

SPH1_2_3 <- list(`3174.count_normalized` ,`3068.count_normalized`, `3039.count_normalized`, 
                 `3379.count_normalized`, `3425.count_normalized`, `3434.count_normalized`, 
                 `3436.count_normalized`) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
colnames(SPH1_2_3)
SPH1_2_3_counts <- SPH1_2_3[, c(2,3, 4,5,8,9,12,13,16,17,20,21,24,25,28,29)]
colnames(SPH1_2_3_counts) <- c( "sgRNA" ,"Gene", "3174d","3174p","3068d","3068p","3039d","3039p",
                                "3379d", "3379p","3425d", "3425p", "3434d", "3434p", "3436d", "3436p")

head(SPH1_2_3_counts)
write.table(SPH1_2_3_counts, "SPH1_2_3_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)
mageck test -k SPH1_2_3_counts.txt -t 3174p,3068p,3039p,3379p,3425p,3434p,3436p -c 3174d,3068d,3039d,3379d,3425d,3434d,3436d -n paired_proximal_distal --paired


paired_proximal_distal.sgrna_summary <- read.delim("paired_proximal_distal.sgrna_summary.txt")
sdata = ReadsgRRA(paired_proximal_distal.sgrna_summary)
View(sdata)
sgRankView(sdata)
sgRankView(sdata, gene = c("Psen1","Plxnb2","App", "Saa1"))+ggsave("sgRNAplotpaired.pdf", width = 5, height = 4)

paired_proximal_distal.gene_summary <- read.delim("paired_proximal_distal.gene_summary.txt")
View(paired_proximal_distal.gene_summary)
gdata = ReadRRA(paired_proximal_distal.gene_summary)
View(gdata)
gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 10, x_cut = 0.05, y_cut = 0.8, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-1,1) +
 ggsave("VolcanoSPH123_paired_all.pdf", width = 5, height = 4)

ggplot(data=paired_proximal_distal.gene_summary, aes(x=neg.lfc, y=-log10(neg.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()
ggplot(data=paired_proximal_distal.gene_summary, aes(x=pos.lfc, y=-log10(pos.p.value), label=id)) + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") + geom_text()+
  geom_hline(yintercept=-log10(0.05), col="red")+
  geom_point() + 
  theme_classic()

#barplot by score 
neg<-which(paired_proximal_distal.gene_summary$neg.lfc<0)
pos<-which(paired_proximal_distal.gene_summary$pos.lfc>0)
colnames(paired_proximal_distal.gene_summary)
neg_data <- paired_proximal_distal.gene_summary[neg, c(1,8)]
head(neg_data)
pos_data <- paired_proximal_distal.gene_summary[pos, c(1,14)]
head(pos_data)
colnames(pos_data) <-  c("id" ,"Score"   )
colnames(neg_data) <-  c("id" ,"Score"   )
plot_data <-rbind(neg_data, pos_data)
#View(plot_data)
plot_data$group <- ifelse(plot_data$lfc < 0, "neg", "pos")
reduced_plot_data <- plot_data %>%
  group_by(group) %>%
  top_n(n = 10, wt = abs(Score))
reduced_plot_data

ggplot(data=reduced_plot_data, aes(x=reorder(id, Score),y= score,  fill=group)) +
  geom_bar(stat="identity")+
  theme_classic() + scale_fill_manual(values=c('#81C341','#D12026'))+ coord_flip()+
  ggsave("SPH123_allpaired_score.pdf", width = 5, height = 4)


#gene set enrichment analysis
library(fgsea)
library(msigdbr)
library(data.table)
library(stringr)
library(dplyr)

#select species and set
msigdbr_show_species()
m_df<- msigdbr(species = "Mus musculus", category = "H")
Hallmarks <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(Hallmarks)
m_df<- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
REACTOME <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(REACTOME)
m_df<- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
GOBP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head(GOBP)

head(gdata)
ranks <- gdata %>% 
  na.omit()%>%
  mutate(ranking=-log10(FDR)/sign(Score))
ranks <- ranks$ranking
names(ranks) <- gdata$id
head(ranks, 10)

BP<- fgsea(pathways = GOBP, 
               stats = ranks,
               minSize=10,
               maxSize=500,
               nperm=1000000)
View(BP)

BP$pathway<-gsub("GOBP_","",BP$pathway)
BP$pathway<-gsub("_"," ",BP$pathway)
BP$pathway <- tolower(BP$pathway)
head(BP)

ggplot(BP %>% filter(abs(NES)>1) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  #geom_bar(stat="identity")+
  geom_point(aes(size=size, colour=NES)) +
  scale_size_area(max_size = 8)+
  scale_colour_gradient( low = '#81C341', high = '#D12026' aesthetics = "colour")+ 
  coord_flip() +
  theme_classic()  +
  ggsave("SPH123_allpaired_GOBP.pdf", width = 10, height = 4)


#---------------------------––-----------------------------––-----------------------------––-----------------------------––-------------------------#

#9. ANALYSIS OF ALBCRE:dCAS9-SPH VS. NOCRE LITTERMATES (Extended Data Fig 1) 
#As number of mice is not the same, cannot compare with paired analysis > first sum all mice per batches, then paired analysis 

#sum all no Cre mice for SPH1, 2 and 3
mageck count -l library2.csv -n SPH1nocre_sum --sample-label "distal,proximal" --fastq 3070distal.bam 3070proximal.bam  
mageck count -l library2.csv -n SPH2nocre_sum --sample-label "distal,proximal" --fastq 3378distal.bam,3381distal.bam 3378proximal.bam,3381proximal.bam 
mageck count -l library2.csv -n SPH3nocre_sum --sample-label "distal,proximal" --fastq 3428distal.bam,3431proximal.bam 3428proximal.bam,3431distal.bam  

`SPH1nocre_sum.count_normalized` <- read.delim("SPH1nocre_sum.count_normalized.txt")
`SPH2nocre_sum.count_normalized` <- read.delim("SPH2nocre_sum.count_normalized.txt")
`SPH3nocre_sum.count_normalized` <- read.delim("SPH3nocre_sum.count_normalized.txt") 

#paired analysis
SPHnocre <- list(`SPH1nocre_sum.count_normalized` ,`SPH2nocre_sum.count_normalized`,`SPH3nocre_sum.count_normalized` ) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
head(SPHnocre)
SPHnocre_counts <- SPHnocre[, c(2,3,4,5,8,9,12,13)]
head(SPHnocre_counts)
colnames(SPHnocre_counts) <- c( "sgRNA" ,"Gene", "SPH1d","SPH1p","SPH2d","SPH2p", "SPH3d","SPH3p")
View(SPHnocre_counts)
write.table(SPHnocre_counts, "SPHnocre_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPHnocre_counts.txt -t SPH1p,SPH2p,SPH3p -c SPH1d,SPH2d,SPH3d -n SPH123nocre_paired --paired
SPH123nocre_paired.gene_summary <- read.delim("SPH123nocre_paired.gene_summary.txt")
View(SPH123nocre_paired.gene_summary)


gdata = ReadRRA(SPH123nocre_paired.gene_summary)
View(gdata)
gdata<- gdata[-which(gdata$id =="ctrl"),]
gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
            model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 2, x_cut = 0.05, y_cut = 0.8, alpha=1)+
  theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-2,2) +
  ggsave("Volcanonocre.pdf", width = 5, height = 4)


#cre mice
mageck count -l library2.csv -n SPH1_sum --sample-label "distal,proximal" --fastq 3174distal.bam,3068distal.bam 3174proximal.bam,3068proximal.bam  
mageck count -l library2.csv -n SPH2_sum --sample-label "distal,proximal" --fastq 3039distal.bam,3379distal.bam,3425distal.bam 3039proximal.bam,3379proximal.bam,3425proximal.bam  
mageck count -l library2.csv -n SPH3_sum --sample-label "distal,proximal" --fastq 3434distal.bam,3436proximal.bam 3434proximal.bam,3436distal.bam  #3436 was switched in lib prep so here reverse

`SPH1_sum.count_normalized` <- read.delim("SPH1_sum.count_normalized.txt")
`SPH2_sum.count_normalized` <- read.delim("SPH2_sum.count_normalized.txt")
`SPH3_sum.count_normalized` <- read.delim("SPH3_sum.count_normalized.txt") 

#paired analysis
SPH <- list(`SPH1_sum.count_normalized` ,`SPH2_sum.count_normalized`,`SPH3_sum.count_normalized` ) %>% 
  lapply(tibble::rownames_to_column) %>% purrr::reduce(full_join, by="sgRNA") %>% 
  mutate_all(~replace(., is.na(.), 0))
head(SPH)
SPH_counts <- SPH[, c(2,3,4,5,8,9,12,13)]
head(SPH_counts)
colnames(SPH_counts) <- c( "sgRNA" ,"Gene", "SPH1d","SPH1p","SPH2d","SPH2p", "SPH3d","SPH3p")
head(SPH_counts)
write.table(SPH_counts, "SPH_counts.txt",append = FALSE, sep= "\t",  row.names = FALSE, quote = FALSE, col.names = TRUE)

mageck test -k SPH_counts.txt -t SPH1p,SPH2p,SPH3p -c SPH1d,SPH2d,SPH3d -n SPH123_paired --paired
SPH123_paired.gene_summary <- read.delim("SPH123_paired.gene_summary.txt")
View(SPH123_paired.gene_summary)

gdata = ReadRRA(SPH123_paired.gene_summary)
gdata$LogFDR = -log10(gdata$FDR)
ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 10, y_cut=0.25, x_cut = 0.05)
VolcanoView(gdata, x = "Score", y = "FDR", Label = "id", top = 2, x_cut = 0.05, y_cut = 0.8, alpha=1)+
theme_classic() + scale_fill_manual(values=c('#81C341', "grey", '#D12026'))+ ylim(0,3)+ xlim(-2,2) +
ggsave("Volcano.pdf", width = 5, height = 4)

