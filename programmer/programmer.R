setwd("/projectnb2/bf528/users/vangogh2022/project_3/programmer/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("EnhancedVolcano")

# Load Packages
library(tidyverse)
library(dplyr)
library(data.table)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)

#Assign sample informatino to variables
sample_info <- read_csv("/project/bf528/project_3/groups/group_1_rna_info.csv")
sample_info_control <- sample_info %>% 
  filter(mode_of_action == "Control")

# Identify sample treatment
treatment_sample = sample_info %>% 
  filter(!mode_of_action %in% "Control")

# Identify sample controls
control_sample = sample_info %>% 
  filter(mode_of_action %in% "Control")

# Read Count Tables from featureCount
SRR1178046 <- fread("filecounts/SRR1178046_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178036 <- fread("filecounts/SRR1178036_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178020 <- fread("filecounts/SRR1178020_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178002 <- fread("filecounts/SRR1178002_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177999 <- fread("filecounts/SRR1177999_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177997 <- fread("filecounts/SRR1177997_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177989 <- fread("filecounts/SRR1177989_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177988 <- fread("filecounts/SRR1177988_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177987 <- fread("filecounts/SRR1177987_Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)

# Combine geneid and count from each file 
SRR1178046 <- SRR1178046[,c(1,7)] %>% 
  dplyr::rename(SRR1178046 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1178046_Aligned.sortedByCoord.out.bam") 
SRR1178036 <- SRR1178036[,c(1,7)] %>% 
  dplyr::rename(SRR1178036 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1178036_Aligned.sortedByCoord.out.bam")
SRR1178020 <- SRR1178020[,c(1,7)] %>% 
  dplyr::rename(SRR1178020 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1178020_Aligned.sortedByCoord.out.bam")
SRR1178002 <- SRR1178002[,c(1,7)] %>% 
  dplyr::rename(SRR1178002 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1178002_Aligned.sortedByCoord.out.bam")
SRR1177999 <- SRR1177999[,c(1,7)] %>% 
  dplyr::rename(SRR1177999 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1177999_Aligned.sortedByCoord.out.bam")
SRR1177997 <- SRR1177997[,c(1,7)] %>% 
  dplyr::rename(SRR1177997 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1177997_Aligned.sortedByCoord.out.bam")
SRR1177989 <- SRR1177989[,c(1,7)] %>% 
  dplyr::rename(SRR1177989 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1177989_Aligned.sortedByCoord.out.bam")
SRR1177988 <- SRR1177988[,c(1,7)] %>% 
  dplyr::rename(SRR1177988 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1177988_Aligned.sortedByCoord.out.bam")
SRR1177987 <- SRR1177987[,c(1,7)] %>% 
  dplyr::rename(SRR1177987 = "/projectnb/bf528/users/vangogh2022/project_3/data_curator/STAR_results/SRR1177987_Aligned.sortedByCoord.out.bam")

# Merging the gids for each sample
id_m <- Reduce(function(...) merge(..., by = "Geneid", all=TRUE), list(SRR1178046, SRR1178036, SRR1178020, SRR1178002, SRR1177999, SRR1177997, SRR1177989, SRR1177988, SRR1177987))
rm(SRR1178046, SRR1178036, SRR1178020, SRR1178002, SRR1177999, SRR1177997, SRR1177989, SRR1177988, SRR1177987)

# Making distribution box plot 
boxplots <- id_m %>% pivot_longer(!Geneid, names_to = "Sample", values_to = "Count")
png(filename = "box_count_distribution.png")
boxplots %>% ggplot(mapping = aes(x = Sample, y = Count)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) + labs(title = "Boxplot of Count Distributions") + ylim(0, 1000)

dev.off()

### Deseq2 RNA Seq Analysis ###

# Load control counts 
control_counts <- read.csv("/project/bf528/project_3/samples/control_counts.csv") %>%
  select(Geneid, as.character(sample_info_control[[1]]))

controls_m <- merge(id_m, control_counts, by = "Geneid") %>% 
  column_to_rownames(var="Geneid")

filtered <- subset(controls_m,rowSums(controls_m==0)==0)

# Splitting by Mode of Action
ahr_meth <- sample_info %>% 
  filter(mode_of_action == "AhR" | mode_of_action == "Control" & vehicle == "CMC_.5_%") %>% 
  select(Run)
car_clot <- sample_info %>% 
  filter(mode_of_action == "CAR/PXR" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% 
  select(Run)
cyto_chloro <- sample_info %>% 
  filter(mode_of_action == "Cytotoxic" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% 
  select(Run)

# Making deseq objects

#ahr
obj_ahr <- DESeqDataSetFromMatrix(countData = filtered %>% select(ahr_meth[[1]]),colData = sample_info[c(1:3, 10:12),],design = ~ mode_of_action)

#CAR/PXR
obj_car <- DESeqDataSetFromMatrix(countData = filtered %>% select(car_clot[[1]]), colData = sample_info[c(4:6, 13:15),], design = ~ mode_of_action)

#Cyto 
obj_cyto <- DESeqDataSetFromMatrix(countData = filtered %>% select(cyto_chloro[[1]]),colData = sample_info[c(7:9, 13:15),],design = ~ mode_of_action)

# Factor for each DESeq object is reset with mode of action
obj_ahr$mode_of_action <- relevel(obj_ahr$mode_of_action, ref='Control')
obj_car$mode_of_action <- relevel(obj_car$mode_of_action, ref='Control')
obj_cyto$mode_of_action <- relevel(obj_cyto$mode_of_action, ref='Control')

# Execute DESeq
obj_ahr <- DESeq(obj_ahr)
obj_car <- DESeq(obj_car)
obj_cyto <- DESeq(obj_cyto)

# Results
res_ahr <- results(obj_ahr, contrast=c('mode_of_action','AhR','Control'))
res_ahr <- lfcShrink(obj_ahr, coef=2)

res_car <- results(obj_car, contrast=c('mode_of_action','CAR/PXR','Control'))
res_car <- lfcShrink(obj_car, coef=2)

res_cyto <- results(obj_cyto, contrast=c('mode_of_action','Cytotoxic','Control'))
res_cyto <- lfcShrink(obj_cyto, coef=2)

# Store counts
ahr_counts <- read.csv("deseq_ahr.csv", header = TRUE, row.names = 1)
car_counts <- read.csv("deseq_car.csv", header = TRUE, row.names = 1)
cyto_counts <- read.csv("deseq_cyto.csv", header = TRUE, row.names = 1)

# Number of genes significant at p-adj < 0.05
#ahr = 313 genes
ahr_counts %>% filter(padj < 0.05) %>% 
  summarize(count = n())

#car = 930 genes
car_counts %>% filter(padj < 0.05) %>% 
  summarize(count = n())

#cyto = 1728 genes
cyto_counts %>% filter(padj < 0.05) %>% 
  summarize(count = n())

# Creating data frame of significant gene counts
sig_genes <- data.frame(mode_of_action=c('AhR','CAR/PXR', 'Cytotoxic'), gene_count=c('313','930','1,728'))


# Top 10 diff expressed genes by p-value
ahr_counts %>% top_n(-10, padj)
car_counts %>% top_n(-10, padj)
cyto_counts %>% top_n(-10, padj)

write.csv("top10_ahr.csv", row.names = TRUE)
write.csv("top10_car.csv", row.names = TRUE)
write.csv("top10_cyto.csv", row.names = TRUE)

# Histograms
png("histo_ahr.png")
ahr_counts %>% filter(padj < 0.05) %>% ggplot(mapping = aes(x = log2FoldChange)) + geom_histogram() + labs(title = "AhR: Histogram Fold Change", y = "Frequency")
dev.off()

png("histo_car.png")
car_counts %>% filter(padj < 0.05) %>% ggplot(mapping = aes(x = log2FoldChange)) + geom_histogram() + labs(title = "CAR/PXR: Histogram Fold Change", x = "Log2 Fold Change", y = "Frequency")
dev.off()

png("histo_cyto.png")
cyto_counts %>% filter(padj < 0.05) %>% ggplot(mapping = aes(x = log2FoldChange)) + geom_histogram() + labs(title = "Cytotoxic: Histogram Fold Change", x = "Log2 Fold Change", y = "Frequency")
dev.off()


# Scatter Plots/Volcano Plots
png("scatter_ahr.png", width = 800, height = 480)
EnhancedVolcano(res_ahr, lab = NA, x = 'log2FoldChange', y = 'pvalue', title = 'AhR: Significant Differentially Expressed Genes', pCutoff = 0.05, FCcutoff = 1.5, legendPosition = 'right')
dev.off()

png("scatter_car.png", width = 800, height = 480)
EnhancedVolcano(res_car, lab = NA, x = 'log2FoldChange',y = 'pvalue', title = 'CAR/PXR: Significant Differentially Expressed Genes',pCutoff = 0.05,FCcutoff = 1.5,legendPosition = 'right')
dev.off()

png("scatter_cyto.png", width = 800, height = 480)
EnhancedVolcano(res_cyto,lab = NA,x = 'log2FoldChange',y = ' pvalue',title = 'Cytotoxic: Significant Differentially Expressed Genes',pCutoff = 0.05,FCcutoff = 1.5,legendPosition = 'right')
dev.off()
