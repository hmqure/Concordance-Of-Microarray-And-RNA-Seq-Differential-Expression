# libraries
library(tidyverse)
library(zeallot)
library(gridExtra)
library(ggpubr)
library(cowplot)

# read files
map <- read_csv('refseq_affy_map.csv', col_names = TRUE)

limma1 <- read_csv('3-methylcholanthrene_limma_results.csv', col_names = TRUE) %>%
  dplyr::rename('PROBEID'='...1') %>%
  dplyr::filter(adj.P.Val < 0.05)

deseq1 <- read_csv('deseq_ahr.csv', col_names = TRUE) %>%
  dplyr::rename('REFSEQ'='...1') %>%
  dplyr::filter(padj < 0.05)

limma2 <- read_csv('chloroform_limma_results.csv', col_names = TRUE) %>%
  dplyr::rename('PROBEID'='...1') %>%
  dplyr::filter(adj.P.Val < 0.05)

deseq2 <- read_csv('deseq_cyto.csv', col_names = TRUE) %>%
  dplyr::rename('REFSEQ'='...1') %>%
  dplyr::filter(padj < 0.05)

limma3 <- read_csv('clotrimazole_limma_results.csv', col_names = TRUE) %>%
  dplyr::rename('PROBEID'='...1') %>%
  dplyr::filter(adj.P.Val < 0.05)

deseq3 <- read_csv('deseq_car.csv', col_names = TRUE) %>%
  dplyr::rename('REFSEQ'='...1') %>%
  dplyr::filter(padj < 0.05)



# write a function that converts affy probes to RefSeq ids or gene symbols

convert <- function(df, map, type) {
  if (type == 'limma'){
    df <- left_join(df, map, by = 'PROBEID') %>%
      dplyr::select(SYMBOL, PROBEID, REFSEQ, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
      distinct(PROBEID, .keep_all = TRUE) %>%
      distinct(SYMBOL, .keep_all = TRUE)
  }
  if (type == 'deseq') {
    df <- left_join(df, map, by = 'REFSEQ') %>%
      dplyr::select(SYMBOL, PROBEID, REFSEQ, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
      distinct(REFSEQ, .keep_all = TRUE) %>%
      distinct(SYMBOL, .keep_all = TRUE)
  }
  return(df)
}

# write a function to calculate concordance score of two deseq results

concordance <- function(limma, deseq) {
  
  N <- 25000 #estimate of number of genes in rat genome
  n1 <- length(unique(limma$SYMBOL))
  n2 <- length(unique(deseq$SYMBOL))
  merged <- dplyr::inner_join(limma, deseq, by = 'SYMBOL')
  merged$rnaFC <- if_else(merged$log2FoldChange > 0, 'Pos', 'Neg')
  merged$micFC <- if_else(merged$logFC > 0, 'Pos', 'Neg')
  merged <- merged[merged$rnaFC == merged$micFC,]
  n0 <- length(unique(merged$SYMBOL))
  backgound <- (N*n0 - n1*n2)/(n0+N-n1-n2)
  concord <- concordance <- (2*abs(backgound))/(n1+n2)
  
  return(concord)
}

# write a function to separate above-median and below-median genes

sep_median_genes <- function(df, type) {
  if (type == 'limma'){
    above <- dplyr::filter(df, AveExpr > median(AveExpr))
    below <- dplyr::filter(df, AveExpr <= median(AveExpr))
  }
  if (type == 'deseq') {
    above <- dplyr::filter(df, baseMean > median(baseMean))
    below <- dplyr::filter(df, baseMean <= median(baseMean))
  }
  return(list(above, below))
}

#convert all dataframes to include gene symbols/affy ids/redseq ids

limma1 <- convert(limma1, map, 'limma')
top10_limma1 <- dplyr::select(limma1, SYMBOL, PROBEID, logFC, P.Value, adj.P.Val) %>% rename('Gene Symbol'=SYMBOL) 
limma2 <- convert(limma2, map, 'limma')
top10_limma2 <- dplyr::select(limma2, SYMBOL, PROBEID, logFC, P.Value, adj.P.Val) %>% rename('Gene Symbol'=SYMBOL) 
limma3 <- convert(limma3, map, 'limma')
top10_limma3 <- dplyr::select(limma3, SYMBOL, PROBEID, logFC, P.Value, adj.P.Val) %>% rename('Gene Symbol'=SYMBOL) 
deseq1 <- convert(deseq1, map, 'deseq')
deseq2 <- convert(deseq2, map, 'deseq')
deseq3 <- convert(deseq3, map, 'deseq')

#write top 10s to csv
write_csv(head(top10_limma1, n=10), 'top10_limma1.csv')
write_csv(head(top10_limma2, n=10), 'top10_limma2.csv')
write_csv(head(top10_limma3, n=10), 'top10_limma3.csv')


#calculate overall concordance
concord1 <- concordance(limma1, deseq1)
concord2 <- concordance(limma2, deseq2)
concord3 <- concordance(limma3, deseq3)

#subset by above and below median genes
c(limma1_above, limma1_below) %<-% sep_median_genes(limma1, 'limma')
c(limma2_above, limma2_below) %<-% sep_median_genes(limma2, 'limma')
c(limma3_above, limma3_below) %<-% sep_median_genes(limma3, 'limma')
c(deseq1_above, deseq1_below) %<-% sep_median_genes(deseq1, 'deseq')
c(deseq2_above, deseq2_below) %<-% sep_median_genes(deseq2, 'deseq')
c(deseq3_above, deseq3_below) %<-% sep_median_genes(deseq3, 'deseq')

#calculate concordance for above and below median genes
concord1_above <- concordance(limma1_above, deseq1_above)
concord2_above <- concordance(limma2_above, deseq2_above)
concord3_above <- concordance(limma3_above, deseq3_above)
concord1_below <- concordance(limma1_below, deseq1_below)
concord2_below <- concordance(limma2_below, deseq2_below)
concord3_below <- concordance(limma3_below, deseq3_below)



limma_deg_counts <- c(length(unique(limma1$SYMBOL)), length(unique(limma2$SYMBOL)), length(unique(limma3$SYMBOL)))
deseq_deg_counts <- c(length(unique(deseq1$SYMBOL)), length(unique(deseq2$SYMBOL)), length(unique(deseq3$SYMBOL))) #c(313, 1728, 930) #from programmer analysis
concordance_overall <- c(concord1, concord2, concord3)
Treatment <- c('3ME', 'CHR', 'CLO')
plotting1<- data.frame(limma_deg_counts, deseq_deg_counts, concordance_overall, as.factor(Treatment))

sc1 <- ggplot(plotting1, aes(x=limma_deg_counts, y=concordance_overall, color=Treatment)) + geom_point(shape=16) + labs(x='Number of DEGs from Microarray', y='Concordance Score')
sc2 <- ggplot(plotting1, aes(x=deseq_deg_counts, y=concordance_overall, color=Treatment)) + geom_point(shape=16) + labs(x='Number of DEGs from RNA-Seq', y='Concordance Score')

concordance_oab <- c(concord1, concord2, concord3, concord1_above, concord2_above, concord3_above, concord1_below, concord2_below, concord3_below)
Status <- rep(c('overall', 'above', 'below'), each=3)
moa <- rep(c('AhR', 'Cyto', 'CAR/PXR'), times=3)
plotting2 <- data.frame(concordance_oab, Status)
plotting2$Status <- as.factor(plotting2$Status)
plotting2$moa <- as.factor(plotting2$moa)

bar <- ggplot(plotting2, aes(fill=Status, y=concordance_oab, x=moa)) + geom_bar(position="dodge", stat="identity") + labs(x='MOA', y='Concordance Score')

sc_figure <- ggarrange(ggarrange(sc1,NULL, sc2, ncol=3, widths = c(1, 0.02, 1), labels = c("a", "", "b"), common.legend = TRUE, legend = 'right'), NULL, bar, nrow=3, heights = c(1, 0.1, 1), labels=c('', '', 'c'))
ggsave('concordance_sc.png', plot=sc_figure, device='png', bg='transparent', height=6, width=8)



