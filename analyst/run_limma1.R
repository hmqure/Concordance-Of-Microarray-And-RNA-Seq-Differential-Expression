library(limma)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_1_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# subset the full expression matrix to just those in this comparison
treatment1 <- samples[samples$chemical=="3-METHYLCHOLANTHRENE",]
control1 <- samples[samples$vehicle=="CMC_.5_%" & samples$chemical=="Control",]
group1<- rbind(treatment1, control1)

treatment2 <- samples[samples$chemical=="CLOTRIMAZOLE",]
control2 <- samples[samples$vehicle=="CORN_OIL_100_%" & samples$chemical=="Control",]
group2<- rbind(treatment2, control2)

treatment3 <- samples[samples$chemical=="CHLOROFORM",]
group3 <- rbind(treatment3, control1) # 3-methyl and chloro both use same vehicle

rmasubset1 <- rma[paste0('X',group1$array_id)]
rmasubset2 <- rma[paste0('X',group2$array_id)]
rmasubset3 <- rma[paste0('X',group3$array_id)]

# construct a design matrix modeling treatment vs control for use by limma
design1 <- model.matrix(~factor(group1$chemical, levels = c("Control", "3-METHYLCHOLANTHRENE")))
colnames(design1) <- c('Intercept','3-METHYLCHOLANTHRENE')
design2 <- model.matrix(~factor(group2$chemical, levels = c("Control", "CLOTRIMAZOLE")))
colnames(design2) <- c('Intercept','CLOTRIMAZOLE')
design3 <- model.matrix(~factor(group3$chemical, levels = c("Control", "CHLOROFORM")))
colnames(design3) <- c('Intercept', 'CHRLORFORM')

# run limma
fit1 <- lmFit(rmasubset1, design1)
fit1 <- eBayes(fit1)
t1 <- topTable(fit1, coef=2, n=nrow(rmasubset1), adjust='BH')
t1 <- dplyr::arrange(t1, adj.P.Val)

fit2 <- lmFit(rmasubset2, design2)
fit2 <- eBayes(fit2)
t2 <- topTable(fit2, coef=2, n=nrow(rmasubset2), adjust='BH')
t2 <- dplyr::arrange(t2, adj.P.Val)

fit3 <- lmFit(rmasubset3, design3)
fit3 <- eBayes(fit3)
t3 <- topTable(fit3, coef=2, n=nrow(rmasubset3), adjust='BH')
t3 <- dplyr::arrange(t3, adj.P.Val)

# write out the results to file
write.csv(t1,'3-methylcholanthrene_limma_results.csv')
write.csv(t2, 'clotrimazole_limma_results.csv')
write.csv(t3, 'chloroform_limma_results.csv')

#plotting

t1$volcano <- "NS"
t1$volcano[t1$logFC > 1.5 & t1$P.Val < 0.05] <- "UP"
t1$volcano[t1$logFC <= -1.5 & t1$P.Val < 0.05] <- "DOWN"
t1$volcano <- factor(t1$volcano)

t2$volcano <- "NS"
t2$volcano[t2$logFC > 1.5 & t2$P.Val < 0.05] <- "UP"
t2$volcano[t2$logFC <= -1.5 & t2$P.Val < 0.05] <- "DOWN"
t2$volcano <- factor(t2$volcano)

t3$volcano <- "NS"
t3$volcano[t3$logFC > 1.5 & t3$P.Val < 0.05] <- "UP"
t3$volcano[t3$logFC <= -1.5 & t3$P.Val < 0.05] <- "DOWN"
t3$volcano <- factor(t3$volcano)


sc1 <- ggplot(data = t1, aes(x=logFC, y=-log(P.Value, base=10), color = volcano)) + geom_point() + scale_color_manual(name = 'Expression Status', values = c('coral', 'gray60', 'cornflowerblue')) + labs(title="3-Methylcholanthrene", x='log2FoldChange', y ='-log10(pvalue)') + theme_classic() + theme(plot.title = element_text(face = "bold", size=10)) 
sc2 <- ggplot(data = t2, aes(x=logFC, y=-log(P.Value, base=10), color = volcano)) + geom_point() + scale_color_manual(name = 'Expression Status', values = c('coral', 'gray60', 'cornflowerblue')) + labs(title="Clotrimazole", x='log2FoldChange', y ='-log10(pvalue)') + theme_classic() + theme(plot.title = element_text(face = "bold", size=10))
sc3 <- ggplot(data = t3, aes(x=logFC, y=-log(P.Value, base=10), color = volcano)) + geom_point() + scale_color_manual(name = 'Expression Status', values = c('coral', 'gray60', 'cornflowerblue')) + labs(title="Chloroform", x='log2FoldChange', y ='-log10(pvalue)') + theme_classic() + theme(plot.title = element_text(face = "bold", size=10))
figure2 <- ggarrange(sc1, sc2, sc3, nrow=1, top="Scatter Plots of log2FoldChange vs Pvalue")
ggsave('scatter_plot.png', plot=figure2, device='png', width=15, bg='transparent')

t1 <- dplyr::filter(t1, adj.P.Val < 0.05)
t2 <- dplyr::filter(t2, adj.P.Val < 0.05)
t3 <- dplyr::filter(t3, adj.P.Val < 0.05)

h1 <- ggplot(data = t1, aes(x=logFC)) + geom_histogram(color="seagreen", fill=NA, bins = 20) + labs(title="3-Methylcholanthrene", x='log2FoldChange', y ='Counts') + theme_classic() + theme(plot.title = element_text(face='bold', size=10))
h2 <- ggplot(data = t2, aes(x=logFC)) + geom_histogram(color="seagreen", fill=NA, bins = 20) + labs(title="Clotrimazole", x='log2FoldChange', y ='Counts') + theme_classic() + theme(plot.title = element_text(face='bold', size=10))
h3 <- ggplot(data = t3, aes(x=logFC)) + geom_histogram(color="seagreen", fill=NA, bins = 20) + labs(title="Chloroform", x='log2FoldChange', y ='Counts') + theme_classic() + theme(plot.title = element_text(face='bold', size=10))
figure1 <- ggarrange(h1, h2, h3, nrow=1, top = "Histograms of log2FoldChange") 
ggsave('histogram.png', plot=figure1, device='png', bg='transparent', width = 15)


figure3 <- ggarrange(sc1, h1, sc2, h2, sc3, h3, labels = 'auto', nrow=3, ncol=2)
ggsave('part5_plots.png', plot = figure3, device = 'png', bg = 'transparent', height = 10, width = 10)



