library(tidyverse)
library(DESeq2)
library(PCAtools)
library(RColorBrewer)
library(enrichR)
library(openxlsx)
library(EnhancedVolcano)
library(paletteer)

# Load Biomart
load("2022-09-30 Biomart Load.RData")

# Load counts and meta data and tidy up
counts <- read.csv("ALS_vs_Control_counts.csv", row.names=1)
samples <- read.csv("ALS_vs_Control_metadata.csv", row.names=1)
samples[sapply(samples, is.character)] <- lapply(samples[sapply(samples, is.character)], 
                                                 as.factor)
rownames(samples) <- make.names(rownames(samples))
colnames(counts) <- make.names(colnames(counts))
identical(rownames(samples),colnames(counts))


# Subset Samples 
table(samples$Disease)


# Overlap metadata and counts
counts <- counts[,colnames(counts) %in% rownames(samples)]
table(rownames(samples) == colnames(counts))
identical(rownames(samples),colnames(counts))


# Relevel samples add seq depth to meta
samples$Disease <- relevel(samples$Disease, "Control")
counts <- round(counts, 0)
plot(colSums(counts) / 1e6, main = "Counts Per Million")
samples$seqdepth <- colSums(counts)/ 1e6
samples$Batch <- as.factor(samples$Batch)
samples$group <- paste(samples$Cell_Line, samples$Batch, sep = "_")
# Generate DESeq object
dds <- DESeqDataSetFromMatrix(counts, colData=samples, design= ~Batch+Disease)

# Gene filtering
# Filter out genes where there are less than 5 samples with normalized counts greater than or equal to 10
dds <- estimateSizeFactors(dds)
keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 5
table(keep)
dds <- dds[keep,]  

# Run DEseq
dds <- DESeq(dds) 

# Run diff expression
resultsNames(dds)
res1 <- results(dds, name = "Disease_ALS_vs_Control", pAdjustMethod = "bonferroni")
final1 <- as.data.frame(res1)
final1 <- final1[order(final1$padj),] 
table(final1$padj<0.05)

# Add gene names
geneids <- as.data.frame(rownames(final1))
colnames(geneids) [1] <- "orig.ident"
geneids$ensembl <- gsub("\\..*","",geneids$orig.ident)
geneids2 <- left_join(geneids, gene, by = c("ensembl" = "ensembl_gene_id"))
geneids2 <- geneids2[!duplicated(geneids2$ensembl),]
geneids3 <- left_join(geneids, geneids2, by = "ensembl")
identical(geneids3[['ensembl']],geneids[['ensembl']])
geneids3$hgnc_symbol[geneids3$hgnc_symbol==""] <- NA
geneids3 <- geneids3 %>% mutate(symbol = coalesce(hgnc_symbol,ensembl))
identical(rownames(final1), geneids3$orig.ident.x)
final1$gene <- geneids3$symbol



# Export variance stablized data
vsd <- vst(dds, blind=TRUE)

# Before batch correction
quantLog <- assay(vsd)
table(rownames(samples) == colnames(quantLog))
geneids <- as.data.frame(rownames(quantLog))
colnames(geneids) [1] <- "orig.ident"
geneids$ensembl <- gsub("\\..*","",geneids$orig.ident)
geneids2 <- left_join(geneids, gene, by = c("ensembl" = "ensembl_gene_id"))
geneids2 <- geneids2[!duplicated(geneids2$ensembl),]
geneids3 <- left_join(geneids, geneids2, by = "ensembl")
identical(geneids3[['ensembl']],geneids[['ensembl']])
geneids3$hgnc_symbol[geneids3$hgnc_symbol==""] <- NA
geneids3 <- geneids3 %>% mutate(symbol = coalesce(hgnc_symbol,ensembl))
identical(rownames(quantLog), geneids3$orig.ident.x)
rownames(quantLog) <- geneids3$symbol
quantLog <- quantLog[!duplicated(rownames(quantLog)),]
quantLog <- quantLog-min(quantLog)

# Generate PCA, removing bottom 10% of variables based on variance
samples$Cell_Line <- as.factor(samples$Cell_Line)
samples$Disease <- as.factor(samples$Disease)
samples$Batch <- as.factor(samples$Batch)

p <- pca(quantLog, metadata = samples, removeVar = 0.1)

# Supplemental Figure 2b
biplot(p, x = "PC1",y="PC2", 
       colby = 'Batch',
       colkey = c('1' = '#D43F3AFF', 
                  '2' = '#EEA236FF',
                  '3' = '#5CB85CFF',
                  '4' = '#357EBDFF'),
       legendPosition = 'right',
       colLegendTitle = "Batch",
       axisLabSize = 10,
       encircle = F,
       lab = NULL) 

# Supplemental Figure 2c
biplot(p, x = "PC1",y="PC2", 
       colby = 'Ribo2',
       colkey = c('Lexogen RiboCop Depletion' = '#374E55FF', 
                  'NEBNext PolyA Mag Isolation' = '#B24745FF',
                  'oligo(dT)-primed RT' = '#79AF97FF'),
       legendPosition = 'right',
       colLegendTitle = "Ribosomal Depletion",
       axisLabSize = 10,
       encircle = F,
       lab = NULL) 

# With batch correction
mat <- assay(vsd)
mm <- model.matrix(~Disease, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Batch, design=mm)
assay(vsd) <- mat
quantLog <- assay(vsd)
table(rownames(samples) == colnames(quantLog))
geneids <- as.data.frame(rownames(quantLog))
colnames(geneids) [1] <- "orig.ident"
geneids$ensembl <- gsub("\\..*","",geneids$orig.ident)
geneids2 <- left_join(geneids, gene, by = c("ensembl" = "ensembl_gene_id"))
geneids2 <- geneids2[!duplicated(geneids2$ensembl),]
geneids3 <- left_join(geneids, geneids2, by = "ensembl")
identical(geneids3[['ensembl']],geneids[['ensembl']])
geneids3$hgnc_symbol[geneids3$hgnc_symbol==""] <- NA
geneids3 <- geneids3 %>% mutate(symbol = coalesce(hgnc_symbol,ensembl))
identical(rownames(quantLog), geneids3$orig.ident.x)
rownames(quantLog) <- geneids3$symbol
quantLog <- quantLog[!duplicated(rownames(quantLog)),]
quantLog <- quantLog-min(quantLog)

p <- pca(quantLog, metadata = samples, removeVar = 0.1)

# Supplemental Figure 2e
biplot(p, x = "PC1",y="PC2", 
       colby = 'Batch',
       colkey = c('1' = '#D43F3AFF', 
                  '2' = '#EEA236FF',
                  '3' = '#5CB85CFF',
                  '4' = '#357EBDFF'),
       legendPosition = 'right',
       colLegendTitle = "Batch",
       axisLabSize = 10,
       encircle = F,
       lab = NULL) 

# Figure 2e
biplot(p, x = "PC1",y="PC2", 
       colby = 'Sample',
       colkey = c('ALS1' = '#FFCC65FF', 
                  'ALS2' = '#FF9932FF',
                  'ALS3' = '#FF6619FF',
                  'ALS4' = '#FF2A00FF',
                  'CTR1' = '#002AFFFF',
                  'CTR2' = '#1965FFFF',
                  'CTR3' = '#3299FFFF',
                  'CTR4' = '#65CCFFFF',
                  'CTR5' = '#99EDFFFF'),
       legendPosition = 'right',
       colLegendTitle = "Cell Line",
       axisLabSize = 10,
       encircle = F,
       lab = NULL) 

# Figure 2f
eigencorplot(p, metavars = c('Batch',
                             'Disease',
                             'Library',
                             'Ribo'),
             components = getComponents(p, seq_len(6)),
             corMultipleTestCorrection = "BH",
             posColKey = "bottom",
             col = paletteer_d("colorBlindness::Blue2Orange12Steps"))


# Plot gene expression
symbol <- "NEFL"
symbol <- "NEFM"
symbol <- "NEFH"
symbol <- "PRPH"
symbol <- "NRGN"
symbol <- "MNX1"
symbol <- "ISL1"
symbol <- "TARDBP"

test <- cbind.data.frame("RNA" = quantLog[symbol,], 
                         "Status" = as.character(samples$Disease),
                         "Cell_Line" = as.character(samples$Sample),
                         "Batch" = as.character(samples$Batch))
ggplot(data = test, mapping = aes(x = Status, y = RNA)) +
  geom_boxplot(outlier.shape=NA, show.legend = FALSE) +
  geom_jitter(aes(color=Cell_Line), width = 0.2, size=2, alpha=0.8, show.legend = TRUE) +
  scale_x_discrete(limits=rev) +
  scale_color_manual(values=c('ALS1' = '#FFCC65FF', 
                              'ALS2' = '#FF9932FF',
                              'ALS3' = '#FF6619FF',
                              'ALS4' = '#FF2A00FF',
                              'CTR1' = '#002AFFFF',
                              'CTR2' = '#1965FFFF',
                              'CTR3' = '#3299FFFF',
                              'CTR4' = '#65CCFFFF',
                              'CTR5' = '#99EDFFFF')) +
  ylab("Normalized Expression") +
  ggtitle(bquote(paste(italic(.(symbol))))) +
  theme_minimal()

# Figure 2g
degs <- final1
keyvals <- ifelse(
  degs$log2FoldChange < -0.5, '#2873C2',
  ifelse(degs$log2FoldChange > 0.5, '#EFC001','gray85'))

keyvals[is.na(keyvals)] <- 'gray85'
names(keyvals)[keyvals == '#EFC001'] <- 'high'
names(keyvals)[keyvals == 'gray85'] <- 'mid'
names(keyvals)[keyvals == '#2873C2'] <- 'low'

table(keyvals)  

lab_italics <- paste0("italic('", degs$gene, "')")
genes <- paste0(
  "italic('",
  c("NEFH","MNX1","NRGN","NLGN1",
    "NRXN1","NRXN3","NLGN4X","GRIA2","GRIA3","GRIK1",
    "GRIN2D","GABRA1","GABRG3","CHRNB4","SLC6A2","SLC18A2",
    "DNAH7","DNAH9","WDR41","ANXA11","GAD1","GAD2","SST",
    "NRG3","BCL2","SLC32A1","SLC1A1", "LINC00682", "TNFSF10","KCNG1","CLU"),
  "')")

EnhancedVolcano(degs, 
                lab = lab_italics,
                selectLab = genes,
                colCustom = keyvals,
                legendPosition = "none",
                ylim = c(0,25),
                xlim = c(-6,6),
                axisLabSize = 12,
                colAlpha = 0.8,
                xlab = bquote(~Log[2] ~ "(fold change)"),
                ylab = bquote(~-Log[10] ~ '(' ~ italic('adjusted p-value') ~ ')'),
                x ="log2FoldChange", 
                y ="padj",
                drawConnectors = T,
                lengthConnectors = unit(0.001, "npc"),
                max.overlaps = 20,
                arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 6,
                labSize = 4,
                gridlines.major = F,
                gridlines.minor = F,
                cutoffLineWidth = 0,
                vline = 0,
                vlineType = "longdash",
                vlineCol = "black",
                titleLabSize = 18,
                subtitleLabSize = 12,
                boxedLabels = FALSE,
                parseLabels = TRUE,
                title = "",
                subtitle = bquote(italic()),
                caption = "",
                raster = TRUE)


# Variance Partitioning
library(variancePartition)
form <- ~ (1|Disease) + (1|Sample) + (1|Batch) + (1|Library) + (1|RNA) +
  (1|Ribo) + (1|Cell_Line) + seqdepth

# Run variancePartition analysis
# Supplemental Figure 2d and 2f
varPart <- fitExtractVarPartModel(quantLog, form, samples, BPPARAM=SnowParam(32))
vp <- sortCols(varPart, FUN = median)
plotVarPart( vp,label.angle = 45) +
  xlab("") +
  ylab("Variance explained (%)") +
  theme(axis.text.x = element_text(size = 10)) +
  ggtitle("Gene expression variation explained by covariates")
