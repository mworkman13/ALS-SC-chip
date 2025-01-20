# Flow vs Static SC-chip analysis
library(tidyverse)
library(DESeq2)
library(PCAtools)
library(RColorBrewer)
library(enrichR)
library(openxlsx)
library(forcats)
library(EnhancedVolcano)


# Load Biomart
load("2022-09-30 Biomart Load.RData")

# Load counts and meta data and tidy up
counts <- read.csv("Flow_vs_Static_counts.csv", row.names=1)
samples <- read.csv("Flow_vs_Static_metadata.csv", row.names=1)
samples[sapply(samples, is.character)] <- lapply(samples[sapply(samples, is.character)], 
                                                 as.factor)
rownames(samples) <- make.names(rownames(samples))
colnames(counts) <- make.names(colnames(counts))
identical(rownames(samples),colnames(counts))

# Relevel samples add seq depth to meta
samples$Treatment <- relevel(samples$Treatment, "Static")
counts <- round(counts, 0)
plot(colSums(counts) / 1e6, main = "Counts Per Million")
samples$seqdepth <- colSums(counts)/ 1e6
samples$Cell_Line <- as.factor(samples$Cell_Line)

# Generate DESeq object
dds <- DESeqDataSetFromMatrix(counts, colData=samples, design= ~Cell_Line+Treatment)

# Filter out low expressed genes
dds <- estimateSizeFactors(dds)
keep <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
table(keep)
dds <- dds[keep,]  

# Run DEseq
dds <- DESeq(dds) 
resultsNames(dds)

# Run diff expression
res1 <- results(dds, name = "Treatment_Flow_vs_Static", pAdjustMethod = "bonferroni")
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
table(final1$padj<0.05)


# Export variance stablized data
quantLog <- assay(vst(dds))
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
p <- pca(quantLog, metadata = samples, removeVar = 0.1)

# Figure 1g
biplot(p, x = "PC1",y="PC2", 
       colby = 'Treatment',
       legendPosition = 'right',
       pointSize = 5,
       colkey = c('Static' = '#888888',
                  'Flow' = '#2874C0'),
       shape = "Cell_Line",
       colLegendTitle = "Condition",
       lab = NULL) 


degs <- subset(final1, padj < 6e-39)
degs_up <- subset(degs, log2FoldChange > 0.25)
degs_down <- subset(degs, log2FoldChange < -0.25)
degs <- final1

# Figure 1h
keyvals <- ifelse(
  degs$log2FoldChange < -0.5, '#888888',
  ifelse(degs$log2FoldChange > 0.5, '#2874C0','gray85'))

keyvals[is.na(keyvals)] <- 'gray85'
names(keyvals)[keyvals == '#2874C0'] <- 'high'
names(keyvals)[keyvals == 'gray85'] <- 'mid'
names(keyvals)[keyvals == '#888888'] <- 'low'

table(keyvals)  

lab_italics <- paste0("italic('", degs$gene, "')")
genes <- paste0(
  "italic('",
  c("BNIP3","FAM162A","PGK1","NDNF","NGFR","PRKCE","DDIT4",
    "HK2","ALDOC","SLC16A3","MT3","PCNA","CDC45","ORC1","CDC6",
    "MCM3","MCM4","MCM5","BRCA1","POLE","MCM2","MCM6",
    "CDK2"),
  "')")


EnhancedVolcano(degs, 
                lab = lab_italics,
                selectLab = genes,
                colCustom = keyvals,
                legendPosition = "none",
                ylim = c(0,215),
                xlim = c(-8,8),
                axisLabSize = 12,
                colAlpha = 0.8,
                xlab = bquote(~Log[2] ~ "(fold change)"),
                ylab = bquote(~-Log[10] ~ '(' ~ italic('adjusted p-value') ~ ')'),
                x ="log2FoldChange", 
                y ="padj",
                drawConnectors = T,
                lengthConnectors = unit(0.01, "npc"),
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
                boxedLabels = F,
                parseLabels = TRUE,
                title = "",
                subtitle = bquote(italic('Flow vs. Static Spinal-cord chips')),
                caption = "",
                raster = T)

# Figure 1i and 1j
setEnrichrSite("Enrichr")

dbs <- c("GO_Molecular_Function_2021", 
         "GO_Cellular_Component_2021", 
         "GO_Biological_Process_2021",
         'KEGG_2021_Human',
         'Reactome_2022')

enriched <- enrichr(degs_up$gene, dbs)
i = 3
plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]


enriched <- enrichr(degs_down$gene, dbs)
i = 3
plotEnrich(enriched[[i]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]






