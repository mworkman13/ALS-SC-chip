library(tidyverse)
library(DESeq2)
library(PCAtools)
library(RColorBrewer)
library(enrichR)
library(openxlsx)
library(EnhancedVolcano)


# Load Biomart
load("2022-09-30 Biomart Load.RData")

# Load counts and meta data and tidy up
counts <- read.csv("Chip_vs_96well_counts.csv", row.names=1)
samples <- read.csv("Chip_vs_96well_metadata.csv", row.names=1)
samples[sapply(samples, is.character)] <- lapply(samples[sapply(samples, is.character)], 
                                                 as.factor)
rownames(samples) <- make.names(rownames(samples))
colnames(counts) <- make.names(colnames(counts))


# Subset Samples 
samples <- subset(samples, Disease == "Control")
samples <- droplevels(samples)
table(samples$Disease)
table(samples$Condition)

# Overlap metadata and counts
counts <- counts[,colnames(counts) %in% rownames(samples)]
table(rownames(samples) == colnames(counts))
counts.order <- match(rownames(samples), colnames(counts))
counts <- counts[,counts.order]
identical(rownames(samples),colnames(counts))

# Relevel samples add seq depth to meta
samples$Condition <- relevel(samples$Condition, "96well")
counts <- round(counts, 0)
plot(colSums(counts) / 1e6, main = "Counts Per Million", ylim = c(0,max(colSums(counts) / 1e6)+5))
samples$seqdepth <- colSums(counts)/ 1e6
samples$Batch <- as.factor(samples$Batch)

samples$group <- ifelse(samples$Condition == "96well", "Plate", "Chip")

# Generate DESeq object
dds <- DESeqDataSetFromMatrix(counts, colData=samples, design= ~Condition)

# # Filter out low expressed genes
# keep <- rowSums(counts(dds)) >= 58
# table(keep)
# dds <- dds[keep,]  

# Alternative filtering
# Filter out genes where there are less than 5 samples with normalized counts greater than or equal to 10
dds <- estimateSizeFactors(dds)
keep <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 5
table(keep)
dds <- dds[keep,]  

# Run DEseq
dds <- DESeq(dds) 
resultsNames(dds)

# Check levels of the group factor
levels(dds$Condition)
# Static SC-chip vs 96 well
res_Stat96 <- as.data.frame(results(dds, contrast = c("Condition", "Static", "96well"), pAdjustMethod = "bonferroni"))
geneids <- as.data.frame(rownames(res_Stat96))
colnames(geneids) [1] <- "orig.ident"
geneids$ensembl <- gsub("\\..*","",geneids$orig.ident)
geneids2 <- left_join(geneids, gene, by = c("ensembl" = "ensembl_gene_id"))
geneids2 <- geneids2[!duplicated(geneids2$ensembl),]
geneids3 <- left_join(geneids, geneids2, by = "ensembl")
identical(geneids3[['ensembl']],geneids[['ensembl']])
geneids3$hgnc_symbol[geneids3$hgnc_symbol==""] <- NA
geneids3 <- geneids3 %>% mutate(symbol = coalesce(hgnc_symbol,ensembl))
identical(rownames(res_Stat96), geneids3$orig.ident.x)
res_Stat96$gene <- geneids3$symbol
# Flow SC-chip vs 96 well
res_Flow96 <- as.data.frame(results(dds, contrast = c("Condition", "Flow", "96well"), pAdjustMethod = "bonferroni"))
geneids <- as.data.frame(rownames(res_Flow96))
colnames(geneids) [1] <- "orig.ident"
geneids$ensembl <- gsub("\\..*","",geneids$orig.ident)
geneids2 <- left_join(geneids, gene, by = c("ensembl" = "ensembl_gene_id"))
geneids2 <- geneids2[!duplicated(geneids2$ensembl),]
geneids3 <- left_join(geneids, geneids2, by = "ensembl")
identical(geneids3[['ensembl']],geneids[['ensembl']])
geneids3$hgnc_symbol[geneids3$hgnc_symbol==""] <- NA
geneids3 <- geneids3 %>% mutate(symbol = coalesce(hgnc_symbol,ensembl))
identical(rownames(res_Flow96), geneids3$orig.ident.x)
res_Flow96$gene <- geneids3$symbol

# Export variance stablized data, correct for cell line
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mm <- model.matrix(~Condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Cell_Line, design=mm)
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

# Plot gene expression
# Supplemental Figure 1b-d
symbol <- "MNX1"
symbol <- "ISL1"
symbol <- "NKX6-1"

test <- cbind.data.frame("RNA" = quantLog[symbol,], 
                         "Status" = as.character(samples$Disease),
                         "Cell_Line" = as.character(samples$Cell_Line),
                         "Condition" = as.character(samples$Condition))
test$Condition <- factor(test$Condition, levels=c('Flow', 'Static', '96well'))
test$Status <- factor(test$Status, levels = c("Control","ALS"))
ggplot(data = test, mapping = aes(x = Condition, y = RNA, fill = Condition)) +
  geom_boxplot(outlier.shape=NA, show.legend = FALSE) +
  geom_jitter( width = 0.2, size=2, alpha=0.5, show.legend = FALSE) +
  scale_fill_manual(values=c("#2874C0","#888888","#69b3a2")) +
  scale_x_discrete(limits=rev) +
  ylab("Normalized Expression") +
  ggtitle(bquote(paste(italic(.(symbol)))))+
  theme_classic()


## Plot ISL1 vs MNX1
test <- cbind.data.frame("ISL1" = quantLog["ISL1",], 
                         "MNX1" = quantLog["MNX1",], 
                         "Status" = as.character(samples$Disease),
                         "Cell_Line" = as.character(samples$Cell_Line),
                         "Condition" = as.character(samples$Condition))
ggplot(test, aes(x=ISL1, y=MNX1, color=Condition)) + 
  geom_point(size=3, alpha = 1) + 
  scale_color_manual(values=c("#69b3a2","#2874C0","#888888"))+
  xlab(expression(paste(italic("ISL1")," Normalized Expression"))) +
  ylab(expression(paste(italic("MNX1")," Normalized Expression"))) +
  theme_classic() +
  geom_smooth(method='lm', se = TRUE, color = "azure4", aes(color=NULL))

cor(test$ISL1, test$MNX1)
cor.test(test$ISL1, test$MNX1)

# Supplemental Figure 1f
# Generate PCA, removing bottom 10% of variables based on variance
samples$Cell_Line <- as.factor(samples$Cell_Line)
samples$Disease <- as.factor(samples$Disease)
samples$Condition <- as.factor(samples$Condition)
levels(samples$Condition)
samples$Condition <- factor(samples$Condition, levels = c("96well","Static","Flow"))

p <- pca(quantLog, metadata = samples, removeVar = 0.1)

biplot(p, x = "PC1",y="PC2", 
       colby = 'Condition',
       legendPosition = 'right',
       colLegendTitle = "Condition",
       encircle = T,
       lab = NULL,
       axisLabSize = 12,
       colkey = c("Flow" =  "#2874C0",
                  "Static" = "#888888", 
                  "96well" = "#69b3a2"))

# Supplemental Figure 1g-k
library(ComplexHeatmap)
library(circlize)
library(dendextend)
# proliferation
genes <- c("MKI67","PCNA","TOP2A","E2F1","E2F2","MCM2","MCM3","MCM4","MCM5","MCM6",
           "CDK2","CDK4")
# Neuron
genes <- c("CHAT","STMN2","TUBB3","MAP2","NEFL","NEFM","NEFH","MAPT",
           "RBFOX3")

# Motor Neuron
genes <- c("ISL1","MNX1","NKX6-1","ISL2","OLIG2")

# Neural progenitor
genes <- c("NES","SOX2","GFAP","NCAM1","DCX")

# Synaptic markers
genes <- c("STX1A", "SNAP25", "SYT", "SV2A", "SYP",
           "DLG1", "DLG2", "DLG3", "DLG4", "SHANK1", "SHANK2", "SHANK3",
           "NLGN1", "NLGN2", "NLGN3", "NLGN4", "NLGN5",
           "NPTX2", "PPP1R9B", "NRGN")

vsd.heat <- quantLog
keep <- rownames(vsd.heat) %in% genes
vsd.heat <- vsd.heat[keep,]  
colnames(vsd.heat) <- samples$Condition
z.matheat <- t(scale(t(vsd.heat), center=TRUE, scale=TRUE))
z.matheat <- z.matheat[complete.cases(z.matheat),]
min(z.matheat)
max(z.matheat)
col_fun <- colorRampPalette(rev(brewer.pal(10, "Spectral")))(10)
ha = HeatmapAnnotation(Cell_Line = samples$Sample,
                       Condition = colnames(z.matheat),
                       annotation_name_gp= gpar(fontsize = 10),
                       gap = unit(2, "points"),
                       simple_anno_size = unit(3, "mm"),
                       show_annotation_name = TRUE,
                       annotation_name_side = "left",
                       annotation_legend_param = list(
                         Condition = list(at = c("96well","Static","Flow"),
                                          title_position = "topcenter")),
                       col = list(Cell_Line = c('CTR1' = '#3E52A1FF',
                                                'CTR2' = '#FFC737FF',
                                                'CTR3' = '#EA3228FF'),
                                  Condition = c("Flow" =  "#2874C0",
                                                "Static" = "#888888", 
                                                "96well" = "#69b3a2")))
hmap1 = Heatmap(z.matheat,
                col = col_fun,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                column_names_gp = gpar(fontsize = 8), 
                cluster_columns = TRUE,
                show_column_names = FALSE,
                show_heatmap_legend = TRUE,
                row_names_side = "left",
                row_dend_side = "right",
                top_annotation = ha,
                heatmap_legend_param = list(title = "Row Z-score", 
                                            color_bar = "continuous",
                                            title_position = "topcenter",
                                            legend_direction = "horizontal",
                                            legend_height = unit(4, "cm")))


draw(hmap1, annotation_legend_side = "right",  heatmap_legend_side = "right")

# Proliferation reordering
column_dend = as.dendrogram(hclust(dist(t(z.matheat))))
plot(column_dend)
column_dend <- rotate(set(column_dend), c(7,8,9,10,11,12,12,14,15,1,2,3,4,5,6))
plot(column_dend)

# NPC reordering
column_dend = as.dendrogram(hclust(dist(t(z.matheat))))
plot(column_dend)
column_dend <- rotate(set(column_dend), c(15,14,13,12,11,10,9,8,6,7,5,4,3,2,1))
plot(column_dend)

# Neuron reordering
column_dend = as.dendrogram(hclust(dist(t(z.matheat))))
plot(column_dend)
column_dend <- rotate(set(column_dend), c(1,2,3,10,11,12,13,14,15,4,5,6,7,8,9))
plot(column_dend)

# MN reordering
column_dend = as.dendrogram(hclust(dist(t(z.matheat))))
plot(column_dend)

# Synaptic reordering
column_dend = as.dendrogram(hclust(dist(t(z.matheat))))
plot(column_dend)
column_dend <- rotate(set(column_dend), c(1,2,3,10,11,12,13,14,15,4,5,6,7,8,9))
plot(column_dend)

hmap1 = Heatmap(z.matheat,
                col = col_fun,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                column_names_gp = gpar(fontsize = 8), 
                cluster_columns = column_dend,
                show_column_names = FALSE,
                show_heatmap_legend = TRUE,
                row_names_side = "right",
                row_dend_side = "left",
                top_annotation = ha,
                heatmap_legend_param = list(title = "Row Z-score", 
                                            color_bar = "continuous",
                                            title_position = "topcenter",
                                            legend_direction = "horizontal",
                                            legend_height = unit(4, "cm")))


draw(hmap1, annotation_legend_side = "right",  heatmap_legend_side = "right")


# GSEA Supplemental Figure 1l and 1m
library(msigdbr)
library(fgsea)

msigdbr_species()

m_df <- msigdbr(species = "Homo sapiens", category = "C5")
m_df <- subset(m_df, gs_subcat != "HPO")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Rank DEGS by logFC and p-value
# Static SC-chip vs 96 well
# Supplemental Figure 1l
degs <- res_Stat96
degs$padj <- ifelse(degs$padj == 0, 1e-305, degs$padj)
degs$ranknum <- -log10(degs$padj) * degs$log2FoldChange
degs <- degs[complete.cases(degs), ]
gsea_degs <- data.frame(gene = degs$gene, logFC = degs$ranknum)
gsea_degs <- gsea_degs %>%
  arrange(desc(logFC))

ranks<- deframe(gsea_degs)

head(ranks)
all(is.finite(ranks))

fgseaRes<- fgsea(fgsea_sets, 
                 stats = ranks,
                 minSize  = 10)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(padj)

fgsea.up <- subset(fgseaRes, NES > 0) %>% arrange(padj)
fgsea.down <- subset(fgseaRes, NES < 0) %>% arrange(padj)
fgsea.plot <- rbind(fgsea.up[1:8,], fgsea.down[1:8,])
ggplot(fgsea.plot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  scale_fill_manual(values = c("grey70", "springgreen4")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA Static SC-chip vs. 96-well") + 
  theme_minimal()

# Flow SC-chip vs 96 well
# Supplemental Figure 1m
degs <- res_Flow96
degs$padj <- ifelse(degs$padj == 0, 1e-305, degs$padj)
degs$ranknum <- -log10(degs$padj) * degs$log2FoldChange
degs <- degs[complete.cases(degs), ]
gsea_degs <- data.frame(gene = degs$gene, logFC = degs$ranknum)
gsea_degs <- gsea_degs %>%
  arrange(desc(logFC))

ranks<- deframe(gsea_degs)

head(ranks)
all(is.finite(ranks))

fgseaRes<- fgsea(fgsea_sets, 
                 stats = ranks,
                 minSize  = 10)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(padj)

fgsea.up <- subset(fgseaRes, NES > 0) %>% arrange(padj)
fgsea.down <- subset(fgseaRes, NES < 0) %>% arrange(padj)
fgsea.plot <- rbind(fgsea.up[1:8,], fgsea.down[1:8,])
ggplot(fgsea.plot, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  scale_fill_manual(values = c("grey70", "springgreen4")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA Flow SC-chip vs. 96-well") + 
  theme_minimal()


