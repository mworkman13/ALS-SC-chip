library(Seurat)
library(tidyverse)
library(future)
library(RColorBrewer)
library(scCustomize)
library(enrichR)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(biomaRt)
library(pheatmap)
library(ggrastr)
library(corrplot)

# Setup parallel
plan()
plan("multicore", workers = 24)
plan()
options(future.globals.maxSize = 64 * 1024 ^ 3)

# Set Seed
set.seed(13)

# Load Seurat object
data <- readRDS("20221222_ALS_chip_integrated_all_2_15PCs_res.0.2_seurat_after_RNAnorm.rds")

# Plot UMAP
DimPlot(data, reduction = 'umap', label = TRUE) + NoLegend()
DefaultAssay(data) <- "RNA"

new.cluster.ids <- c("NSCs", 
                     "Neurons", 
                     "NPCs",
                     "Neurons", 
                     "Neurons", 
                     "Neurons", 
                     "BMEC-like cells",
                     "Neurons",
                     "PHOX2B MNs", 
                     "MNX1 MNs",
                     "Neurons",
                     "NPCs",
                     "NSCs",
                     "Neurons")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
levels(data) <- c("NSCs","NPCs","Neurons","MNX1 MNs",
                  "PHOX2B MNs","BMEC-like cells")

# Figure 3b
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
p <- DimPlot(data, label = FALSE, raster = TRUE, pt.size = 2) + NoLegend() + NoAxes()
LabelClusters(p, id = "ident", size = 4, fontface = "bold", color = "grey30")

# Figure 3c
DimPlot(data, label = FALSE, raster = TRUE, pt.size = 2, split.by = "genotype") + NoLegend() + NoAxes()

# Supplemental Figure 3a
DimPlot(data, label = FALSE, raster = TRUE, pt.size = 2, split.by = "orig.ident") + NoLegend() + NoAxes()

# QC metrics
qc <- data.frame("UMI" = tapply(data$nCount_RNA, data$orig.ident, median),
                 "Features" = tapply(data$nFeature_RNA, data$orig.ident, median),
                 "Mito" = tapply(data$percent.mt, data$orig.ident, median),
                 "Ribo" = tapply(data$percent.ribo, data$orig.ident, median),
                 "Nuclei" = table(data$orig.ident))

# Figure 3e
heat.marks <- c("MKI67","TOP2A","VIM",
                "SOX2","NES","MAP2","TUBB3",
                "SNAP25","RBFOX3",
                "ISL1","CHAT","MNX1","PHOX2B")

avgexp = AggregateExpression(data, return.seurat = T, add.ident = 'orig.ident')

DoHeatmap(avgexp, features = heat.marks, size = 3, draw.lines = FALSE) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(10)) +
  NoLegend()



# Subset MNs
data$cell_type <- Idents(data)
mn.sub <- subset(data, cell_type == "MNX1 MNs" | cell_type == "PHOX2B MNs")
DefaultAssay(mn.sub) <- "integrated"
mn.sub <- FindNeighbors(mn.sub, dims = 1:30, reduction = "pca")
mn.sub <- FindClusters(mn.sub, resolution = 0.11, cluster.name = "unintegrated_clusters")
mn.sub <- RunUMAP(mn.sub, dims = 1:30, reduction = "pca", min.dist = 0.4,
                  spread = 0.5, reduction.name = "umap.unintegrated")

# Figure 4a
DimPlot(mn.sub, reduction = "umap.unintegrated", group.by = c("seurat_clusters"))

# Figure 4b
DimPlot(mn.sub, reduction = "umap.unintegrated", split.by = "genotype") + NoAxes()

# Normalize and Scale
DefaultAssay(mn.sub) <- "RNA"
mn.sub <- NormalizeData(mn.sub)
mn.sub <- FindVariableFeatures(mn.sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mn.sub)
mn.sub <- ScaleData(mn.sub, features = all.genes)

# Figure 4c
FeaturePlot_scCustom(mn.sub, features = "CHAT", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "SLC18A3", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "SLC5A7", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4d
FeaturePlot_scCustom(mn.sub, features = "MNX1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "PHOX2B", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4e
FeaturePlot_scCustom(mn.sub, features = "BCL6", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "ESRRG", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "GFRA1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4f
FeaturePlot_scCustom(mn.sub, features = "PRPH", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "PHOX2A", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4g
FeaturePlot_scCustom(mn.sub, features = "MECOM", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4h
FeaturePlot_scCustom(mn.sub, features = "FOXP1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()

# Figure 4i
FeaturePlot_scCustom(mn.sub, features = "TBX20", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()
FeaturePlot_scCustom(mn.sub, features = "NKX6-1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap.unintegrated") + NoAxes() + NoLegend()


# Figure 4j
degs_orig = data.frame()
Idents(mn.sub) <- mn.sub$seurat_clusters
levels(mn.sub)
for (i in levels(mn.sub)){
  # vector output
  model <- FindMarkers(mn.sub, test.use = "wilcox",
                       ident.1 = "ALS", ident.2 = "CTL", 
                       group.by = 'genotype',
                       subset.ident = i,
                       min.pct = 0.1,
                       logfc.threshold = 0,
                       verbose = TRUE)
  model <- rownames_to_column(model, "gene")
  model$cluster <- paste0(i)
  
  # add vector to a dataframe
  df <- data.frame(model)
  degs_orig <- rbind(degs_orig,df)
}

degs <- subset(degs_orig, p_val_adj < 0.05)
degs_up <- subset(degs, avg_log2FC > 0)
degs_down <- subset(degs, avg_log2FC < 0)


degs_orig$volcano <- ifelse(degs_orig$p_val_adj<0.05, degs_orig$cluster, "ns")
cols <- c("0" = "#F8766D", 
          "1" = "#B79F00",
          "2" = "#00BA38", 
          "3" = "#00BFC4", 
          "4" = "#619CFF", 
          "5" = "#F564E3", 
          "ns" = "grey") 
sizes <- c("0" = 2, 
           "1" = 2,
           "2" = 2, 
           "3" = 2, 
           "4" = 2, 
           "5" = 2, 
           "ns" = 1.5) 
alphas <- c("0" = 2, 
            "1" = 2,
            "2" = 2, 
            "3" = 2, 
            "4" = 2, 
            "5" = 2, 
            "ns" = 1)
degs_volcano <- degs_orig[sample(nrow(degs_orig)),]
ggplot(degs_volcano, aes(x = avg_log2FC,
                         y = -log10(p_val_adj),
                         fill = volcano,    
                         size = volcano,
                         alpha = volcano)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-6, 6, 2)),       
                     limits = c(-6, 6)) +
  theme_classic()

# Figure 4k
setEnrichrSite("Enrichr")
dbs <- c("GO_Molecular_Function_2023", 
         "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

up.ALS <- enrichr(degs_up$gene, dbs)
down.ALS <- enrichr(degs_down$gene, dbs)

i = 2
plotEnrich(up.ALS[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
up.ALS.list <- up.ALS[[i]]

plotEnrich(down.ALS[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
down.ALS.list <- down.ALS[[i]]

## Monocle3 pseudotime analysis
# Convert to cell data set and cluster
cds <- as.cell_data_set(data)
cds <- estimate_size_factors(cds)
cds <- cluster_cells(cds, resolution = 0.00001, verbose = TRUE)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

plot_cells(cds, color_cells_by = "cell_type", show_trajectory_graph = FALSE)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == c(1,2))
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds, use_partition = FALSE, verbose = TRUE)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

# Choose starting point
cds <- order_cells(cds)

# Figure 3f
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE, rasterize = TRUE) + NoAxes()

# Add metadata
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

# Figure 3g
pseudo_genes <- c("ISL1")
pseudo_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pseudo_genes,
                          colData(cds)$monocle3_partitions %in% c("1","2")]
plot_genes_in_pseudotime(pseudo_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5,
                         vertical_jitter = 0.1)
FeaturePlot_scCustom(data, features = "ISL1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap") + NoAxes() + NoLegend()

pseudo_genes <- c("CHAT")
pseudo_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pseudo_genes,
                          colData(cds)$monocle3_partitions %in% c("1","2")]
plot_genes_in_pseudotime(pseudo_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5,
                         vertical_jitter = 0.1)
FeaturePlot_scCustom(data, features = "CHAT", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap") + NoAxes() + NoLegend()

# Supplemental Figure 3b
pseudo_genes <- c("SOX2")
pseudo_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pseudo_genes,
                          colData(cds)$monocle3_partitions %in% c("1","2")]
plot_genes_in_pseudotime(pseudo_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5,
                         vertical_jitter = 0.1)
FeaturePlot_scCustom(data, features = "SOX2", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap") + NoAxes() + NoLegend()

# Supplemental Figure 3c
pseudo_genes <- c("SNAP25")
pseudo_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pseudo_genes,
                          colData(cds)$monocle3_partitions %in% c("1","2")]
plot_genes_in_pseudotime(pseudo_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5,
                         vertical_jitter = 0.1)
FeaturePlot_scCustom(data, features = "SNAP25", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap") + NoAxes() + NoLegend()

# Supplemental Figure 3d
pseudo_genes <- c("MAP2")
pseudo_lineage_cds <- cds[rowData(cds)$gene_short_name %in% pseudo_genes,
                          colData(cds)$monocle3_partitions %in% c("1","2")]
plot_genes_in_pseudotime(pseudo_lineage_cds,
                         color_cells_by="cell_type",
                         min_expr=0.5,
                         vertical_jitter = 0.1)
FeaturePlot_scCustom(data, features = "MAP2", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     reduction = "umap") + NoAxes() + NoLegend()

# Differential expression across the pseudotime
# cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=32)

genes <- c("VIM","NKX6-1","SOX2","NES","MKI67",
           "HES5","OLIG2","PAX6","NEUROG2","MAP2","HOXD3","LHX4","HOXA5","SNAP25","ONECUT1","TUBB3",
           "RBFOX3","STMN2","ACHE","LHX3","NEFH","CHAT","SLC5A7","SLC18A3","ISL1","MNX1","PHOX2B")

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

# Supplemental Figure 3e
Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)


# Andersen et al developing spinal cord correlation
pasca_skeletal <- read.csv("pasca_skeletal.csv", sep="")
pasca_visceral <- read.csv("pasca_visceral.csv", sep="")

skeletal_mod <- pasca_skeletal$gene
visceral_mod <- pasca_visceral$gene

data <- AddModuleScore(data,
                       features = list(skeletal_mod),
                       name="skeletal")
data <- AddModuleScore(data,
                       features = list(visceral_mod),
                       name="visceral")

# Plot scores
meta.data.module <- data.frame(data@meta.data)

# Supplemental Figure 3g
FeaturePlot_scCustom(data, features = "skeletal1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     na_cutoff = 0.15654) + NoAxes()
# Supplemental Figure 3h
ggplot(meta.data.module, aes(x=cell_type,y=skeletal1,fill = cell_type))+
  geom_boxplot(outlier.shape = NA)+
  theme_minimal()

# Supplemental Figure 3i
FeaturePlot_scCustom(data, features = "visceral1", raster = FALSE,
                     colors_use = viridis_plasma_light_high,
                     na_cutoff = 0.34127) + NoAxes()
# Supplemental Figure 3j
ggplot(meta.data.module, aes(x=cell_type,y=visceral1,fill = cell_type))+
  geom_boxplot(outlier.shape = NA)+
  theme_minimal()

# histogram
hist(data$skeletal1)
hist(data$visceral1)
# QQ-plot
library(car)
library(multcomp)
qqPlot(data$skeletal1, id = FALSE)
qqPlot(data$visceral1, id = FALSE)

res_aov <- aov(skeletal1 ~ cell_type,
               data = meta.data.module)

summary(res_aov)

# Tukey HSD test:
post_test <- glht(res_aov,
                  linfct = mcp(cell_type = "Tukey")
)

summary(post_test)
par(mar = c(3, 14, 3, 3))
plot(post_test)

res_aov <- aov(visceral1 ~ cell_type,
               data = meta.data.module)

summary(res_aov)

# Tukey HSD test:
post_test <- glht(res_aov,
                  linfct = mcp(cell_type = "Tukey")
)

summary(post_test)
par(mar = c(3, 14, 3, 3))
plot(post_test)

# Correlation to mouse and human in vivo spinal cord
# Code adapted from Gautier...Gitler, Neuron 2023
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(humanx)
}

# Reference datasets from http://spinalcordatlas.org/more_info.html
gautier_mns <- readRDS('GSE228778_gautier_mns.rds')
mouse <- readRDS("lepichon_gitler_integrated.RDS")

# Plot UMAP
mouse$motor_neurons <- NA
mouse$motor_neurons <- ifelse(mouse$seurat_clusters %in% c("0", "3", "1", "23","24"), 
                              "Skeletal Motor Neurons", mouse$motor_neurons)
mouse$motor_neurons <- ifelse(mouse$seurat_clusters %in% c("22", "14", "15", "11","2","5"), 
                              yes = "Cholinergic Interneurons", no = mouse$motor_neurons)
mouse$motor_neurons <- ifelse(mouse$seurat_clusters %in% c("10", "12", "13", "16","8",
                                                           "6","18","9","21","4","19",
                                                           "20","7","17"), 
                              yes = "Visceral Motor Neurons", no = mouse$motor_neurons)

mouse_mn_markers <- rownames(mouse)
mouse_markers_to_human <- convertMouseGeneList(mouse_mn_markers)

mouse_mn_avg_exp <- AggregateExpression(mouse, assays = c('RNA'), group.by = "motor_neurons") %>% as.data.frame()
colnames(mouse_mn_avg_exp) <- c("Cholinergic Interneurons-M",
                                "Skeletal Motor Neurons-M",
                                "Visceral Motor Neurons-M")
mouse_mn_avg_exp <- mouse_mn_avg_exp[,-1]

unique_mouse <- mouse_markers_to_human %>% group_by(MGI.symbol) %>% filter(n() == 1)
unique_human <- mouse_markers_to_human %>% group_by(HGNC.symbol) %>% filter(n() == 1)
mouse_genes_to_human_dedup <- inner_join(unique_mouse, unique_human)

gautier_mn_avg_exp <- AggregateExpression(gautier_mns, assays = c('RNA'), group.by = "motor_neuron") %>% as.data.frame()
colnames(gautier_mn_avg_exp) <- c("MNs - H (Gautier)")

gautier_subtype_avg_exp <- AggregateExpression(gautier_mns, assays = c('RNA'), group.by = "motor_neuron_subtype") %>% as.data.frame()
colnames(gautier_subtype_avg_exp) <- c("Alpha MNs - H (Gautier)", "Gamma MNs - H (Gautier)")

dimn_avg_exp <- AggregateExpression(data, assays = c('RNA'), group.by = "cell_type") %>% as.data.frame()

mouse_mn_avg_exp$MGI.symbol <- rownames(mouse_mn_avg_exp)
mn_avg_exp_combined <- inner_join(mouse_mn_avg_exp, mouse_genes_to_human_dedup, by = "MGI.symbol")

gautier_subtype_avg_exp$HGNC.symbol <- rownames(gautier_subtype_avg_exp)
mn_avg_exp_combined <- right_join(gautier_subtype_avg_exp, mn_avg_exp_combined, by = "HGNC.symbol")

gautier_mn_avg_exp$HGNC.symbol <- rownames(gautier_mn_avg_exp)
mn_avg_exp_combined <- right_join(gautier_mn_avg_exp, mn_avg_exp_combined, by = "HGNC.symbol")

dimn_avg_exp$HGNC.symbol <- rownames(dimn_avg_exp)
mn_avg_exp_combined <- right_join(dimn_avg_exp, mn_avg_exp_combined, by = "HGNC.symbol")

# Set remaining NA to 0
mn_avg_exp_combined[is.na(mn_avg_exp_combined)] = 0
mn_avg_exp_combined <- dplyr::select(mn_avg_exp_combined, -MGI.symbol)


# Only keep mouse marker genes that mapped to human orthologs
mn_marker_avg_exp_combined <- inner_join(mn_avg_exp_combined, mouse_markers_to_human, by = "HGNC.symbol")
row.names(mn_marker_avg_exp_combined) <- mn_marker_avg_exp_combined$HGNC.symbol
mouse.mn.corr <- read.csv("MouseMNmarkersgenes_GautierTableS1E.csv", sep="")
mn_marker_avg_exp_combined <- inner_join(mn_avg_exp_combined, mouse.mn.corr, by = "HGNC.symbol")

mn_marker_avg_exp_combined <- dplyr::select(mn_marker_avg_exp_combined, -HGNC.symbol)
cor.mat <- as.matrix(cor(as.matrix(mn_marker_avg_exp_combined), method = "pearson"))


col.pal <- rev(brewer.pal(n=8, name="RdYlBu"))
col.pal <- colorRampPalette(c(col.pal))

# matrix of the correlations
corrplot(cor.mat, method="color",  
         col=col.pal(200),
         col.lim = c(0,1),
         is.corr = F,
         type="upper", order="original", 
         hclust.method = "ward.D2",
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=TRUE)


# 96 well plate integration and correlation
# Load 96 well data and normalize
plate2NK3 = Read10X_h5('./plate/3_hg38_2NK3-96_outs/filtered_feature_bc_matrix.h5')
plate6DHM = Read10X_h5('./plate/3_hg38_6DHM-96_outs/filtered_feature_bc_matrix.h5')
plate7TN6 = Read10X_h5('./plate/3_hg38_7TN6-96_outs/filtered_feature_bc_matrix.h5')
plate8FF6 = Read10X_h5('./plate/3_hg38_8FF6-96_outs/filtered_feature_bc_matrix.h5')
plateEDi029 = Read10X_h5('./plate/hg38_188-96_outs/filtered_feature_bc_matrix.h5')
plate9FP9 = Read10X_h5('./plate/hg38_9FP9-96_outs/filtered_feature_bc_matrix.h5')



plate2NK3 <- CreateSeuratObject(counts=plate2NK3, project = "plate2NK3", min.cells = 0, min.features=0)
plate2NK3 <- AddMetaData(plate2NK3, metadata = "Plate", col.name = "culture")
plate2NK3 <- AddMetaData(plate2NK3, metadata = "3", col.name = "batch")
plate2NK3 <- AddMetaData(plate2NK3, metadata = "CTR", col.name = "status")
plate2NK3 <- AddMetaData(plate2NK3, metadata = "2NK3", col.name = "cell_line")
plate2NK3 <- AddMetaData(plate2NK3, metadata = "CTR7", col.name = "label")

plate6DHM <- CreateSeuratObject(counts=plate6DHM, project = "plate6DHM", min.cells = 0, min.features=0)
plate6DHM <- AddMetaData(plate6DHM, metadata = "Plate", col.name = "culture")
plate6DHM <- AddMetaData(plate6DHM, metadata = "3", col.name = "batch")
plate6DHM <- AddMetaData(plate6DHM, metadata = "CTR", col.name = "status")
plate6DHM <- AddMetaData(plate6DHM, metadata = "6DHM", col.name = "cell_line")
plate6DHM <- AddMetaData(plate6DHM, metadata = "CTR8", col.name = "label")

plate7TN6 <- CreateSeuratObject(counts=plate7TN6, project = "plate7TN6", min.cells = 0, min.features=0)
plate7TN6 <- AddMetaData(plate7TN6, metadata = "Plate", col.name = "culture")
plate7TN6 <- AddMetaData(plate7TN6, metadata = "3", col.name = "batch")
plate7TN6 <- AddMetaData(plate7TN6, metadata = "ALS", col.name = "status")
plate7TN6 <- AddMetaData(plate7TN6, metadata = "7TN6", col.name = "cell_line")
plate7TN6 <- AddMetaData(plate7TN6, metadata = "ALS7", col.name = "label")

plate8FF6 <- CreateSeuratObject(counts=plate8FF6, project = "plate8FF6", min.cells = 0, min.features=0)
plate8FF6 <- AddMetaData(plate8FF6, metadata = "Plate", col.name = "culture")
plate8FF6 <- AddMetaData(plate8FF6, metadata = "3", col.name = "batch")
plate8FF6 <- AddMetaData(plate8FF6, metadata = "ALS", col.name = "status")
plate8FF6 <- AddMetaData(plate8FF6, metadata = "8FF6", col.name = "cell_line")
plate8FF6 <- AddMetaData(plate8FF6, metadata = "ALS8", col.name = "label")

plate9FP9 <- CreateSeuratObject(counts=plate9FP9, project = "plate9FP9", min.cells = 0, min.features=0)
plate9FP9 <- AddMetaData(plate9FP9, metadata = "Plate", col.name = "culture")
plate9FP9 <- AddMetaData(plate9FP9, metadata = "2", col.name = "batch")
plate9FP9 <- AddMetaData(plate9FP9, metadata = "ALS", col.name = "status")
plate9FP9 <- AddMetaData(plate9FP9, metadata = "9FP9", col.name = "cell_line")
plate9FP9 <- AddMetaData(plate9FP9, metadata = "ALS1", col.name = "label")

plateEDi029 <- CreateSeuratObject(counts=plateEDi029, project = "plateEDi029", min.cells = 0, min.features=0)
plateEDi029 <- AddMetaData(plateEDi029, metadata = "Plate", col.name = "culture")
plateEDi029 <- AddMetaData(plateEDi029, metadata = "2", col.name = "batch")
plateEDi029 <- AddMetaData(plateEDi029, metadata = "CTR", col.name = "status")
plateEDi029 <- AddMetaData(plateEDi029, metadata = "EDi029", col.name = "cell_line")
plateEDi029 <- AddMetaData(plateEDi029, metadata = "CTR4", col.name = "label")


# Merge Seurat objects
plate <- merge(plate2NK3, y = c(plate6DHM,plate7TN6,
                                plate8FF6,plate9FP9, plateEDi029))
plate <- JoinLayers(plate,  assay = "RNA")

# Remove objects from memory
rm(plate2NK3,plate6DHM,plate7TN6,plate8FF6,plateEDi029,plate9FP9)


##Calculate the percentage of mitochondrial genes and store it in percent.mito
plate[["percent.mt"]] <- PercentageFeatureSet(plate, pattern = "^MT-")
plate[["percent.ribo"]] <- PercentageFeatureSet(plate, pattern = "^RP[SL]")

# Normalize and scale data with SCTransform
plate <- SCTransform(plate, verbose = TRUE)

# Run PCA
plate <- RunPCA(plate, verbose = TRUE)
# Determine percent of variation associated with each PC
pct <- plate[["pca"]]@stdev / sum(plate[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2

# UMAP embedding
plate <- RunUMAP(plate, dims = 1:co2, min.dist = 0.3, verbose = TRUE)

# Cluster assignment
plate <- FindNeighbors(plate, dims = 1:co2,
                       verbose = TRUE)
plate <- FindClusters(plate, resolution = 0.2,
                      verbose = TRUE)

# Plot UMAP
DimPlot(plate, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

# Log normalize and scale
DefaultAssay(plate) <- "RNA"
plate <- NormalizeData(plate)
plate <- FindVariableFeatures(plate, selection.method = "vst")
all.genes <- rownames(plate)
plate <- ScaleData(plate, features = all.genes)

# Map query
DefaultAssay(data)<- "RNA"
DefaultAssay(plate) <- "RNA"
data.anchors <- FindTransferAnchors(reference = data, query = plate, dims = 1:30,
                                    reference.reduction = "pca",
                                    normalization.method = "LogNormalize")
# set embeddings
data[["umap.new"]] <- CreateDimReducObject(embeddings = data[["umap"]]@cell.embeddings, key = "UMAPnew_")

# set UMAP models
umap.new.model <- list()
umap.new.model$n_epochs <- 500
umap.new.model$alpha <-1
umap.new.model$method <- "umap"
umap.new.model$negative_sample_rate <- 5
umap.new.model$n_neighbors <- 30
umap.new.model$gamma <- 1
umap.new.model$approx_pow <- 0
umap.new.model$metric$cosine <- list()
umap.new.model$embedding <- data[["umap.new"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap.new.model$a <- ab_param["a"]
umap.new.model$b <- ab_param["b"]
data[["umap.new"]]@misc$model <- umap.new.model

plate.query <- MapQuery(anchorset = data.anchors, reference = data, query = plate,
                        refdata = list(cell_type = "cell_type"), 
                        reference.reduction = "pca", reduction.model = "umap.new")

plate.query$predicted.cell_type <- factor(plate.query$predicted.cell_type,
                                          levels = c("NSCs","NPCs","Neurons","MNX1 MNs",
                                                     "PHOX2B MNs","BMEC-like cells"))
DimPlot(plate.query, reduction = "ref.umap", group.by = "predicted.cell_type", label = FALSE,
        label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

p1 <- DimPlot(data, reduction = "umap", label = F, raster = F, pt.size = 0.1) +
  NoLegend() + NoAxes() + ggtitle("SC-chip (reference labels)")
p2 <- DimPlot(plate.query, reduction = "ref.umap", group.by = "predicted.cell_type",
              label = F, raster = FALSE, pt.size = 0.1) +
  NoLegend() + NoAxes() + ggtitle("96 well plate (predicted labels)")

p1+p2


plate.query$orig.ident <- factor(plate.query$orig.ident,
                                 levels = c("plateEDi029",
                                            "plate2NK3",
                                            "plate6DHM",
                                            "plate9FP9",
                                            "plate7TN6",
                                            "plate8FF6"))
DimPlot(plate.query, label = FALSE, raster = TRUE, pt.size = 2, split.by = "orig.ident") + NoLegend() + NoAxes()


plate.query$cell_type <- plate.query$predicted.cell_type
Idents(plate.query) <- plate.query$cell_type
Idents(plate.query)
Idents(data)
data$culture <- "Chip"

# Supplemental Figure 4e
DefaultAssay(plate.query) <- "RNA"
heat.marks <- c("MKI67","TOP2A","VIM",
                "SOX2","NES","MAP2","TUBB3",
                "SNAP25","RBFOX3",
                "ISL1","CHAT","MNX1","PHOX2B")


avgexp = AggregateExpression(plate.query, return.seurat = T,
                             group.by = c("cell_type", "orig.ident"))
avgexp <- subset(avgexp, cell_type != "BMEC-like cells")
DoHeatmap(avgexp, features = heat.marks, size = 2, draw.lines = FALSE) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(10)) +
  NoLegend()

# Supplemental Figure 4f
FeaturePlot_scCustom(plate.query, features = "ISL1", raster = FALSE,
                     colors_use = viridis_plasma_light_high, order = T,
                     reduction = "ref.umap") + NoAxes()


degs = data.frame()
for (i in levels(plate.query)){
  # vector output
  model <- FindMarkers(plate.query, test.use = "wilcox",
                       ident.1 = "ALS", ident.2 = "CTR", 
                       group.by = 'status',
                       subset.ident = i,
                       min.pct = 0.1,
                       logfc.threshold = 0.25,
                       verbose = TRUE)
  model <- rownames_to_column(model, "gene")
  model$cluster <- paste0(i)
  
  # add vector to a dataframe
  df <- data.frame(model)
  degs <- rbind(degs,df)
}

degs <- subset(degs, p_val_adj < 0.05)
degs_up <- subset(degs, avg_log2FC > 0)
degs_down <- subset(degs, avg_log2FC < 0)


setEnrichrSite("Enrichr")
dbs <- c("GO_Molecular_Function_2023", 
         "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

up <- degs_up[degs_up$cluster=="MNX1 MNs" | degs_up$cluster=="PHOX2B MNs",]


enriched <- enrichr(up$gene, dbs)
i = 2
plotEnrich(enriched[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]

down <- degs_down[degs_down$cluster=="MNX1 MNs" | degs_down$cluster=="PHOX2B MNs",]
enriched <- enrichr(down$gene, dbs)
i = 3
plotEnrich(enriched[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]




all.samples <- merge(data, y= plate.query, collapse = TRUE)

Idents(all.samples)
DefaultAssay(all.samples) <- "RNA"
all.samples <- JoinLayers(all.samples)

all.samples <- NormalizeData(all.samples)
all.samples <- FindVariableFeatures(all.samples, selection.method = "vst")
all.genes <- rownames(all.samples)
all.samples <- ScaleData(all.samples, features = all.genes)
table(all.samples$cell_type)

degs = data.frame()
for (i in levels(all.samples)){
  # vector output
  model <- FindMarkers(all.samples, test.use = "wilcox",
                       ident.1 = "ALS", ident.2 = "CTR", 
                       group.by = 'status',
                       subset.ident = i,
                       min.pct = 0.1,
                       logfc.threshold = 0.25,
                       verbose = TRUE)
  model <- rownames_to_column(model, "gene")
  model$cluster <- paste0(i)
  
  # add vector to a dataframe
  df <- data.frame(model)
  degs <- rbind(degs,df)
}

degs <- subset(degs, p_val_adj < 0.05)
degs_up <- subset(degs, avg_log2FC > 0)
degs_down <- subset(degs, avg_log2FC < 0)


up <- degs_up[degs_up$cluster=="MNX1 MNs" | degs_up$cluster=="PHOX2B MNs",]
enriched <- enrichr(up$gene, dbs)
i = 3
plotEnrich(enriched[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]

down <- degs_down[degs_down$cluster=="MNX1 MNs" | degs_down$cluster=="PHOX2B MNs",]
enriched <- enrichr(down$gene, dbs)
i = 3
plotEnrich(enriched[[i]], showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value")
temp <- enriched[[i]]


# Andersen et al correlation
all.samples <- AddModuleScore(all.samples,
                              features = list(visceral_mod),
                              name="visceral")
VlnPlot(all.samples, features = "visceral1", pt.size = 0, slot = "data",
        group.by = "cell_type", split.by = "culture")

all.samples <- AddModuleScore(all.samples,
                              features = list(skeletal_mod),
                              name="skeletal")
VlnPlot(all.samples, features = "skeletal1", pt.size = 0, slot = "data",
        group.by = "cell_type", split.by = "culture")


meta.data.module <- data.frame(all.samples@meta.data)
meta.data.module <- subset(meta.data.module, cell_type == "MNX1 MNs" | cell_type == "PHOX2B MNs")

ggplot(meta.data.module, aes(x=cell_type,y=skeletal1,fill = culture))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=alpha(c("#2E58A4","#B69D71"),0.8)) +
  theme_minimal()


ggplot(meta.data.module, aes(x=cell_type,y=visceral1,fill = culture))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=alpha(c("#2E58A4","#B69D71"),0.8)) +
  theme_minimal()

ggplot(meta.data.module, aes(x=cell_type,y=visceral1,fill = culture))+
  geom_violin(outlier.shape = NA)+
  scale_fill_manual(values=alpha(c("#2E58A4","#B69D71"),0.8)) +
  theme_minimal()
