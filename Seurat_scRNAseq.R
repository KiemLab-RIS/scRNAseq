#
# 04/27/2022
# Kiem Lab
# Fred Hutchinson Cancer Center
# Seurat processing for scRNAseq data with cellular barcoding (CellTag).
# 
#
library('Seurat')
library('ggplot2')
library('dplyr')

basedir <- 'Filepath for the Seurat filtered counts matrix.'
resdir <- 'Filepath to save the result files.'
setwd(basedir)

# The raw fastq files were processed using the Cell Ranger Suiteto generate the filtered_feature_bc_matrix.
# Load both the data sets and initialize the Seurat object with the raw data.
data1  <- Read10X(data.dir = paste0(basedir,"17088_CD34GFP_1/outs/filtered_feature_bc_matrix/"))
data.1 <- CreateSeuratObject(counts = data1, project = "17088_CD34GFP_1", min.cells = 3, min.features = 200)
data.1

data2  <- Read10X(data.dir = paste0(basedir,"17088_CD34GFP_2/outs/filtered_feature_bc_matrix/"))
data.2 <- CreateSeuratObject(counts = data2, project = "17088_CD34GFP_2", min.cells = 3, min.features = 200) 
data.2

# Merge the two data sets.
scadden.combined <- merge(data.1, y = data.2, add.cell.ids = c("D1", "D2"), project = "17088_CD34GFP_merged") 
scadden.combined
head(colnames(scadden.combined))
table(scadden.combined$orig.ident)
scadden.combined@meta.data[, "protocol"] <- "CD34"

# Read the processed CellTag barcodes and add them to the Seurat object as meta data.
# The csv file has the cell barcode and the associated clone ID.
cd34.cloneID <- as.sparse(read.csv(file = paste0(basedir,"Apr_22_17088_CD34_Combined_10x_clone_barcodes.csv"), 
                                   sep = ",", header = TRUE, row.names = 1,stringsAsFactors = FALSE))
head(cd34.cloneID)
dim(cd34.cloneID)
dim(scadden.combined)
barcodes = scadden.combined@assays$RNA@counts@Dimnames[[2]]

#
# Match the barcodes in the csv to the  barcodes of the Seurat object.
# Some cells filtered out at start of file.
# 
clone.names = colnames(cd34.cloneID)
clone.names = paste0(clone.names,"-1")

# Match the cell barcodes and the dimensions of both the data sets.
clone.names[1]
barcodes[1]
t = clone.names %in% barcodes
table(t)
cd34.cloneID.a =  cd34.cloneID[,t]
class(cd34.cloneID)
class(cd34.cloneID.a)
nm = names(cd34.cloneID.a)
cd34.clone.m = matrix(data = cd34.cloneID.a,nrow = 1,ncol=length(cd34.cloneID.a),dimnames = list("cloneID",nm))
class(cd34.clone.m)
dim(cd34.clone.m)

# Create a new metadat slot cloneID to add the corresponding clone ID to the data.
scadden.combined$cloneID = cd34.clone.m[1,]


# QC metrics
scadden.combined[["percent.mt"]] <- PercentageFeatureSet(scadden.combined, pattern = "^NDUF")
scadden.combined[["percent.rps"]] <- PercentageFeatureSet(scadden.combined, pattern = "^RPS")

# Print the violin plot for the metrics.
pdf(paste0(resdir,"17088_CD34GFP_merged_QC_pre_filtering.pdf"))
VlnPlot(scadden.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps"), ncol = 4)
dev.off()

# We filter the cells based on the unique feature counts and mitochondrial counts.
dim(scadden.combined)
scadden.combined <- subset(scadden.combined, subset = nFeature_RNA > 200 & nCount_RNA > 2000 &
                             percent.mt > 0.1 & percent.mt < 1.1 &
                             percent.rps > 4 & percent.rps < 17) 
dim(scadden.combined)

pdf(paste0(resdir,"17088_CD34GFP_merged_QC_post_filtering.pdf"))
VlnPlot(scadden.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps"), ncol = 4)
dev.off()

# Normalization
scadden.combined <- NormalizeData(scadden.combined)

# Feature selection
scadden.combined <- FindVariableFeatures(scadden.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scadden.combined), 10)
top10

# Scaling data
all.genes <- rownames(scadden.combined)
scadden.combined <- ScaleData(scadden.combined, features = all.genes)

# Perform linear dimensional reduction.
scadden.combined <- RunPCA(scadden.combined, features = VariableFeatures(object = scadden.combined), npcs = 50) 
print(scadden.combined[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(scadden.combined, reduction = "pca")
DimHeatmap(scadden.combined, dims = 1, cells = 500, balanced = TRUE)

# Cell-cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

scadden.combined <- CellCycleScoring(scadden.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(scadden.combined[[]])

# Regress the cell-cycle scores during data scaling.
scadden.combined$CC.Difference <- scadden.combined$S.Score - scadden.combined$G2M.Score
scadden.combined <- ScaleData(scadden.combined, vars.to.regress = c("CC.Difference"), features = rownames(scadden.combined))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
scadden.combined <- RunPCA(scadden.combined, features = VariableFeatures(scadden.combined), nfeatures.print = 10, npcs = 50)

# Cluster the cells
scadden.combined <- FindNeighbors(scadden.combined, dims = 1:10)
# The resolution of 0.4 was decided after trying out many different variables that best represented the data.
scadden.combined <- FindClusters(scadden.combined, resolution = 0.4) 
head(Idents(scadden.combined), 5)

# Run non-linear dimensional reduction
scadden.combined <- RunUMAP(scadden.combined, dims = 1:10)

pdf(paste0(resdir,"17088_CD34GFP_merged_umap.pdf")) 
DimPlot(scadden.combined, reduction = "umap", pt.size = 0.8)
dev.off()

# Cell sorting function for Feature plots.
# Cell sorting ensures that the cells with higher counts are plotted over the ones with zero counts.
gene.cells <- WhichCells(object = scadden.combined, expression = protocol %in% c("CD34"), invert = FALSE)
raw_data <- GetAssayData(object = scadden.combined, slot = "counts")
gene.raw.data <- as.matrix(raw_data[, WhichCells(object = scadden.combined, cells = gene.cells)])
scadden_gene_names <- rownames(gene.raw.data)

gene_list <- list("APOC1","CD52","CD79B","CKS1B","DNTT","ELANE","FCER1A","GATA2",
                  "JCHAIN","KLF1","MLLT3","MPO","TYMS","EGFP") 

ft_plots <- function(x){
  
  rnum <- (which(scadden_gene_names == x))
  g_rawdata <- gene.raw.data[rnum,]
  g_rawdata_sort <- sort(g_rawdata)
  g_celluse <- names(g_rawdata_sort)
  print(x)
  pdf(paste0(resdir,"17088_CD34GFP_merged_0.8_",x,".pdf"))
  print(FeaturePlot(scadden.combined, features = x, reduction = 'spring', pt.size = 0.8, cells = g_celluse))#
  dev.off()
}

feature_plots <- lapply(gene_list,ft_plots)

# Finding differentially expressed features (cluster biomarkers)
# Find markers for every cluster compared to all remaining cells, report only the positive ones.
# The markers along with the resolution was used to inform the clusters.
scadden.combined.markers <- FindAllMarkers(scadden.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scadden.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#
# Spring plot
#
pdf(paste0(resdir,"17088_CD34GFP_merged_spring.pdf"))
DimPlot(scadden.combined, pt.size = 0.8, reduction = "spring", label = T)
dev.off()


# Plots
theme_bare <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines")
)
pdf(paste0(resdir,"17088_CD34GFP_merged_umap_cell_type_scoring.pdf"))
DimPlot(scadden.combined, label = T, pt.size = 0.8, reduction = "umap", 
        group.by = "CellType", cols = c('HSC'='mediumorchid1','Erythroid 1'='darkgoldenrod3',
                                        'Erythroid 2'='darkolivegreen4','Myeloid 1'='turquoise3', 'Myeloid 2'='cadetblue2',
                                        'Basophil'='indianred1','Lymphoid'='springgreen3','NULL'='cornflowerblue','Undecided'='maroon1'))
dev.off()

# EGFP feature plot
pdf(paste0(resdir,"Feature_plot_17088_CD34GFP_merged_GFP.pdf"))
FeaturePlot(scadden.combined, features = "EGFP", reduction = 'umap', pt.size = 0.8)
dev.off()


## Visualizing the CellTag clonal data.
for (i in 1:50) {
  
  cells = names(which(scadden.combined$cloneID == i))
  print(paste0("Clone_",i))
  pdf(paste0(resdir,"clone_",i,".pdf"))
  print(DimPlot(scadden.combined, reduction = 'umap', pt.size = 0.8, 
                cells.highlight = cells, sizes.highlight = 2.0) + NoLegend())
  dev.off()
}


# Save the processed Seurat object with the clone IDs, cell type scoring and the spring plot dimensional reduction values.
saveRDS(scadden.combined, file = paste0(resdir,"17088_CD34GFP_merged_CC_0.4.rds"))
