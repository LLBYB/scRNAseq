library(Seurat)
library(monocle3)
library(harmony)
library(patchwork)
library(ggplot2)
library(scater)
library(ComplexHeatmap)

data <- Read10X(data.dir = "/home/3002/")
rrna.genes <- rownames(data)[grep("^rr.$",rownames(data))]
rrna.genes

data <- CreateSeuratObject(counts = data, project = "MTB", min.cells = 3, min.features = 10)
# data$batch <- "batch1"
data[["percent.rrna"]] <- PercentageFeatureSet(data, pattern = "^rr.$")
data <- data[!rownames(data) %in% c(rrna.genes) , ]
data;data
data <- subset(data, subset = nCount_RNA >200 & nCount_RNA < 3000)

n_UMI <- c(15,50, 100, 200, 400)
data_subsets <- list()

for (umi in n_UMI) {
  subset_name <- paste0("UMI>", umi)
  data_subsets[[subset_name]] <- subset(
    x = data,
    subset = nCount_RNA > umi
  )
}

umi_results <- data.frame(
  Threshold = character(),
  CellCount = integer(),
  Median_UMI = numeric(),
  Mean_UMI = numeric(),
  stringsAsFactors = FALSE
)


for (threshold in names(data_subsets)) {
  current_subset <- data_subsets[[threshold]]
  cell_count <- ncol(current_subset)
  median_umi <- median(current_subset$nCount_RNA, na.rm = TRUE)
  Mean_umi <- mean(current_subset$nCount_RNA, na.rm = TRUE)
  umi_results <- rbind(
    umi_results,
    data.frame(
      Threshold = threshold,
      CellCount = cell_count,
      Median_UMI = median_umi,
      Mean_UMI = Mean_umi,
      stringsAsFactors = FALSE
    )
  )
}
umi_results

gene_results <- data.frame(
  Threshold = character(),
  CellCount = integer(),
  Median_gene = numeric(),
  Mean_gene = numeric(),
  stringsAsFactors = FALSE
)


for (threshold in names(data_subsets)) {

  current_subset <- data_subsets[[threshold]]
  cell_count <- ncol(current_subset)
  Median_gene <- median(current_subset$nFeature_RNA, na.rm = TRUE)
  Mean_gene <- mean(current_subset$nFeature_RNA, na.rm = TRUE)
  gene_results <- rbind(
    gene_results,
    data.frame(
      Threshold = threshold,
      CellCount = cell_count,
      Median_gene = Median_gene,
      Meang_ene = Mean_gene,
      stringsAsFactors = FALSE
    )
  )
}
gene_results


combined_data <- data.frame()
for (threshold in names(data_subsets)) {
  current_data <- data_subsets[[threshold]]
  temp_df <- data.frame(
    UMI = current_data$nCount_RNA,
    Threshold = factor(threshold, levels = names(data_subsets)),
    Cell = colnames(current_data)
  )

  combined_data <- rbind(combined_data, temp_df)
}


combined_violin <- ggplot(combined_data, aes(x = Threshold, y = UMI, fill = Threshold)) +
  geom_violin(
    draw_quantiles = 0.5,
    scale = "width",
    trim = FALSE,
    alpha = 0.7
  ) +
  labs(
    title = "UMI",
    x = "UMI",
    y = "UMIs per cell"
  ) +
  theme_minimal(base_family = "sans") + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5,lineend = "square"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(trans = "log10",    
                     breaks = 10^(0:5),   
                     labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                     limits = c(1, 100000)  
                     ,expand = expansion(mult = c(0, 0.05))
  ) 

print(combined_violin)


data <- subset(data, subset = nCount_RNA >200 & nCount_RNA < 3000)
#Ctrl:200<UMI<3000;INH:200<UMI<3000;PZA:160<UMI<3000;EMB:90<UMI<3000;RIF:40<UMI<3000;RFP:40<UMI<3000;INH_7d:200<UMI<3000;
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data <- SCTransform(
  data,
  vars.to.regress = "percent.rrna",
)

data <- RunPCA(object = data)
# data <- RunHarmony(data, group.by.vars = "batch") #多组样品去除批次效应
VizDimLoadings(data,dims = 1:10,reduction = "pca",nfeatures = 30)
DimPlot(data,reduction = "pca",
        pt.size = 1,  
        cols = "#4A90E2")
data <- FindNeighbors(object = data, dims = 1:5)
data <- FindClusters(object = data, resolution = 0.2)
data <- RunUMAP(object = data, dims = 1:10)
DimPlot(object = data, reduction = 'umap',label = TRUE,pt.size = 0.8)
data.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

VlnPlot(data, features = c("pks13"), pt.size = 0)
FeaturePlot(
  sce, 
  features = c("pks13"), 
  pt.size = 0.5,
  cols = c( "blue","yellow","red" )  
)

Cellratio <- prop.table(table(Idents(data), data$stim), margin = 2)
Cellratio <- as.data.frame(Cellratio)
ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 



expression_matrix <- GetAssayData(data, assay = 'SCT',slot = 'counts')
cell_metadata <- data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data@assays$SCT@data))
rownames(gene_annotation) <- gene_annotation[,1]

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 10)
plot_pc_variance_explained(cds)
cds <- align_cds(cds, alignment_group = "orig.ident")
cds <- reduce_dimension(cds, preprocess_method = "PCA") 
cds <- cluster_cells(cds, resolution=1e-5) 
cds <- learn_graph(
  cds,
  use_partition = TRUE,
  close_loop = FALSE,
  learn_graph_control = list(
    euclidean_distance_ratio = 3,  
    geodesic_distance_ratio = 1/3,   
    minimal_branch_len = 10,    
    orthogonal_proj_tip = TRUE,   
    prune_graph = TRUE,             
    ncenter =100               
  ),
  verbose = TRUE
)

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = TRUE, group_label_size = 3.5)

get_earliest_principal_node <- function(cds, time_bin="1"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds <- order_cells(cds, 
                   root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups=FALSE,label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)

plot_cells(
  cds,
  genes = "pks13", 
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE,
  cell_size = 0.7
)


plot_cells(
  cds,
  genes = "pks13",
  color_cells_by = "celltype",
  label_cell_groups = TRUE,
  cell_size = 0.7
) 

markergenes=c("desA1","desA2","mce1F","mce1C","ppsC","fadD29","PPE20",
              "fbpC","umaA","pcaA","hadC","mmaA4","fabD","kasA","kasB",
              "accD6","fas","accD5","accE5","accA3",
              "accD4","pks13","fadD32","fbpA","embA","embB","mmpL3","mmpL10")

cds_subset <- cds[rowData(cds)$gene_short_name %in% markergenes, ]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "celltype")

pseudotime_vals <- pseudotime(cds)
ordered_cells <- names(sort(pseudotime_values, na.last = TRUE))
expr_matrix <- exprs(cds)[rowData(cds)$gene_short_name %in%markergenes, ordered_cells]
expr_matrix_dense <- as.matrix(expr_matrix)
expr_scaled <- t(scale(t(expr_matrix)))
expr_scaled[is.na(expr_scaled)] <- 0
gene_clusters <- kmeans(expr_scaled, centers =2)$cluster 
gene_order <- order(gene_clusters)
expr_scaled_ordered <- expr_scaled[gene_order, ]

gene_clusters_ordered <- gene_clusters[gene_order]
ha_time <- HeatmapAnnotation(
  Pseudotime = anno_barplot(
    pseudotime(cds)[ordered_cells],
    gp = gpar(fill = "lightblue"),
    border = FALSE,
    axis = TRUE,
    axis_param = list(side = "left", labels_rot = 0)
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  height = unit(1, "cm")
)

row_ha <- rowAnnotation(
  Cluster = as.factor(gene_clusters_ordered),
  col = list(Cluster = setNames(brewer.pal(3, "Set2"), 1:3)), 
  show_legend = TRUE,
  show_annotation_name = TRUE
)


Heatmap(
  expr_scaled_ordered,
  name = "Expression\nZ-score", 
  col = colorRamp2(c(-2, 0,2), c("#00bfff", "#E2E2E2","#FF0000")),
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  show_row_names = TRUE, 
  show_column_names = FALSE, 
  row_names_gp = gpar(fontsize = 8), 
  column_title = paste("Cells ordered by pseudotime"),
  top_annotation = ha_time,
  use_raster = FALSE,
  raster_quality = 2,
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_height = unit(3, "cm")
  )
)








