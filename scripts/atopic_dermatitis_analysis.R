############################################################
# GSE147424 - Atopic dermatitis scRNA-seq

setwd("~/Desktop/GSE147424_RAW")

set.seed(1234)

############################################################
# 0. Libraries
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scDblFinder)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(CellChat)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(harmony)

############################################################
# 1. Output folders
############################################################

dir.create("plots", showWarnings = FALSE)
dir.create("Doublet_results", showWarnings = FALSE)
dir.create("DEG_GO_KEGG_by_celltype", showWarnings = FALSE)
dir.create("Global_DEG_GO_KEGG", showWarnings = FALSE)
dir.create("CellChat_results", showWarnings = FALSE)

dir.create("Metascape_inputs", showWarnings = FALSE)
dir.create("Metascape_inputs/By_celltype", showWarnings = FALSE)
dir.create("Metascape_inputs/Global", showWarnings = FALSE)
############################################################
# 2. Helper functions
############################################################

save_gg <- function(plot, filename, width = 8, height = 6) {
  print(plot)
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height
  )
}

save_heatmap <- function(ht, filename, width = 10, height = 8) {
  ComplexHeatmap::draw(ht)
  pdf(filename, width = width, height = height)
  ComplexHeatmap::draw(ht)
  dev.off()
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]", "_", x)
}

run_enrichment <- function(
    deg_table,
    prefix,
    outdir,
    showCategory = 15
) {
  
  deg_sig <- deg_table %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  
  message(prefix, " significant DEG count: ", nrow(deg_sig))
  
  if (nrow(deg_sig) < 5) {
    message(prefix, ": fewer than 5 significant DEGs, skipping GO/KEGG")
    return(NULL)
  }
  
  entrez <- bitr(
    deg_sig$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  entrez_ids <- unique(entrez$ENTREZID)
  
  message(prefix, " ENTREZ count: ", length(entrez_ids))
  
  if (length(entrez_ids) < 5) {
    message(prefix, ": too few ENTREZ IDs, skipping GO/KEGG")
    return(NULL)
  }
  
  ##########################################################
  # GO Biological Process
  ##########################################################
  
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    
    write.csv(
      as.data.frame(ego),
      file = file.path(outdir, paste0(prefix, "_GO_BP.csv")),
      row.names = FALSE
    )
    
    p_go <- dotplot(ego, showCategory = showCategory) +
      ggtitle(paste("GO BP -", prefix))
    
    save_gg(
      p_go,
      filename = file.path(outdir, paste0(prefix, "_GO_BP.pdf")),
      width = 9,
      height = 7
    )
    
  } else {
    message(prefix, ": no significant GO terms")
  }
  
  ##########################################################
  # KEGG
  ##########################################################
  
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    
    write.csv(
      as.data.frame(ekegg),
      file = file.path(outdir, paste0(prefix, "_KEGG.csv")),
      row.names = FALSE
    )
    
    p_kegg <- dotplot(ekegg, showCategory = showCategory) +
      ggtitle(paste("KEGG -", prefix))
    
    save_gg(
      p_kegg,
      filename = file.path(outdir, paste0(prefix, "_KEGG.pdf")),
      width = 9,
      height = 7
    )
    
  } else {
    message(prefix, ": no significant KEGG terms")
  }
}

############################################################
# 3. Load object
############################################################

combined <- readRDS("~/Desktop/GSE147424_RAW/combined.rds")

dim(combined)
table(combined$sample_id)
table(combined$condition)
head(combined@meta.data)

############################################################
# 4. QC
############################################################

# Bu veri setinde gerçek MT genleri görünmüyor.
grep("^MT-", rownames(combined), value = TRUE, ignore.case = TRUE)

# percent.mt güvenilir olmadığı için kaldırıyoruz.
if ("percent.mt" %in% colnames(combined@meta.data)) {
  combined$percent.mt <- NULL
}

p_qc_sample <- VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "sample_id",
  ncol = 2,
  pt.size = 0.2
)

save_gg(
  p_qc_sample,
  "plots/QC_before_filter_by_sample.pdf",
  width = 12,
  height = 6
)

p_qc_condition <- VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "condition",
  ncol = 2,
  pt.size = 0.2
)

save_gg(
  p_qc_condition,
  "plots/QC_before_filter_by_condition.pdf",
  width = 10,
  height = 6
)

summary(combined$nFeature_RNA)
summary(combined$nCount_RNA)

cat("Before filtering:", ncol(combined), "cells\n")

############################################################
# 5. QC filtering
############################################################

combined_filtered <- subset(
  combined,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000
)

cat("After nFeature filtering:", ncol(combined_filtered), "cells\n")

table(combined_filtered$sample_id)
table(combined_filtered$condition)

p_qc_after <- VlnPlot(
  combined_filtered,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "sample_id",
  ncol = 2,
  pt.size = 0.2
)

save_gg(
  p_qc_after,
  "plots/QC_after_filter_by_sample.pdf",
  width = 12,
  height = 6
)

p_qc_after_condition <- VlnPlot(
  combined_filtered,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "condition",
  ncol = 2,
  pt.size = 0.2
)

save_gg(
  p_qc_after_condition,
  "plots/QC_after_filter_by_condition.pdf",
  width = 10,
  height = 6
)

p_scatter <- FeatureScatter(
  combined_filtered,
  feature1 = "nFeature_RNA",
  feature2 = "nCount_RNA"
) 

save_gg(
  p_scatter,
  "plots/QC_scatter_nFeature_vs_nCount.pdf",
  width = 7,
  height = 6
)

############################################################
# 6. Doublet detection with scDblFinder
############################################################

set.seed(1234)

cat(
  "Doublet analizi öncesi toplam hücre:",
  ncol(combined_filtered),
  "\n"
)

# Not:
# Bu veri clean.data / processed matrix gibi görünüyor.
# Yine de pipeline gereği scDblFinder çalıştırıyoruz.
# Sonuçları dikkatli yorumlamak gerekir.

sce_doublet <- as.SingleCellExperiment(
  combined_filtered,
  assay = "RNA"
)

sce_doublet <- scDblFinder(
  sce_doublet,
  samples = "sample_id",
  verbose = TRUE
)

combined_filtered$doublet_score <- sce_doublet$scDblFinder.score
combined_filtered$doublet_class <- sce_doublet$scDblFinder.class

table(combined_filtered$doublet_class)

table(
  combined_filtered$sample_id,
  combined_filtered$doublet_class
)

table(
  combined_filtered$condition,
  combined_filtered$doublet_class
)

############################################################
# Doublet summary tables
############################################################

doublet_summary_sample <- as.data.frame(
  table(
    sample_id = combined_filtered$sample_id,
    doublet_class = combined_filtered$doublet_class
  )
)

write.csv(
  doublet_summary_sample,
  file = "Doublet_results/Doublet_summary_by_sample.csv",
  row.names = FALSE
)

doublet_summary_condition <- as.data.frame(
  table(
    condition = combined_filtered$condition,
    doublet_class = combined_filtered$doublet_class
  )
)

write.csv(
  doublet_summary_condition,
  file = "Doublet_results/Doublet_summary_by_condition.csv",
  row.names = FALSE
)

############################################################
# Doublet score violin plot
############################################################

p_doublet_score <- VlnPlot(
  combined_filtered,
  features = "doublet_score",
  group.by = "sample_id",
  pt.size = 0.1
)

save_gg(
  p_doublet_score,
  "plots/Doublet_score_by_sample.pdf",
  width = 12,
  height = 6
)

############################################################
# Doublet / singlet fraction barplot
############################################################

doublet_plot_data <- combined_filtered@meta.data %>%
  dplyr::count(
    sample_id,
    doublet_class,
    name = "n"
  ) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    fraction = n / sum(n)
  ) %>%
  dplyr::ungroup()

write.csv(
  doublet_plot_data,
  file = "Doublet_results/Doublet_fraction_by_sample.csv",
  row.names = FALSE
)

p_doublet_fraction <- ggplot(
  doublet_plot_data,
  aes(
    x = sample_id,
    y = fraction,
    fill = doublet_class
  )
) +
  geom_col() +
  theme_classic() +
  labs(
    x = "Sample",
    y = "Cell fraction",
    fill = "Doublet class",
    title = "Doublet / singlet fraction by sample"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_gg(
  p_doublet_fraction,
  "plots/Doublet_fraction_by_sample.pdf",
  width = 10,
  height = 6
)

############################################################
# Keep singlets only
############################################################

combined_filtered <- subset(
  combined_filtered,
  subset = doublet_class == "singlet"
)

cat(
  "Doublet çıkarıldıktan sonra toplam hücre:",
  ncol(combined_filtered),
  "\n"
)

table(combined_filtered$sample_id)
table(combined_filtered$condition)

saveRDS(
  combined_filtered,
  file = "GSE147424_AD_QC_singlets.rds"
)

############################################################
# 7. Normalize / Variable features / Scale
############################################################

combined_filtered <- NormalizeData(
  combined_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

combined_filtered <- FindVariableFeatures(
  combined_filtered,
  selection.method = "vst",
  nfeatures = 2500
)

top20 <- head(VariableFeatures(combined_filtered), 20)

p_var <- VariableFeaturePlot(combined_filtered)

p_var_labeled <- LabelPoints(
  plot = p_var,
  points = top20,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

save_gg(
  p_var,
  "plots/VariableFeaturePlot.pdf",
  width = 7,
  height = 6
)

save_gg(
  p_var_labeled,
  "plots/VariableFeaturePlot_top20_labeled.pdf",
  width = 7,
  height = 6
)

save_gg(
  p_var + p_var_labeled,
  "plots/VariableFeaturePlot_side_by_side.pdf",
  width = 14,
  height = 6
)

combined_filtered <- ScaleData(
  combined_filtered,
  features = VariableFeatures(combined_filtered),
  vars.to.regress = "nCount_RNA"
)

############################################################
# 8. PCA / Harmony / UMAP / tSNE / Clustering
############################################################

combined_filtered <- RunPCA(
  combined_filtered,
  features = VariableFeatures(combined_filtered),
  verbose = TRUE
)

p_elbow <- ElbowPlot(
  combined_filtered,
  ndims = 40
)

save_gg(
  p_elbow,
  "plots/ElbowPlot_40PCs.pdf",
  width = 7,
  height = 6
)

############################################################
# PCA heatmaps
############################################################

p_heat_1_10 <- DimHeatmap(
  combined_filtered,
  dims = 1:10,
  cells = 500,
  balanced = TRUE,
  ncol = 5
)

save_gg(
  p_heat_1_10,
  "plots/PCA_DimHeatmap_PC1_to_PC10.pdf",
  width = 15,
  height = 11
)

############################################################
# PCA plots
############################################################

p_pca_sample <- DimPlot(
  combined_filtered,
  reduction = "pca",
  group.by = "sample_id",
  pt.size = 0.3
)

p_pca_condition <- DimPlot(
  combined_filtered,
  reduction = "pca",
  group.by = "condition",
  pt.size = 0.3
)

save_gg(
  p_pca_sample + p_pca_condition,
  "plots/PCA_sample_condition_side_by_side.pdf",
  width = 16,
  height = 6
)

############################################################
# Harmony + clustering
############################################################

combined_filtered <- RunHarmony(
  object = combined_filtered,
  group.by.vars = "sample_id"
)

dims_use <- 1:20

combined_filtered <- FindNeighbors(
  combined_filtered,
  reduction = "harmony",
  dims = dims_use
)

combined_filtered <- FindClusters(
  combined_filtered,
  resolution = 0.5
)

combined_filtered <- RunUMAP(
  combined_filtered,
  reduction = "harmony",
  dims = dims_use
)

combined_filtered <- RunTSNE(
  combined_filtered,
  reduction = "harmony",
  dims = dims_use,
  check_duplicates = FALSE
)

############################################################
# UMAP / tSNE
############################################################

p_umap_condition <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.3
)

p_umap_sample <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "sample_id",
  pt.size = 0.3
)

p_umap_cluster <- DimPlot(
  combined_filtered,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.3
)

save_gg(
  p_umap_condition + p_umap_sample + p_umap_cluster,
  "plots/UMAP_condition_sample_cluster_combined.pdf",
  width = 21,
  height = 6
)

p_tsne_condition <- DimPlot(
  combined_filtered,
  reduction = "tsne",
  group.by = "condition",
  pt.size = 0.3
)

p_tsne_cluster <- DimPlot(
  combined_filtered,
  reduction = "tsne",
  label = TRUE,
  pt.size = 0.3
)

save_gg(
  p_tsne_condition + p_tsne_cluster,
  "plots/tSNE_condition_cluster_side_by_side.pdf",
  width = 16,
  height = 6
)

table(combined_filtered$seurat_clusters, combined_filtered$sample_id)
table(combined_filtered$seurat_clusters, combined_filtered$condition)

saveRDS(
  combined_filtered,
  file = "GSE147424_AD_filtered_clustered.rds"
)

############################################################
# 9. Cluster markers
############################################################

cluster_markers <- FindAllMarkers(
  combined_filtered,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25
)

write.csv(
  cluster_markers,
  file = "GSE147424_cluster_markers.csv",
  row.names = FALSE
)

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top_markers,
  file = "GSE147424_top10_markers_per_cluster.csv",
  row.names = FALSE
)

############################################################
# 10. Manual skin marker check
############################################################

skin_markers <- c(
  "KRT14", "KRT5", "KRT1", "KRT10", "KRT16", "KRT6A",
  "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",
  "CD3D", "CD3E", "TRAC", "IL7R", "CD4", "CD8A",
  "NKG7", "GNLY", "GZMB", "PRF1",
  "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN",
  "LYZ", "CD68", "C1QA", "C1QB", "LST1",
  "FCER1A", "CD1A", "CD207", "CLEC10A", "LAMP3",
  "PECAM1", "VWF", "CLDN5",
  "RGS5", "ACTA2", "TAGLN"
)

skin_markers_present <- skin_markers[
  skin_markers %in% rownames(combined_filtered)
]

p_dot <- DotPlot(
  combined_filtered,
  features = skin_markers_present
) +
  RotatedAxis()

save_gg(
  p_dot,
  "plots/Skin_marker_DotPlot.pdf",
  width = 14,
  height = 7
)

p_feature <- FeaturePlot(
  combined_filtered,
  features = skin_markers_present[1:min(12, length(skin_markers_present))],
  ncol = 3
)

print(p_feature)

ggsave(
  filename = "plots/Skin_marker_FeaturePlot_first12.png",
  plot = p_feature,
  width = 12,
  height = 10,
  dpi = 300
)


############################################################
# 11. SingleR annotation - HPCA
############################################################

refHPCA <- HumanPrimaryCellAtlasData()

sce_test <- GetAssayData(
  object = combined_filtered,
  assay = "RNA",
  layer = "data"
)

prediction_HPCA_main <- SingleR(
  test = sce_test,
  ref = refHPCA,
  clusters = Idents(combined_filtered),
  labels = refHPCA$label.main
)

prediction_HPCA_fine <- SingleR(
  test = sce_test,
  ref = refHPCA,
  clusters = Idents(combined_filtered),
  labels = refHPCA$label.fine
)

combined_filtered$HPCA_main <- prediction_HPCA_main$labels[
  Idents(combined_filtered)
]

combined_filtered$HPCA_fine <- prediction_HPCA_fine$labels[
  Idents(combined_filtered)
]

p_hpca_main <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "HPCA_main",
  label = TRUE,
  pt.size = 0.4
) +
  NoLegend()

p_hpca_fine <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "HPCA_fine",
  label = TRUE,
  pt.size = 0.4
) +
  NoLegend()

save_gg(
  p_hpca_main + p_hpca_fine,
  "plots/UMAP_HPCA_main_fine_side_by_side.pdf",
  width = 18,
  height = 7
)

table(combined_filtered$HPCA_main)
table(combined_filtered$HPCA_fine)

############################################################
# 12. SingleR annotation - Monaco immune
############################################################

refMonaco <- MonacoImmuneData()

prediction_Monaco_main <- SingleR(
  test = sce_test,
  ref = refMonaco,
  clusters = Idents(combined_filtered),
  labels = refMonaco$label.main
)

prediction_Monaco_fine <- SingleR(
  test = sce_test,
  ref = refMonaco,
  clusters = Idents(combined_filtered),
  labels = refMonaco$label.fine
)

combined_filtered$Monaco_main <- prediction_Monaco_main$labels[
  Idents(combined_filtered)
]

combined_filtered$Monaco_fine <- prediction_Monaco_fine$labels[
  Idents(combined_filtered)
]

p_monaco_main <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "Monaco_main",
  label = TRUE,
  pt.size = 0.4
) +
  NoLegend()

p_monaco_fine <- DimPlot(
  combined_filtered,
  reduction = "umap",
  group.by = "Monaco_fine",
  label = TRUE,
  pt.size = 0.4
) +
  NoLegend()

save_gg(
  p_monaco_main + p_monaco_fine,
  "plots/UMAP_Monaco_main_fine_side_by_side.pdf",
  width = 18,
  height = 7
)

saveRDS(
  combined_filtered,
  file = "GSE147424_AD_with_celltypes.rds"
)

############################################################
# 13. Cell-type-specific DEG + GO/KEGG
############################################################

seurat_obj <- combined_filtered

celltype_col <- "HPCA_main"
condition_col <- "condition"

meta <- seurat_obj@meta.data
celltypes <- na.omit(unique(meta[[celltype_col]]))

for (ct in celltypes) {
  
  message("Processing cell type: ", ct)
  
  cells_ct <- rownames(meta)[
    meta[[celltype_col]] == ct &
      !is.na(meta[[celltype_col]])
  ]
  
  obj_ct <- subset(seurat_obj, cells = cells_ct)
  
  tab <- table(obj_ct@meta.data[[condition_col]])
  print(tab)
  
  if (!all(c("Lesional", "Non_lesional") %in% names(tab))) {
    message(ct, ": missing Lesional or Non_lesional")
    next
  }
  
  if (any(tab[c("Lesional", "Non_lesional")] < 10)) {
    message(ct, ": too few cells")
    next
  }
  
  Idents(obj_ct) <- obj_ct@meta.data[[condition_col]]
  
  deg <- FindMarkers(
    obj_ct,
    ident.1 = "Lesional",
    ident.2 = "Non_lesional",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  deg$gene <- rownames(deg)
  
  prefix <- paste0(
    "Lesional_vs_Nonlesional_",
    safe_name(ct)
  )
  
  write.csv(
    deg,
    file = file.path(
      "DEG_GO_KEGG_by_celltype",
      paste0(prefix, "_DEG.csv")
    ),
    row.names = FALSE
  )
  
  ############################################################
  # Metascape input - cell-type specific
  ############################################################
  
  deg_sig <- deg %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  
  if (nrow(deg_sig) > 0) {
    
    genes_all <- unique(deg_sig$gene)
    
    genes_lesional_up <- deg_sig %>%
      filter(avg_log2FC > 0) %>%
      pull(gene) %>%
      unique()
    
    genes_nonlesional_up <- deg_sig %>%
      filter(avg_log2FC < 0) %>%
      pull(gene) %>%
      unique()
    
    writeLines(
      genes_all,
      con = file.path(
        "Metascape_inputs/By_celltype",
        paste0("Metascape_", prefix, "_all_DEGs.txt")
      )
    )
    
    writeLines(
      genes_lesional_up,
      con = file.path(
        "Metascape_inputs/By_celltype",
        paste0("Metascape_", prefix, "_Lesional_up.txt")
      )
    )
    
    writeLines(
      genes_nonlesional_up,
      con = file.path(
        "Metascape_inputs/By_celltype",
        paste0("Metascape_", prefix, "_Nonlesional_up.txt")
      )
    )
  }
  
  run_enrichment(
    deg_table = deg,
    prefix = prefix,
    outdir = "DEG_GO_KEGG_by_celltype",
    showCategory = 15
  )
}

############################################################
# 14. Global DEG + GO/KEGG
############################################################

Idents(seurat_obj) <- "condition"

deg_global <- FindMarkers(
  seurat_obj,
  ident.1 = "Lesional",
  ident.2 = "Non_lesional",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

deg_global$gene <- rownames(deg_global)

write.csv(
  deg_global,
  file = "Global_DEG_GO_KEGG/Global_DEG_Lesional_vs_Nonlesional.csv",
  row.names = FALSE
)

############################################################
# Metascape input - global DEG
############################################################

deg_global_sig <- deg_global %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

if (nrow(deg_global_sig) > 0) {
  
  global_genes_all <- unique(deg_global_sig$gene)
  
  global_lesional_up <- deg_global_sig %>%
    filter(avg_log2FC > 0) %>%
    pull(gene) %>%
    unique()
  
  global_nonlesional_up <- deg_global_sig %>%
    filter(avg_log2FC < 0) %>%
    pull(gene) %>%
    unique()
  
  writeLines(
    global_genes_all,
    con = "Metascape_inputs/Global/Metascape_Global_all_DEGs.txt"
  )
  
  writeLines(
    global_lesional_up,
    con = "Metascape_inputs/Global/Metascape_Global_Lesional_up.txt"
  )
  
  writeLines(
    global_nonlesional_up,
    con = "Metascape_inputs/Global/Metascape_Global_Nonlesional_up.txt"
  )
}

run_enrichment(
  deg_table = deg_global,
  prefix = "Global_Lesional_vs_Nonlesional",
  outdir = "Global_DEG_GO_KEGG",
  showCategory = 20
)

############################################################
# 15. CellChat
############################################################

seurat_lesional <- subset(seurat_obj, subset = condition == "Lesional")
seurat_nonlesional <- subset(seurat_obj, subset = condition == "Non_lesional")

cat("Lesional hücre sayısı:", ncol(seurat_lesional), "\n")
table(seurat_lesional$HPCA_main)

cat("Non-lesional hücre sayısı:", ncol(seurat_nonlesional), "\n")
table(seurat_nonlesional$HPCA_main)


############################################################
# CellChat - Lesional
############################################################

set.seed(1234)

cellchat_lesional <- createCellChat(
  object = seurat_lesional,
  group.by = "HPCA_main",
  assay = "RNA"
)

CellChatDB <- CellChatDB.human
cellchat_lesional@DB <- CellChatDB

cellchat_lesional <- subsetData(cellchat_lesional)
cellchat_lesional <- identifyOverExpressedGenes(cellchat_lesional)
cellchat_lesional <- identifyOverExpressedInteractions(cellchat_lesional)
cellchat_lesional <- computeCommunProb(cellchat_lesional, type = "triMean")
cellchat_lesional <- filterCommunication(cellchat_lesional, min.cells = 10)

df_lesional <- subsetCommunication(cellchat_lesional)

write.csv(
  df_lesional,
  file = "CellChat_results/CellChat_Lesional_communications.csv",
  row.names = FALSE
)

cellchat_lesional <- computeCommunProbPathway(cellchat_lesional)
cellchat_lesional <- aggregateNet(cellchat_lesional)

groupSize_lesional <- as.numeric(table(cellchat_lesional@idents))

netVisual_circle(
  cellchat_lesional@net$count,
  vertex.weight = groupSize_lesional,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Lesional - Number of interactions"
)

pdf("CellChat_results/CellChat_Lesional_circle_count.pdf", width = 8, height = 8)
netVisual_circle(
  cellchat_lesional@net$count,
  vertex.weight = groupSize_lesional,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Lesional - Number of interactions"
)
dev.off()

cellchat_lesional <- netAnalysis_computeCentrality(
  cellchat_lesional,
  slot.name = "netP"
)

p1 <- netAnalysis_signalingRole_heatmap(
  cellchat_lesional,
  pattern = "outgoing"
)

p2 <- netAnalysis_signalingRole_heatmap(
  cellchat_lesional,
  pattern = "incoming"
)

p1 + p2

pdf(
  "CellChat_results/CellChat_Lesional_outgoing_incoming.pdf",
  width = 18,
  height = 8
)

ComplexHeatmap::draw(p1 + p2)

dev.off()

saveRDS(
  cellchat_lesional,
  file = "CellChat_results/cellchat_lesional.rds"
)


############################################################
# CellChat - Non-lesional
############################################################

set.seed(1234)

cellchat_nonlesional <- createCellChat(
  object = seurat_nonlesional,
  group.by = "HPCA_main",
  assay = "RNA"
)

cellchat_nonlesional@DB <- CellChatDB.human

cellchat_nonlesional <- subsetData(cellchat_nonlesional)
cellchat_nonlesional <- identifyOverExpressedGenes(cellchat_nonlesional)
cellchat_nonlesional <- identifyOverExpressedInteractions(cellchat_nonlesional)
cellchat_nonlesional <- computeCommunProb(cellchat_nonlesional, type = "triMean")
cellchat_nonlesional <- filterCommunication(cellchat_nonlesional, min.cells = 10)

df_nonlesional <- subsetCommunication(cellchat_nonlesional)

write.csv(
  df_nonlesional,
  file = "CellChat_results/CellChat_Nonlesional_communications.csv",
  row.names = FALSE
)

cellchat_nonlesional <- computeCommunProbPathway(cellchat_nonlesional)
cellchat_nonlesional <- aggregateNet(cellchat_nonlesional)

groupSize_nonlesional <- as.numeric(table(cellchat_nonlesional@idents))

netVisual_circle(
  cellchat_nonlesional@net$count,
  vertex.weight = groupSize_nonlesional,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Non-lesional - Number of interactions"
)

pdf("CellChat_results/CellChat_Nonlesional_circle_count.pdf", width = 8, height = 8)
netVisual_circle(
  cellchat_nonlesional@net$count,
  vertex.weight = groupSize_nonlesional,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Non-lesional - Number of interactions"
)
dev.off()

cellchat_nonlesional <- netAnalysis_computeCentrality(
  cellchat_nonlesional,
  slot.name = "netP"
)

p3 <- netAnalysis_signalingRole_heatmap(
  cellchat_nonlesional,
  pattern = "outgoing"
)

p4 <- netAnalysis_signalingRole_heatmap(
  cellchat_nonlesional,
  pattern = "incoming"
)

p3 + p4

pdf(
  "CellChat_results/CellChat_Nonlesional_outgoing_incoming.pdf",
  width = 18,
  height = 8
)

ComplexHeatmap::draw(p3 + p4)

dev.off()

saveRDS(
  cellchat_nonlesional,
  file = "CellChat_results/cellchat_nonlesional.rds"
)


############################################################
# 16. CellChat comparison
############################################################

object.list <- list(
  Non_lesional = cellchat_nonlesional,
  Lesional = cellchat_lesional
)

cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

p_count <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)

p_weight <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

ggsave(
  filename = "CellChat_results/CellChat_compare_count_weight.pdf",
  plot = p_count + p_weight,
  width = 12,
  height = 6
)

netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "count"
)

pdf("CellChat_results/CellChat_difference_count.pdf", width = 8, height = 8)
netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "count"
)
dev.off()

netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight"
)

pdf("CellChat_results/CellChat_difference_weight.pdf", width = 8, height = 8)
netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight"
)
dev.off()

saveRDS(
  cellchat_merged,
  file = "CellChat_results/cellchat_merged_Lesional_Nonlesional.rds"
)


############################################################
# 17. Summary tables and final object
############################################################

cat("\n===== CONDITION =====\n")
table(seurat_obj$condition)

cat("\n===== SAMPLE =====\n")
table(seurat_obj$sample_id)

cat("\n===== CLUSTERS x CONDITION =====\n")
table(seurat_obj$seurat_clusters, seurat_obj$condition)

cat("\n===== HPCA MAIN x CONDITION =====\n")
table(seurat_obj$HPCA_main, seurat_obj$condition)

cat("\n===== HPCA FINE =====\n")
table(seurat_obj$HPCA_fine)

saveRDS(
  seurat_obj,
  file = "GSE147424_AD_final_analysis_object.rds"
)

############################################################
# Metascape input summary
############################################################

metascape_files <- list.files(
  "Metascape_inputs",
  pattern = "\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

metascape_summary <- data.frame(
  file = metascape_files,
  gene_count = sapply(
    metascape_files,
    function(x) length(readLines(x))
  )
)

write.csv(
  metascape_summary,
  file = "Metascape_inputs/Metascape_input_summary.csv",
  row.names = FALSE
)

print(metascape_summary)

message("Analysis completed successfully.")