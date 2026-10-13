# IGF1R
# Stomach Cancer Analysis Code
# ============================================================
# 0. 环境与路径
# ============================================================
PROJECT_DIR <- "<PROJECT_DIR>"   # 替换为实际项目根目录
setwd(PROJECT_DIR)

library(Seurat)
library(SeuratObject)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(enrichR)
library(harmony)
library(limma)
library(infercnv)

set.seed(123)

out_dir <- file.path(PROJECT_DIR, "results_IGF1R_USP30")
fig_dir <- file.path(PROJECT_DIR, "figures_IGF1R_USP30")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. 数据读取与 Seurat 对象构建
# ============================================================
csv_files <- list.files("GSE183904_RAW", pattern = "\\.csv\\.gz$", full.names = TRUE)
seurat_list <- list()

for (f in csv_files) {
  sample_name <- gsub("\\.csv\\.gz$", "", basename(f))
  df <- fread(f, header = TRUE, showProgress = FALSE)
  gene_names <- df[[1]]
  df[, 1 := NULL]
  mat <- as.matrix(df)
  rownames(mat) <- gene_names
  rm(df, gene_names); gc()
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  obj <- RenameCells(obj, add.cell.id = sample_name)
  seurat_list[[sample_name]] <- obj
  rm(mat, obj); gc()
}

merged_seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)
rm(seurat_list); gc()

# ============================================================
# 2. 质控、SCTransform、Harmony 整合与聚类
# ============================================================
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
merged_seurat <- subset(
  merged_seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25
)

obj_list <- SplitObject(merged_seurat, split.by = "orig.ident")
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- SCTransform(
    obj_list[[i]],
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
}

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
all_samples <- merge(obj_list[[1]], y = obj_list[-1])

all_samples <- RunPCA(all_samples, features = features, npcs = 50)
DefaultAssay(all_samples) <- "SCT"
all_samples <- RunHarmony(
  all_samples,
  group.by.vars = "orig.ident",
  reduction = "pca",
  reduction.save = "harmony"
)
all_samples <- RunUMAP(all_samples, reduction = "harmony", dims = 1:30)
all_samples <- FindNeighbors(all_samples, reduction = "harmony", dims = 1:30)
all_samples <- FindClusters(all_samples, resolution = 0.5)

# ============================================================
# 3. 细胞类型手动注释
# ============================================================
marker_list <- list(
  Epithelial   = c("EPCAM", "KRT8", "KRT18", "CDH1"),
  Fibroblast   = c("COL1A1", "DCN", "LUM", "COL1A2"),
  T_cell       = c("CD3D", "CD3E", "CD2"),
  Myeloid      = c("CD68", "CD14", "FCGR3A", "ITGAM"),
  B_cell       = c("CD79A", "MS4A1", "CD19"),
  Endothelial  = c("PECAM1", "VWF", "CLDN5"),
  Mast         = c("TPSAB1", "CPA3"),
  Plasma       = c("JCHAIN", "MZB1", "SDC1")
)

DefaultAssay(all_samples) <- "SCT"
avg_expr <- AverageExpression(
  all_samples,
  assays = "SCT",
  group.by = "seurat_clusters",
  slot = "data"
)$SCT

cluster_ids <- colnames(avg_expr)
annotation <- setNames(rep("Unknown", length(cluster_ids)), cluster_ids)

for (cl in cluster_ids) {
  scores <- sapply(names(marker_list), function(type) {
    genes <- intersect(marker_list[[type]], rownames(avg_expr))
    if (length(genes) == 0) return(-Inf)
    mean(avg_expr[genes, cl])
  })
  annotation[cl] <- names(which.max(scores))
}

names(annotation) <- gsub("^g", "", names(annotation))
cell_annot <- annotation[as.character(all_samples$seurat_clusters)]
names(cell_annot) <- colnames(all_samples)
all_samples <- AddMetaData(all_samples, metadata = cell_annot, col.name = "manual_celltype")

# ============================================================
# 4. inferCNV 恶性上皮鉴定
# ============================================================
epi_cells <- subset(all_samples, manual_celltype == "Epithelial")
ref_cells <- subset(all_samples, manual_celltype == "Endothelial")

epi_counts <- GetAssayData(epi_cells, assay = "RNA", layer = "counts")
ref_counts <- GetAssayData(ref_cells, assay = "RNA", layer = "counts")

common_genes <- intersect(rownames(epi_counts), rownames(ref_counts))
combined_counts <- cbind(epi_counts[common_genes, ], ref_counts[common_genes, ])

cnv_cells <- CreateSeuratObject(counts = combined_counts, assay = "RNA")
cnv_cells$cell_group <- ifelse(
  colnames(cnv_cells) %in% colnames(epi_cells),
  "Epithelial", "Normal"
)

annotations <- data.frame(
  row.names = colnames(cnv_cells),
  cell_group = cnv_cells$cell_group
)
write.table(annotations, "inferCNV_annotations.txt",
            sep = "\t", quote = FALSE, col.names = FALSE)

gene_order_raw <- read.table("gencode_v21_gen_pos.complete.txt",
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_order_raw) <- c("gene_info", "chr", "start", "end")
gene_order_raw$symbol <- gsub("\\|.*", "", gene_order_raw$gene_info)

my_gene <- rownames(combined_counts)
matched_idx <- match(my_gene, gene_order_raw$symbol)
custom_order <- data.frame(
  gene = my_gene[!is.na(matched_idx)],
  chr  = gene_order_raw$chr[matched_idx[!is.na(matched_idx)]],
  start = gene_order_raw$start[matched_idx[!is.na(matched_idx)]],
  end   = gene_order_raw$end[matched_idx[!is.na(matched_idx)]]
)
write.table(custom_order, "custom_gene_order.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = combined_counts,
  annotations_file = "inferCNV_annotations.txt",
  gene_order_file = "custom_gene_order.txt",
  ref_group_names = "Normal"
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "./inferCNV_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 8,
  output_format = "pdf"
)

malignant_cells <- colnames(epi_cells)
all_samples$cell_group <- ifelse(
  colnames(all_samples) %in% malignant_cells,
  "Malignant_Epithelial", "Other_Cells"
)

# ============================================================
# 5. IGF1R 表达与恶性上皮亚群分析
# ============================================================
# 恶性上皮亚群
DefaultAssay(epi_cells) <- "RNA"
epi_cells <- NormalizeData(epi_cells)
epi_cells <- FindVariableFeatures(epi_cells, selection.method = "vst", nfeatures = 2000)
epi_cells <- ScaleData(epi_cells)
epi_cells <- RunPCA(epi_cells, features = VariableFeatures(epi_cells), npcs = 30)
epi_cells <- RunUMAP(epi_cells, dims = 1:15)
epi_cells <- FindNeighbors(epi_cells, dims = 1:15)
epi_cells <- FindClusters(epi_cells, resolution = 0.4)

# 功能评分
epi_cells <- CellCycleScoring(
  epi_cells,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

stem_genes <- c("PROM1", "CD44", "ALDH1A1", "SOX2", "NANOG", "OCT4",
                "KLF4", "LGR5", "CD133", "ALCAM", "CD24")
stem_genes <- intersect(stem_genes, rownames(epi_cells))
epi_cells <- AddModuleScore(epi_cells, features = list(Stem = stem_genes),
                            name = "Stemness")

emt_genes <- c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1",
               "ZEB1", "ZEB2", "CDH1", "CTNNB1", "TGFB1", "MMP2", "MMP9")
emt_genes <- intersect(emt_genes, rownames(epi_cells))
epi_cells <- AddModuleScore(epi_cells, features = list(EMT = emt_genes),
                            name = "EMT")

epi_cells$Stemness_Score <- epi_cells$Stemness1
epi_cells$EMT_Score <- epi_cells$EMT1

# 保存
saveRDS(all_samples, "all_samples_final.rds")
saveRDS(epi_cells, "epi_cells_final.rds")

# ============================================================
# 6. 免疫检查点分析：PD-L1、PDCD1、CD163、T cell exhaustion score
# ============================================================
all_samples <- readRDS("all_samples_final.rds")
epi_cells <- readRDS("epi_cells_final.rds")

# 重建细胞分组
malignant_cells <- colnames(epi_cells)
all_samples$cell_group <- ifelse(
  colnames(all_samples) %in% malignant_cells,
  "Malignant_Epithelial", "Other_Cells"
)

# 按患者计算恶性上皮 IGF1R 平均表达，划分 IGF1R_high / IGF1R_low
patient_igf1r <- FetchData(all_samples, vars = c("orig.ident", "IGF1R", "cell_group")) %>%
  filter(cell_group == "Malignant_Epithelial") %>%
  group_by(orig.ident) %>%
  summarise(mean_IGF1R = mean(IGF1R), .groups = "drop")

median_igf1r <- median(patient_igf1r$mean_IGF1R)
patient_igf1r$patient_IGF1R_group <- ifelse(
  patient_igf1r$mean_IGF1R > median_igf1r,
  "IGF1R_high", "IGF1R_low"
)

all_samples$patient_IGF1R_group <- patient_igf1r$patient_IGF1R_group[
  match(all_samples$orig.ident, patient_igf1r$orig.ident)
]

# ---------- PD-L1 (CD274) 在恶性上皮 ----------
malig_cells <- subset(all_samples, cell_group == "Malignant_Epithelial")
p_pdl1 <- VlnPlot(
  malig_cells, features = "CD274", group.by = "patient_IGF1R_group",
  cols = c("#E41A1C", "#377EB8"), pt.size = 0.05
) +
  labs(title = "PD-L1 (CD274) in malignant epithelial cells",
       y = "Log-norm expression", x = "") +
  stat_compare_means(
    comparisons = list(c("IGF1R_high", "IGF1R_low")),
    method = "wilcox.test", label = "p.signif"
  ) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "PDL1_malignant_IGF1R_high_vs_low.pdf"),
       p_pdl1, width = 6, height = 5)

# ---------- PDCD1 在 T 细胞 ----------
t_cells <- subset(all_samples, manual_celltype == "T_cell")
p_pdcd1 <- VlnPlot(
  t_cells, features = "PDCD1", group.by = "patient_IGF1R_group",
  cols = c("#E41A1C", "#377EB8"), pt.size = 0.05
) +
  labs(title = "PDCD1 expression in T cells",
       y = "Log-norm expression", x = "") +
  stat_compare_means(
    comparisons = list(c("IGF1R_high", "IGF1R_low")),
    method = "wilcox.test", label = "p.signif"
  ) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "PDCD1_Tcell_IGF1R_high_vs_low.pdf"),
       p_pdcd1, width = 6, height = 5)

# ---------- CD163 在髓系细胞 ----------
myeloid_cells <- subset(all_samples, manual_celltype == "Myeloid")
p_cd163 <- VlnPlot(
  myeloid_cells, features = "CD163", group.by = "patient_IGF1R_group",
  cols = c("#E41A1C", "#377EB8"), pt.size = 0.05
) +
  labs(title = "CD163 in myeloid cells",
       y = "Log-norm expression", x = "") +
  stat_compare_means(
    comparisons = list(c("IGF1R_high", "IGF1R_low")),
    method = "wilcox.test", label = "p.signif"
  ) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "CD163_myeloid_IGF1R_high_vs_low.pdf"),
       p_cd163, width = 6, height = 5)

# ---------- T cell exhaustion score ----------
exhaustion_genes <- list(c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX"))
t_cells <- AddModuleScore(t_cells, features = exhaustion_genes, name = "Exhaustion_Score")

df_exhaust <- FetchData(t_cells, vars = c("Exhaustion_Score1", "patient_IGF1R_group"))
colnames(df_exhaust) <- c("score", "group")

p_exhaust <- ggplot(df_exhaust, aes(x = group, y = score, fill = group)) +
  geom_violin(scale = "width", adjust = 1.5, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.08, fill = "white", alpha = 0.7, outlier.shape = NA, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("IGF1R_high" = "#E41A1C", "IGF1R_low" = "#377EB8")) +
  stat_compare_means(
    comparisons = list(c("IGF1R_high", "IGF1R_low")),
    method = "wilcox.test", label = "p.signif"
  ) +
  labs(title = "T cell exhaustion score", x = "", y = "Module score") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(fig_dir, "Tcell_exhaustion_score_IGF1R_high_vs_low.pdf"),
       p_exhaust, width = 6, height = 6)

# ============================================================
# 7. IGF1R × USP30 患者级分析
# ============================================================
# 以下代码整合自患者级分析脚本，路径已替换，删除 Bulk 验证部分

# ---------- 7.1 读取对象与准备 ----------
epi_cells <- readRDS("epi_cells_final.rds")
all_samples <- readRDS("all_samples_final.rds")

epi_cells$patient_id <- epi_cells$orig.ident
epi_cells$group_c1   <- epi_cells$cluster1_vs_rest

# 取原始 RNA counts
get_counts <- function(obj, assay = "RNA") {
  a <- obj[[assay]]
  if (inherits(a, "Assay5")) {
    LayerData(a, layer = "counts")
  } else {
    GetAssayData(obj, assay = assay, slot = "counts")
  }
}
epi_cells <- JoinLayers(epi_cells, assay = "RNA")
cnt <- get_counts(epi_cells, assay = "RNA")

# 患者 × 亚群 pseudobulk
df_expr <- data.frame(
  cell       = colnames(cnt),
  patient_id = epi_cells$patient_id,
  group_c1   = epi_cells$group_c1,
  IGF1R      = as.numeric(cnt["IGF1R", ]),
  USP30      = as.numeric(cnt["USP30", ]),
  stringsAsFactors = FALSE
)

pseudobulk <- df_expr %>%
  group_by(patient_id, group_c1) %>%
  summarise(
    IGF1R   = sum(IGF1R),
    USP30   = sum(USP30),
    n_cells = n(),
    .groups = "drop"
  )

pb_all_raw <- df_expr %>%
  group_by(patient_id) %>%
  summarise(
    IGF1R   = sum(IGF1R),
    USP30   = sum(USP30),
    n_cells = n(),
    .groups = "drop"
  )

pb_c1_raw    <- pseudobulk %>% filter(group_c1 == "Cluster1")
pb_other_raw <- pseudobulk %>% filter(group_c1 == "Other_Malignant")

# Spearman + bootstrap CI
spearman_ci <- function(x, y, B = 10000, seed = 123) {
  set.seed(seed)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3) {
    return(list(rho = NA, ci = c(NA, NA), p = NA, n = n))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  boot <- replicate(B, {
    idx <- sample(seq_len(n), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  list(
    rho = unname(ct$estimate),
    ci  = as.numeric(quantile(boot, c(0.025, 0.975), na.rm = TRUE)),
    p   = ct$p.value,
    n   = n
  )
}

run_spearman <- function(df, label) {
  r <- spearman_ci(df$IGF1R, df$USP30)
  data.frame(
    analysis     = label,
    n_patients   = r$n,
    rho          = r$rho,
    ci_low       = r$ci[1],
    ci_high      = r$ci[2],
    p            = r$p,
    median_cells = median(df$n_cells),
    min_cells    = min(df$n_cells),
    max_cells    = max(df$n_cells)
  )
}

res_raw <- rbind(
  run_spearman(pb_all_raw,   "Malignant_epi_all_rawcounts"),
  run_spearman(pb_c1_raw,    "Malignant_epi_Cluster1_rawcounts"),
  run_spearman(pb_other_raw, "Malignant_epi_Other_rawcounts")
)

# 主分析：nCount_RNA CPM 归一化
meta_epi <- data.frame(
  cell       = colnames(epi_cells),
  patient_id = epi_cells$patient_id,
  group_c1   = epi_cells$group_c1,
  stringsAsFactors = FALSE
)
meta_epi$total_nCount <- epi_cells@meta.data$nCount_RNA[
  match(meta_epi$cell, colnames(epi_cells))
]

pb_norm <- meta_epi %>%
  mutate(
    IGF1R = as.numeric(cnt["IGF1R", ]),
    USP30 = as.numeric(cnt["USP30", ])
  ) %>%
  group_by(patient_id, group_c1) %>%
  summarise(
    IGF1R        = sum(IGF1R),
    USP30        = sum(USP30),
    total_nCount = sum(total_nCount),
    n_cells      = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    IGF1R_cpm = IGF1R / total_nCount * 1e6,
    USP30_cpm = USP30 / total_nCount * 1e6
  )

pb_all_norm <- meta_epi %>%
  mutate(
    IGF1R = as.numeric(cnt["IGF1R", ]),
    USP30 = as.numeric(cnt["USP30", ])
  ) %>%
  group_by(patient_id) %>%
  summarise(
    IGF1R        = sum(IGF1R),
    USP30        = sum(USP30),
    total_nCount = sum(total_nCount),
    n_cells      = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    IGF1R_cpm = IGF1R / total_nCount * 1e6,
    USP30_cpm = USP30 / total_nCount * 1e6
  )

pb_c1_norm    <- pb_norm %>% filter(group_c1 == "Cluster1")
pb_other_norm <- pb_norm %>% filter(group_c1 == "Other_Malignant")

run_spearman_norm <- function(df, label) {
  r <- spearman_ci(df$IGF1R_cpm, df$USP30_cpm)
  data.frame(
    analysis     = label,
    n_patients   = r$n,
    rho          = r$rho,
    ci_low       = r$ci[1],
    ci_high      = r$ci[2],
    p            = r$p,
    median_cells = median(df$n_cells),
    min_cells    = min(df$n_cells),
    max_cells    = max(df$n_cells)
  )
}

res_norm <- rbind(
  run_spearman_norm(pb_all_norm,   "Malignant_epi_all_CPM"),
  run_spearman_norm(pb_c1_norm,    "Malignant_epi_Cluster1_CPM"),
  run_spearman_norm(pb_other_norm, "Malignant_epi_Other_CPM")
)

# 敏感性分析：最小细胞数阈值
run_spearman_thresh <- function(pb, label, min_n, value_col = c("CPM", "raw")) {
  value_col <- match.arg(value_col)
  df <- pb %>% filter(n_cells >= min_n)
  if (nrow(df) < 3) {
    return(data.frame(
      analysis = label, min_cells = min_n, n_patients = nrow(df),
      rho = NA, ci_low = NA, ci_high = NA, p = NA
    ))
  }
  if (value_col == "CPM") {
    r <- spearman_ci(df$IGF1R_cpm, df$USP30_cpm)
  } else {
    r <- spearman_ci(df$IGF1R, df$USP30)
  }
  data.frame(
    analysis   = label,
    min_cells  = min_n,
    n_patients = r$n,
    rho        = r$rho,
    ci_low     = r$ci[1],
    ci_high    = r$ci[2],
    p          = r$p
  )
}

thresh_tab <- rbind(
  run_spearman_thresh(pb_all_norm,   "all",      0,  "CPM"),
  run_spearman_thresh(pb_all_norm,   "all",      10, "CPM"),
  run_spearman_thresh(pb_all_norm,   "all",      20, "CPM"),
  run_spearman_thresh(pb_c1_norm,    "Cluster1", 0,  "CPM"),
  run_spearman_thresh(pb_c1_norm,    "Cluster1", 5,  "CPM"),
  run_spearman_thresh(pb_c1_norm,    "Cluster1", 10, "CPM"),
  run_spearman_thresh(pb_c1_norm,    "Cluster1", 20, "CPM"),
  run_spearman_thresh(pb_other_norm, "Other",    0,  "CPM"),
  run_spearman_thresh(pb_other_norm, "Other",    10, "CPM"),
  run_spearman_thresh(pb_other_norm, "Other",    20, "CPM")
)

# Leave-one-patient-out
loo_spearman <- function(df, label, value_col = c("CPM", "raw")) {
  value_col <- match.arg(value_col)
  if (value_col == "CPM") {
    x <- df$IGF1R_cpm; y <- df$USP30_cpm
  } else {
    x <- df$IGF1R; y <- df$USP30
  }
  if (length(x) < 4) {
    return(data.frame(
      analysis = label, n = length(x),
      rho_full = NA, rho_loo_min = NA, rho_loo_max = NA
    ))
  }
  rho_full <- suppressWarnings(cor(x, y, method = "spearman"))
  rho_loo <- sapply(seq_along(x), function(i) {
    suppressWarnings(cor(x[-i], y[-i], method = "spearman"))
  })
  data.frame(
    analysis    = label,
    n           = length(x),
    rho_full    = rho_full,
    rho_loo_min = min(rho_loo, na.rm = TRUE),
    rho_loo_max = max(rho_loo, na.rm = TRUE)
  )
}

loo_tab <- rbind(
  loo_spearman(pb_all_norm,   "all_CPM",      "CPM"),
  loo_spearman(pb_c1_norm,    "Cluster1_CPM", "CPM"),
  loo_spearman(pb_other_norm, "Other_CPM",    "CPM")
)

# 图 1：CPM 归一化患者级 pseudobulk 散点
plot_norm <- function(df, title, size_by_cells = FALSE) {
  r <- spearman_ci(df$IGF1R_cpm, df$USP30_cpm)
  sub <- sprintf(
    "Spearman rho = %.3f, 95%% CI [%.3f, %.3f], P = %.3g, n = %d",
    r$rho, r$ci[1], r$ci[2], r$p, r$n
  )
  p <- ggplot(df, aes(x = IGF1R_cpm, y = USP30_cpm))
  if (size_by_cells) {
    p <- p +
      geom_point(aes(size = n_cells), alpha = 0.8, color = "#2C7FB8") +
      scale_size_continuous(range = c(1.5, 6), name = "n cells")
  } else {
    p <- p + geom_point(size = 2.5, alpha = 0.85, color = "#2C7FB8")
  }
  p +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
    geom_vline(xintercept = median(df$IGF1R_cpm), linetype = "dashed",
               color = "grey60") +
    geom_hline(yintercept = median(df$USP30_cpm), linetype = "dashed",
               color = "grey60") +
    labs(
      title = title,
      subtitle = sub,
      x = "IGF1R (CPM, patient-level pseudobulk)",
      y = "USP30 (CPM, patient-level pseudobulk)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.subtitle = element_text(size = 9.5),
      plot.title = element_text(size = 12, face = "bold")
    )
}

p_all_norm   <- plot_norm(pb_all_norm,   "Malignant epithelial cells (all)")
p_c1_norm    <- plot_norm(pb_c1_norm,    "Malignant epithelial Cluster 1",
                          size_by_cells = TRUE)
p_other_norm <- plot_norm(pb_other_norm, "Other malignant epithelial cells")

ggsave(file.path(fig_dir, "Fig1a_malignant_all_CPM_scatter.pdf"),
       p_all_norm, width = 5.5, height = 5)
ggsave(file.path(fig_dir, "Fig1b_Cluster1_CPM_scatter.pdf"),
       p_c1_norm, width = 5.5, height = 5)
ggsave(file.path(fig_dir, "Fig1c_OtherMalignant_CPM_scatter.pdf"),
       p_other_norm, width = 5.5, height = 5)

p_main_norm <- (p_all_norm | p_c1_norm | p_other_norm) +
  plot_annotation(
    tag_levels = "A",
    title = "Patient-level IGF1R–USP30 correlation (CPM-normalized pseudobulk)"
  )
ggsave(file.path(fig_dir, "Fig1_main_CPM_panel.pdf"),
       p_main_norm, width = 15, height = 5.2)

# 补图 S1：min_cells 敏感性
thresh_tab$analysis <- factor(
  thresh_tab$analysis,
  levels = c("Cluster1", "all", "Other")
)
p_sens <- ggplot(
  thresh_tab,
  aes(x = factor(min_cells), y = rho, color = analysis, group = analysis)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.5) +
  geom_point(size = 2.8) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Cluster1" = "#D73027",
      "all"      = "#4575B4",
      "Other"    = "#91BFDB"
    ),
    labels = c(
      "Cluster1" = "Cluster 1",
      "all"      = "All malignant",
      "Other"    = "Other malignant"
    )
  ) +
  labs(
    x = "Minimum cells per patient",
    y = "Spearman rho (95% CI)",
    color = NULL,
    title = "Sensitivity to minimum-cell threshold",
    subtitle = "Patient-level pseudobulk, CPM-normalized"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
ggsave(file.path(fig_dir, "FigS1_min_cells_sensitivity.pdf"),
       p_sens, width = 6.5, height = 4.5)

# 补图 S2：LOO 稳定性
loo_dist <- function(df, label, value_col = "CPM") {
  if (value_col == "CPM") {
    x <- df$IGF1R_cpm; y <- df$USP30_cpm
  } else {
    x <- df$IGF1R; y <- df$USP30
  }
  rho_full <- suppressWarnings(cor(x, y, method = "spearman"))
  rho_loo <- sapply(seq_along(x), function(i) {
    suppressWarnings(cor(x[-i], y[-i], method = "spearman"))
  })
  data.frame(
    analysis = label,
    rho_loo  = rho_loo,
    rho_full = rho_full
  )
}

loo_df <- rbind(
  loo_dist(pb_all_norm,   "All malignant"),
  loo_dist(pb_c1_norm,    "Cluster 1"),
  loo_dist(pb_other_norm, "Other malignant")
)
loo_df$analysis <- factor(
  loo_df$analysis,
  levels = c("Cluster 1", "All malignant", "Other malignant")
)
p_loo <- ggplot(loo_df, aes(x = analysis, y = rho_loo, fill = analysis)) +
  geom_boxplot(width = 0.5, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.7) +
  geom_hline(data = loo_df, aes(yintercept = rho_full),
             linetype = "dashed", color = "grey30") +
  scale_fill_manual(
    values = c(
      "Cluster 1"       = "#D73027",
      "All malignant"   = "#4575B4",
      "Other malignant" = "#91BFDB"
    )
  ) +
  labs(
    x = NULL,
    y = "Spearman rho (leave-one-patient-out)",
    title = "Leave-one-patient-out stability",
    subtitle = "Dashed line = full-data rho"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
ggsave(file.path(fig_dir, "FigS2_LOO_stability.pdf"),
       p_loo, width = 6, height = 4.5)

# ---------- 7.2 各细胞类型 IGF1R / USP30 表达与阳性率 ----------
all_samples$celltype_plot <- as.character(all_samples$manual_celltype)
all_samples$celltype_plot[
  all_samples$manual_celltype == "Epithelial"
] <- "Malignant_Epithelial"

c1_cells <- colnames(epi_cells)[epi_cells$cluster1_vs_rest == "Cluster1"]
all_samples$celltype_plot[colnames(all_samples) %in% c1_cells] <-
  "Malignant_Cluster1"

DefaultAssay(all_samples) <- "RNA"
all_samples <- JoinLayers(all_samples, assay = "RNA")

expr_mat <- LayerData(all_samples[["RNA"]], layer = "data")
rna_counts <- LayerData(all_samples[["RNA"]], layer = "counts")

expr_df <- data.frame(
  cell     = colnames(all_samples),
  celltype = all_samples$celltype_plot,
  IGF1R    = as.numeric(expr_mat["IGF1R", ]),
  USP30    = as.numeric(expr_mat["USP30", ]),
  stringsAsFactors = FALSE
)

cnt_df <- data.frame(
  cell      = colnames(all_samples),
  celltype  = all_samples$celltype_plot,
  IGF1R_pos = as.numeric(rna_counts["IGF1R", ] > 0),
  USP30_pos = as.numeric(rna_counts["USP30", ] > 0),
  stringsAsFactors = FALSE
)

summary_df <- expr_df %>%
  group_by(celltype) %>%
  summarise(
    n_cells    = n(),
    mean_IGF1R = mean(IGF1R),
    mean_USP30 = mean(USP30),
    .groups = "drop"
  )

pos_df <- cnt_df %>%
  group_by(celltype) %>%
  summarise(
    pct_IGF1R_RNA = mean(IGF1R_pos) * 100,
    pct_USP30_RNA = mean(USP30_pos) * 100,
    both_pos_RNA  = mean(IGF1R_pos > 0 & USP30_pos > 0) * 100,
    .groups = "drop"
  )

comb2 <- merge(summary_df, pos_df, by = "celltype")

# 图 2a：IGF1R 细胞类型 DotPlot
ct_order <- c(
  "Malignant_Cluster1", "Malignant_Epithelial",
  "Fibroblast", "Endothelial", "Myeloid",
  "T_cell", "B_cell", "Plasma", "Mast"
)
ct_order <- intersect(ct_order, comb2$celltype)

dot_igf1r <- data.frame(
  celltype  = comb2$celltype,
  mean_expr = comb2$mean_IGF1R,
  pct_pos   = comb2$pct_IGF1R_RNA
) %>%
  filter(celltype %in% ct_order)
dot_igf1r$celltype <- factor(dot_igf1r$celltype, levels = ct_order)

p2a <- ggplot(dot_igf1r, aes(x = "IGF1R", y = celltype)) +
  geom_point(aes(size = pct_pos, color = mean_expr)) +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(2, 10)) +
  labs(
    x = NULL, y = NULL,
    size = "% positive", color = "Mean expr",
    title = "IGF1R expression across cell types",
    subtitle = "Cluster 1 shows highest IGF1R positivity"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "Fig2a_IGF1R_celltype_dotplot.pdf"),
       p2a, width = 4.5, height = 4.5)

# Cluster1 内 IGF1R+ vs IGF1R− 的 USP30
c1_df <- meta_epi %>% filter(group_c1 == "Cluster1")
c1_df$IGF1R_pos <- as.numeric(cnt["IGF1R", c1_df$cell] > 0)
c1_df$USP30_pos <- as.numeric(cnt["USP30", c1_df$cell] > 0)

binom_ci <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  bt <- binom.test(x, n)$conf.int
  bt * 100
}

bar_df <- c1_df %>%
  group_by(IGF1R_pos) %>%
  summarise(
    n     = n(),
    n_pos = sum(USP30_pos),
    pct   = mean(USP30_pos) * 100,
    .groups = "drop"
  ) %>%
  mutate(group = ifelse(IGF1R_pos == 1, "IGF1R+", "IGF1R-"))

bar_df$ci_low  <- mapply(function(x, n) binom_ci(x, n)[1],
                         bar_df$n_pos, bar_df$n)
bar_df$ci_high <- mapply(function(x, n) binom_ci(x, n)[2],
                         bar_df$n_pos, bar_df$n)

p2b <- ggplot(bar_df, aes(x = group, y = pct, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
    vjust = -0.4, size = 3.5
  ) +
  scale_fill_manual(values = c("IGF1R-" = "#91BFDB", "IGF1R+" = "#D73027")) +
  labs(
    x = NULL,
    y = "USP30 positive cells (%)",
    title = "USP30 positivity in Cluster 1",
    subtitle = "Patient-level single-cell counts; Fisher test in text"
  ) +
  ylim(0, max(bar_df$ci_high, na.rm = TRUE) * 1.3) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
ggsave(file.path(fig_dir, "Fig2b_USP30_positivity_by_IGF1R.pdf"),
       p2b, width = 4.5, height = 4.5)

p_fig2 <- (p2a | p2b) + plot_annotation(tag_levels = "A")
ggsave(file.path(fig_dir, "Fig2_main_panel.pdf"),
       p_fig2, width = 9.5, height = 4.8)

# ---------- 7.3 组织来源标签与患者级 pseudobulk 按细胞类型 ----------
primary_tumor <- c(
  "GSM5573467_sample2", "GSM5573468_sample3", "GSM5573470_sample5",
  "GSM5573472_sample7", "GSM5573473_sample8", "GSM5573475_sample10",
  "GSM5573477_sample12", "GSM5573478_sample13", "GSM5573479_sample14",
  "GSM5573480_sample15", "GSM5573481_sample16", "GSM5573482_sample17",
  "GSM5573483_sample18", "GSM5573487_sample22", "GSM5573489_sample24",
  "GSM5573491_sample26", "GSM5573492_sample27", "GSM5573493_sample28",
  "GSM5573494_sample29", "GSM5573495_sample30", "GSM5573497_sample32",
  "GSM5573498_sample33", "GSM5573499_sample34", "GSM5573501_sample36",
  "GSM5573504_sample39", "GSM5573505_sample40"
)
primary_normal <- c(
  "GSM5573466_sample1", "GSM5573469_sample4", "GSM5573471_sample6",
  "GSM5573474_sample9", "GSM5573476_sample11", "GSM5573486_sample21",
  "GSM5573488_sample23", "GSM5573490_sample25", "GSM5573496_sample31",
  "GSM5573500_sample35"
)
peritoneum_samples <- c(
  "GSM5573484_sample19", "GSM5573485_sample20",
  "GSM5573503_sample38", "GSM5573502_sample37"
)

all_samples$tissue_source <- ifelse(
  all_samples$orig.ident %in% primary_tumor, "Primary_Tumor",
  ifelse(
    all_samples$orig.ident %in% primary_normal, "Primary_Normal",
    ifelse(
      all_samples$orig.ident %in% peritoneum_samples, "Peritoneum", "Other"
    )
  )
)

cell_types <- unique(all_samples$manual_celltype)
expr_data  <- GetAssayData(all_samples, assay = "RNA", layer = "data")

expression_summary <- data.frame()
for (ct in cell_types) {
  cells_ct <- colnames(all_samples)[all_samples$manual_celltype == ct]
  if (length(cells_ct) < 50) next
  igf1r_vals <- expr_data["IGF1R", cells_ct]
  usp30_vals <- expr_data["USP30", cells_ct]
  expression_summary <- rbind(
    expression_summary,
    data.frame(
      cell_type      = ct,
      n_cells        = length(cells_ct),
      IGF1R_mean     = round(mean(igf1r_vals), 4),
      IGF1R_pos_rate = round(mean(igf1r_vals > 0) * 100, 2),
      USP30_mean     = round(mean(usp30_vals), 4),
      USP30_pos_rate = round(mean(usp30_vals > 0) * 100, 2),
      both_pos_rate  = round(mean(igf1r_vals > 0 & usp30_vals > 0) * 100, 2)
    )
  )
}
expression_summary <- expression_summary %>% arrange(desc(both_pos_rate))
write.csv(expression_summary,
          file.path(out_dir, "IGF1R_USP30_celltype_summary.csv"),
          row.names = FALSE)

plot_data <- expression_summary %>%
  select(cell_type, IGF1R_pos_rate, USP30_pos_rate, both_pos_rate) %>%
  pivot_longer(-cell_type, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(
      metric,
      levels = c("IGF1R_pos_rate", "USP30_pos_rate", "both_pos_rate")
    )
  )
p_bar <- ggplot(
  plot_data,
  aes(x = reorder(cell_type, value), y = value, fill = metric)
) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "IGF1R_pos_rate" = "#E41A1C",
      "USP30_pos_rate" = "#377EB8",
      "both_pos_rate"  = "#4DAF4A"
    ),
    labels = c(
      "IGF1R_pos_rate" = "IGF1R+",
      "USP30_pos_rate" = "USP30+",
      "both_pos_rate"  = "Both+"
    ),
    name = ""
  ) +
  labs(
    title = "Positivity rate by cell type",
    x = "", y = "Positivity (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
ggsave(file.path(fig_dir, "IGF1R_USP30_positivity_barplot.pdf"),
       p_bar, width = 9, height = 5)

# 患者级 pseudobulk 按细胞类型：仅原发肿瘤
tumor_samples <- all_samples[, all_samples$tissue_source == "Primary_Tumor"]
celltypes_to_test <- expression_summary$cell_type[
  expression_summary$USP30_pos_rate > 3 &
    expression_summary$n_cells > 500
]

pseudobulk_list <- list()
for (ct in celltypes_to_test) {
  cells_ct <- colnames(tumor_samples)[tumor_samples$manual_celltype == ct]
  if (length(cells_ct) < 100) next
  sub_obj <- tumor_samples[, cells_ct]
  sub_obj$orig.ident <- gsub("_", "-", sub_obj$orig.ident)
  pb <- AggregateExpression(
    sub_obj,
    group.by      = "orig.ident",
    assays        = "RNA",
    slot          = "counts",
    return.seurat = FALSE
  )$RNA
  if (!all(c("IGF1R", "USP30") %in% rownames(pb))) next
  n_cells_df <- data.frame(
    patient = names(table(sub_obj$orig.ident)),
    n_cells = as.numeric(table(sub_obj$orig.ident))
  )
  total_counts <- Matrix::colSums(pb)
  pb_df <- data.frame(
    patient    = colnames(pb),
    cell_type  = ct,
    IGF1R_raw  = as.numeric(pb["IGF1R", ]),
    USP30_raw  = as.numeric(pb["USP30", ]),
    total      = as.numeric(total_counts)
  ) %>%
    left_join(n_cells_df, by = "patient") %>%
    mutate(
      IGF1R_cpm    = IGF1R_raw / total * 1e6,
      USP30_cpm    = USP30_raw / total * 1e6,
      IGF1R_logcpm = log1p(IGF1R_cpm),
      USP30_logcpm = log1p(USP30_cpm)
    )
  pseudobulk_list[[ct]] <- pb_df
}
pseudobulk_all <- bind_rows(pseudobulk_list)

cor_results <- data.frame()
for (ct in unique(pseudobulk_all$cell_type)) {
  df_ct <- pseudobulk_all %>% filter(cell_type == ct, n_cells >= 10)
  if (nrow(df_ct) < 5) next
  test <- cor.test(df_ct$IGF1R_logcpm, df_ct$USP30_logcpm, method = "spearman")
  n       <- nrow(df_ct)
  z       <- atanh(test$estimate)
  se      <- 1 / sqrt(n - 3)
  ci_low  <- tanh(z - 1.96 * se)
  ci_high <- tanh(z + 1.96 * se)
  cor_results <- rbind(
    cor_results,
    data.frame(
      cell_type  = ct,
      n_patients = n,
      rho        = round(test$estimate, 3),
      CI_low     = round(ci_low, 3),
      CI_high    = round(ci_high, 3),
      p_value    = test$p.value
    )
  )
}
if (nrow(cor_results) > 0) {
  cor_results$p_adj <- p.adjust(cor_results$p_value, method = "BH")
  cor_results <- cor_results %>% arrange(p_value)
  write.csv(cor_results,
            file.path(out_dir, "IGF1R_USP30_patient_pseudobulk_correlation_by_celltype.csv"),
            row.names = FALSE)
}

if (nrow(cor_results) > 0) {
  cor_results$cell_type <- factor(
    cor_results$cell_type,
    levels = c("Epithelial", "Fibroblast", "Endothelial")
  )
  cor_results$label <- sprintf(
    "rho = %.2f\nFDR = %.1e",
    cor_results$rho,
    cor_results$p_adj
  )
  p_forest <- ggplot(cor_results, aes(x = rho, y = cell_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(
      aes(xmin = CI_low, xmax = CI_high),
      height = 0.15, color = "gray40", linewidth = 0.8
    ) +
    geom_point(size = 4, color = "#E41A1C") +
    geom_text(
      aes(label = label),
      nudge_y = 0.35,
      size = 3.5, lineheight = 1.0
    ) +
    scale_x_continuous(
      limits = c(0, 1.15),
      breaks = seq(0, 1, 0.25)
    ) +
    labs(
      title = "IGF1R-USP30 correlation by cell type (patient-level pseudobulk)",
      subtitle = "GSE183904 primary tumor samples",
      x = "Spearman rho (95% CI)",
      y = ""
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  ggsave(file.path(fig_dir, "IGF1R_USP30_forest_plot.pdf"),
         p_forest, width = 9, height = 4.5)
}

# ---------- 7.4 全样本患者级 pseudobulk ----------
pb_allcells <- AggregateExpression(
  tumor_samples,
  group.by      = "orig.ident",
  assays        = "RNA",
  slot          = "counts",
  return.seurat = FALSE
)$RNA
total_allcells <- Matrix::colSums(pb_allcells)
pb_allcells_df <- data.frame(
  patient   = colnames(pb_allcells),
  IGF1R_raw = as.numeric(pb_allcells["IGF1R", ]),
  USP30_raw = as.numeric(pb_allcells["USP30", ]),
  total     = as.numeric(total_allcells)
) %>%
  mutate(
    IGF1R_cpm    = IGF1R_raw / total * 1e6,
    USP30_cpm    = USP30_raw / total * 1e6,
    IGF1R_logcpm = log1p(IGF1R_cpm),
    USP30_logcpm = log1p(USP30_cpm)
  )

test_allcells <- cor.test(
  pb_allcells_df$IGF1R_logcpm,
  pb_allcells_df$USP30_logcpm,
  method = "spearman"
)

p_allcells <- ggplot(pb_allcells_df, aes(x = IGF1R_logcpm, y = USP30_logcpm)) +
  geom_point(alpha = 0.7, size = 3, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, fill = "#3498db", alpha = 0.2) +
  labs(
    title = "All cells pseudobulk (primary tumor, n = 26 patients)",
    subtitle = sprintf(
      "Spearman rho = %.3f, P = %.2e",
      test_allcells$estimate,
      test_allcells$p.value
    ),
    x = "IGF1R (logCPM, patient pseudobulk)",
    y = "USP30 (logCPM, patient pseudobulk)"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(fig_dir, "Scatter_IGF1R_USP30_AllCells_primaryTumor.pdf"),
       p_allcells, width = 5, height = 5)

# 肿瘤 vs 正常患者级 pseudobulk
tissue_plot_data <- data.frame()
for (tissue in c("Primary_Tumor", "Primary_Normal")) {
  samples_tissue <- all_samples[, all_samples$tissue_source == tissue]
  if (ncol(samples_tissue) < 50) next
  pb_tissue <- AggregateExpression(
    samples_tissue,
    group.by      = "orig.ident",
    assays        = "RNA",
    slot          = "counts",
    return.seurat = FALSE
  )$RNA
  total_tissue <- Matrix::colSums(pb_tissue)
  tissue_plot_data <- rbind(
    tissue_plot_data,
    data.frame(
      tissue    = tissue,
      patient   = colnames(pb_tissue),
      IGF1R_raw = as.numeric(pb_tissue["IGF1R", ]),
      USP30_raw = as.numeric(pb_tissue["USP30", ]),
      total     = as.numeric(total_tissue)
    )
  )
}
tissue_plot_data <- tissue_plot_data %>%
  mutate(
    IGF1R_logcpm = log1p(IGF1R_raw / total * 1e6),
    USP30_logcpm = log1p(USP30_raw / total * 1e6)
  )
tissue_plot_data$tissue <- factor(
  tissue_plot_data$tissue,
  levels = c("Primary_Normal", "Primary_Tumor"),
  labels = c("Normal (n=10)", "Tumor (n=26)")
)
p_tissue <- ggplot(
  tissue_plot_data,
  aes(x = IGF1R_logcpm, y = USP30_logcpm, color = tissue)
) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(
    values = c("Normal (n=10)" = "#4DBBD5", "Tumor (n=26)" = "#E64B35")
  ) +
  labs(
    title = "IGF1R-USP30 correlation: Tumor vs Normal",
    x = "IGF1R (logCPM, patient pseudobulk)",
    y = "USP30 (logCPM, patient pseudobulk)",
    color = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )
ggsave(file.path(fig_dir, "Scatter_IGF1R_USP30_Tumor_vs_Normal.pdf"),
       p_tissue, width = 6, height = 5)

# ---------- 7.5 Cluster1 补充细胞级分析 ----------
c1_obj <- subset(epi_cells, cluster1_vs_rest == "Cluster1")
df_c1_cell <- FetchData(c1_obj, vars = c("IGF1R", "USP30"))
test_c1_cell <- cor.test(
  df_c1_cell$IGF1R,
  df_c1_cell$USP30,
  method = "spearman"
)
p_c1_cell <- ggplot(df_c1_cell, aes(x = IGF1R, y = USP30)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#2c3e50") +
  geom_smooth(method = "lm", se = TRUE, fill = "#3498db", alpha = 0.2) +
  labs(
    title = "IGF1R vs USP30 in Cluster 1 (malignant epithelial)",
    subtitle = sprintf(
      "Spearman rho = %.3f, P = %.3g",
      test_c1_cell$estimate,
      test_c1_cell$p.value
    ),
    x = "IGF1R expression",
    y = "USP30 expression"
  ) +
  theme_classic(base_size = 14)
ggsave(file.path(fig_dir, "IGF1R_USP30_Cluster1_scatter.pdf"),
       p_c1_cell, width = 6, height = 5)

df2 <- FetchData(epi_cells, vars = c("USP30", "group_c1"))
colnames(df2) <- c("USP30", "Group")
p_vln <- ggplot(df2, aes(x = Group, y = USP30, fill = Group)) +
  geom_violin(scale = "width", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(
    values = c("Cluster1" = "#E41A1C", "Other_Malignant" = "#999999")
  ) +
  labs(
    title = "USP30 expression: Cluster1 vs Other Malignant",
    y = "USP30 log-normalized expression",
    x = ""
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "USP30_Cluster1_vs_Others.pdf"),
       p_vln, width = 6, height = 5)

df_pair <- FetchData(epi_cells, vars = c("USP30", "group_c1", "orig.ident"))
paired <- df_pair %>%
  group_by(orig.ident, group_c1) %>%
  summarise(
    mean_USP30 = mean(USP30),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group_c1,
    values_from = c(mean_USP30, n)
  ) %>%
  filter(n_Cluster1 >= 10, n_Other_Malignant >= 10)
test_paired <- wilcox.test(
  paired$mean_USP30_Cluster1,
  paired$mean_USP30_Other_Malignant,
  paired = TRUE
)

# ============================================================
# 8. 完成
# ============================================================
cat("\n分析完成。\n")
cat("输出目录：", out_dir, "\n")
cat("图目录：", fig_dir, "\n")
