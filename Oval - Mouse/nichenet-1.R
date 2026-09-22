setwd("/home/Data_Drive_8TB/kykim/Nichenet/")
renv::init()
source("/home/Data_Drive_8TB/kykim/Nichenet/renv/activate.R")
renv::project()

####
save.image("my_session.RData")
save(cellchat, file = "cellchat_object.RData")
load("my_session.RData")
####

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")}
BiocManager::install(c("BiocGenerics", "BiocNeighbors", "ComplexHeatmap"))
BiocManager::install("Biobase")
BiocManager::install("biomaRt")
BiocManager::install("org.Mm.eg.db", ask = FALSE)
BiocManager::install("limma")
remotes::install_github("saeyslab/nichenetr@v2.0.4")
remotes::install_version("Seurat", version = "4.0.2")
reticulate::py_install("anndata", pip=TRUE)

install.packages("Cairo")
install.packages("ragg")
install.packages("NMF")
install.packages("reticulate")
install.packages("pheatmap")
install.packages("gridExtra")
install.packages("ggraph")
install.packages("Seurat")

library(patchwork)
library(dplyr)
library(Cairo)
library(tibble)
library(tidyr)
library(ragg)
library(reticulate)
library(Seurat)
library(nichenetr)
library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(igraph)
library(ggraph)
library(limma)
library(nichenetr)
library(dplyr)
library(stringr)


ad <- import("anndata")
data <- ad$read_h5ad("/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/CD8_Module.h5ad")
cell_names <- py_to_r(data$obs_names$tolist())
gene_names <- py_to_r(data$var_names$tolist())
expr_matrix <- t(as.matrix(data$X))
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- cell_names
data <- CreateSeuratObject(counts = expr_matrix, meta.data = data$obs)

data@meta.data$age_numeric <- str_extract(data@meta.data$age, "\\d+")
data@meta.data$age_numeric <- as.numeric(data@meta.data$age_numeric)
data@meta.data$age_group <- ifelse(data@meta.data$age_numeric > 20, "Aged", "Young")
data@meta.data$age_group <- factor(data@meta.data$age_group, levels = c("Young", "Aged"))
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)

## Reference data import
organism <- "mouse"
base_dir <- "/home/Data_Drive_8TB/kykim/Nichenet/nichenet_files"
lr_network <- readRDS(file.path(base_dir, "lr_network_mouse_21122021.rds"))
gr_network <- readRDS(file.path(base_dir, "gr_network_mouse_21122021.rds"))
signaling_network <- readRDS(file.path(base_dir, "signaling_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(base_dir, "ligand_target_matrix_nsga2r_final_mouse.rds"))
ligand_tf_matrix <- readRDS(file.path(base_dir, "ligand_tf_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(file.path(base_dir, "weighted_networks_nsga2r_final_mouse.rds"))

lr_network <- lr_network %>% distinct(from, to)

Idents(data) <- data@meta.data$module_cluster
receiver = c("module2","module7")
expr_mat_receiver <- GetAssayData(data, assay = "RNA", slot = "data")
rm_receiver <- Matrix::rowMeans(expr_mat_receiver)
expressed_genes_receiver <- names(rm_receiver[rm_receiver > 0.1])

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
all_clusters <- levels(Idents(data))
sender_celltypes <- setdiff(all_clusters, receiver)

get_expressed_genes <- function(celltype, seurat_obj, pct_threshold) {
  cells <- colnames(seurat_obj)[seurat_obj@meta.data$module_cluster== celltype]
  if(length(cells) == 0) {
    warning(paste("No cells found for", celltype))
    return(character(0))
  }
  expr_matrix <- LayerData(seurat_obj, layer = "data")[, cells, drop = FALSE]
  if(!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if(ncol(expr_matrix) == 1) {
    pct_expressed <- as.numeric(expr_matrix > 0)
    names(pct_expressed) <- rownames(expr_matrix)
  } else {
    pct_expressed <- rowSums(expr_matrix > 0) / ncol(expr_matrix)
  }
  genes <- names(pct_expressed[pct_expressed >= pct_threshold])
  return(genes)
}

list_expressed_genes_sender <- sender_celltypes %>% unique() %>%  lapply(function(ct) get_expressed_genes(ct, data, 0.05))
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

condition_oi <-  "Aged"
condition_reference <- "Young"
seurat_obj_receiver <- subset(data, subset = module_cluster == receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "age_group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

list_expressed_genes_receiver <- receiver %>% lapply(function(ct) get_expressed_genes(ct, data, 0.05))
expressed_genes_receiver <- list_expressed_genes_receiver %>% 
  unlist() %>% 
  unique()

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>%arrange(-aupr_corrected) %>%mutate(rank = rank(-aupr_corrected))
#ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/histo_lig_acticity.png", type = "png", width = 2000, height = 1000, units = "px", dpi = 235)
print(p_hist_lig_activity)
dev.off()

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
vis_ligand_aupr <- ligand_activities %>% 
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% 
  dplyr::select(aupr_corrected) %>%  # dplyr:: 추가
  arrange(aupr_corrected) %>% 
  as.matrix(ncol = 1)

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/AUPR.png", type = "png", width = 500, height = 1000, units = "px", dpi = 180)
(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  
dev.off()


## L-R result Network ##
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 50) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.60) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/Prioritized_ligands.png", type = "png", width = 3000, height = 1000, units = "px", dpi = 150)
make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
dev.off()


## Prior Prioritized_Potential_ligands interaction Network ##
get_weighted_ligand_receptor_links <- function(ligands, receptors, lr_network, lr_sig) {
  lr_net <- as.data.frame(lr_network)
  names(lr_net)[1:2] <- c("ligand", "receptor")
  
  ligand_receptor_links <- lr_net %>% 
    dplyr::filter(!is.na(ligand) & !is.na(receptor)) %>%
    dplyr::filter(ligand %in% ligands & receptor %in% receptors)
  
  if(nrow(ligand_receptor_links) == 0) {
    return(data.frame(ligand = character(), receptor = character(), weight = numeric()))
  }
  
  lr_sig_df <- as.data.frame(lr_sig)
  if(ncol(lr_sig_df) >= 3) {
    names(lr_sig_df)[1:3] <- c("ligand", "receptor", "weight")
    ligand_receptor_links <- ligand_receptor_links %>% dplyr::left_join(lr_sig_df[,1:3], by = c("ligand", "receptor"))
  }
  
  if(!"weight" %in% colnames(ligand_receptor_links)) {
    ligand_receptor_links$weight <- 1
  } else {
    ligand_receptor_links$weight[is.na(ligand_receptor_links$weight)] <- 1
  }
  
  return(ligand_receptor_links)
}

prepare_ligand_receptor_visualization <- function(ligand_receptor_links_df, ligands, order_hclust = "both") {
  vis_matrix <- ligand_receptor_links_df %>% 
    dplyr::filter(!is.na(ligand) & !is.na(receptor)) %>%
    dplyr::select(ligand, receptor, weight) %>% 
    tidyr::pivot_wider(names_from = receptor, values_from = weight, values_fill = 0) %>% 
    dplyr::filter(!is.na(ligand)) %>%
    tibble::column_to_rownames("ligand") %>% 
    as.matrix()
  
  if(any(is.na(rownames(vis_matrix)))) {
    vis_matrix <- vis_matrix[!is.na(rownames(vis_matrix)), , drop = FALSE]
  }
  if(any(is.na(colnames(vis_matrix)))) {
    vis_matrix <- vis_matrix[, !is.na(colnames(vis_matrix)), drop = FALSE]
  }
  
  valid_ligands <- intersect(ligands[!is.na(ligands)], rownames(vis_matrix))
  vis_matrix <- vis_matrix[valid_ligands, , drop = FALSE]
  
  if((order_hclust == "both" || order_hclust == "ligand") && nrow(vis_matrix) > 1) vis_matrix <- vis_matrix[hclust(dist(vis_matrix))$order, , drop = FALSE]
  if((order_hclust == "both" || order_hclust == "receptor") && ncol(vis_matrix) > 1) vis_matrix <- vis_matrix[, hclust(dist(t(vis_matrix)))$order, drop = FALSE]
  return(vis_matrix)
}
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df, best_upstream_ligands, order_hclust = "both") 

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/Ligand_Receptor.png", type = "png", width = 2500, height = 1500, units = "px", dpi = 225)
p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), y_name = "Receptors", x_name = "Ligands", color = "mediumvioletred", legend_title = "Prior interaction potential") + scale_fill_gradient2(low = "whitesmoke", high = "mediumvioletred")
print(p_ligand_receptor)
dev.off()

Idents(data) <- "module_cluster"
p_dotplot <- DotPlot(subset(data, idents = sender_celltypes), features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() + 
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.05, vjust = 1), panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), plot.margin = margin(10, 10, 30, 10))

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/Ligand_expression.png", type = "png", width = 2000, height = 1500, units = "px", dpi = 180)
print(p_dotplot)
dev.off()

