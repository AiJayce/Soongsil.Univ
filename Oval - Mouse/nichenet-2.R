setwd("/home/Data_Drive_8TB/kykim/Nichenet/")
renv::init()
source("/home/Data_Drive_8TB/kykim/Nichenet/renv/activate.R")
renv::project()

####
save.image("my_session.RData")
save(prior_table_combined, file = "prior_table_combined.RData")
load("/home/Data_Drive_8TB/kykim/Nichenet/prior_table_combined.RData")
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
BiocManager::install("spatstat.core")
devtools::install_github("saeyslab/nichenetr")
reticulate::py_install("anndata", pip=TRUE)
install.packages("Cairo")
install.packages("ragg")
install.packages("NMF")
install.packages("reticulate")
install.packages("pheatmap")
install.packages("gridExtra")
install.packages("ggraph")
install.packages("Seurat")
install.packages("ggplot2", version = "3.4.4")

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
library(RColorBrewer)
library(magrittr)
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
  cells <- colnames(seurat_obj)[seurat_obj@meta.data$module_cluster == celltype]
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
all_ligands <- unique(lr_network$from)
all_receptors <- unique(lr_network$to)
expressed_ligands <- intersect(all_ligands, expressed_genes_sender)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
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
#ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

lr_network_filtered <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
info_tables <- generate_info_tables(
  data,
  celltype_colname = "module_cluster",
  senders_oi = sender_celltypes,
  receivers_oi = receiver,
  lr_network = lr_network_filtered,
  condition_colname = "age_group",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  scenario = "case_control")

prior_table <- generate_prioritization_tables(info_tables$sender_receiver_info,
                                              info_tables$sender_receiver_de,
                                              ligand_activities,
                                              info_tables$lr_condition_de,
                                              scenario = "case_control")
DE_table <- FindAllMarkers( subset(data, subset = age_group == "Aged"),  min.pct = 0,  logfc.threshold = 0,  return.thresh = 1,  features = unique(unlist(lr_network_filtered)))
expression_info <- get_exprs_avg(data, "module_cluster", condition_colname = "age_group", condition_oi = condition_oi, features = unique(unlist(lr_network_filtered)))
condition_markers <- FindMarkers(object = data, ident.1 = condition_oi, ident.2 = condition_reference,group.by = "age_group", min.pct = 0, logfc.threshold = 0,features = unique(unlist(lr_network_filtered))) %>% rownames_to_column("gene")
processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_filtered,senders_oi = sender_celltypes, receivers_oi = receiver)
processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_filtered)
processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_filtered)
prioritizing_weights = c("de_ligand" = 1,
                         "de_receptor" = 1,
                         "activity_scaled" = 1,
                         "exprs_ligand" = 1,
                         "exprs_receptor" = 1,
                         "ligand_condition_specificity" = 1,
                         "receptor_condition_specificity" = 1)
prior_table <- generate_prioritization_tables(processed_expr_table,  processed_DE_table,  ligand_activities,  processed_condition_markers,  prioritizing_weights)

### Filtering Seuratobject < 10 contents
min_cells = 10
count_table <- table(data$module_cluster, data$age_group)
min_counts <- apply(count_table, 1, min)
valid_celltypes <- names(min_counts)[min_counts >= 10]
data <- subset(data, subset = module_cluster %in% valid_celltypes)
table(data$module_cluster, data$age_group)
counts <- table(data@meta.data$module_cluster, data@meta.data$age_group)
keep_ct <- rownames(counts)[apply(counts, 1, function(x) all(x > 0))]
data <- subset(data, subset = module_cluster %in% keep_ct)
expr_mat <- data@assays$RNA@layers$data
if (is.null(rownames(expr_mat))) rownames(expr_mat) <- rownames(data)
colnames(expr_mat) <- colnames(data)
data[["RNA_temp"]] <- Seurat::CreateAssayObject(data = expr_mat)
sender_celltypes <- sender_celltypes[sapply(sender_celltypes, function(ct){
  any(table(data@meta.data$age_group[data@meta.data$module_cluster == ct]) >= min_cells)
})]

receiver_celltypes <- c("module2","module7")
nichenet_outputs <- lapply(c("module2","module7"), function(receiver_ct){
  output <- nichenet_seuratobj_aggregate(
    receiver = receiver_celltypes,
    seurat_obj = data,
    assay = "RNA_temp",  # slot 인자 제거
    condition_colname = "age_group",
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    sender = sender_celltypes,
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    weighted_networks = weighted_networks,
    expression_pct = 0.05
  )
  output$ligand_activities$receiver <- receiver_ct
  return(output)
})

info_tables <- lapply(nichenet_outputs, function(output) {
  lr_network_filtered <- lr_network %>% 
    dplyr::select(from, to) %>% 
    dplyr::filter(from %in% output$ligand_activities$test_ligand &
                    to %in% output$background_expressed_genes) 
  
  generate_info_tables(data, 
                       celltype_colname = "module_cluster", 
                       senders_oi = sender_celltypes, 
                       receivers_oi = unique(output$ligand_activities$receiver), 
                       lr_network_filtered = lr_network_filtered, 
                       condition_colname = "age_group", 
                       condition_oi = condition_oi, 
                       condition_reference = condition_reference, 
                       scenario = "case_control") 
})
info_tables_combined <- purrr::pmap(info_tables, bind_rows)
ligand_activities_combined <- purrr::map_dfr(nichenet_outputs, "ligand_activities")

prior_table_combined <- generate_prioritization_tables(
  sender_receiver_info = info_tables_combined$sender_receiver_info %>% distinct,
  sender_receiver_de = info_tables_combined$sender_receiver_de,
  ligand_activities = ligand_activities_combined,
  lr_condition_de = info_tables_combined$lr_condition_de %>% distinct,
  scenario = "case_control")
prior_table_oi <- prior_table_combined %>% slice_max(prioritization_score, n = 50)

#### binding ####
all_celltypes <- unique(data@meta.data$module_cluster)
sender <- setdiff(all_celltypes, c("module2","module7"))
senders_receivers <- prior_table_oi %>%
  dplyr::select(sender, receiver) %>%
  lapply(as.character) %>%  # factor -> character
  unlist() %>%
  unique() %>%
  sort()


n <- length(senders_receivers)
pal <- rep(RColorBrewer::brewer.pal(12, 'Set3'), length.out = n)
celltype_colors <- pal %>% set_names(senders_receivers)
pdf("/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/circos_plot.pdf", width = 20, height = 15, bg = "white")
print(make_circos_lr(prior_table_oi,
                     colors_sender = celltype_colors,
                     colors_receiver = celltype_colors))
dev.off()



receiver_oi <- c("module2","module7")
legend_adjust <- c(0.5, 0.5)
legend_adjust <- c(1, 1)
legend_adjust <- c("right", "top")

pdf("/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/Mush_plot.pdf",
    width = 5, height = 5, bg = "white")

my_make_mushroom_plot(
  prior_table_combined %>% filter(receiver == receiver_oi),
  top_n = 5, 
  true_color_range = TRUE,
  show_rankings = TRUE,
  show_all_datapoints = TRUE
) +
  theme(
    text = element_text(size = 6),
    
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    
    strip.text = element_text(size = 6),
    
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    
    axis.title.x = element_text(
      hjust = 0.25,
      margin = margin(b = 15)
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0,
      vjust = 0,
      margin = margin(b = -1)
    )
  )

dev.off()

#Mushroomplot Defind (Custom fuction)
{
#Mushroomplot fuction Define
my_make_mushroom_plot <- function(prioritization_table, top_n = 30, show_rankings = FALSE, 
                                    show_all_datapoints = FALSE, true_color_range = TRUE, use_absolute_rank = FALSE, 
                                    size = "scaled_avg_exprs", color = "scaled_p_val_adapted", 
                                    ligand_fill_colors = c("#DEEBF7", "#08306B"), 
                                    receptor_fill_colors = c("#FEE0D2", "#A50F15"), 
                                    unranked_ligand_fill_colors = c(scales::alpha("#FFFFFF", alpha = 0.2), 
                                                                    scales::alpha("#252525", alpha = 0.2)), 
                                    unranked_receptor_fill_colors = c(scales::alpha("#FFFFFF", alpha = 0.2), 
                                                                      scales::alpha("#252525", alpha = 0.2)), ...) 
  {
    size_ext <- c("ligand", "receptor")
    color_ext <- c("ligand", "receptor")
    if (size == "pct_expressed") 
      size_ext <- c("sender", "receiver")
    if (color == "pct_expressed") 
      color_ext <- c("sender", "receiver")
    cols_to_use <- c("sender", "ligand", "receptor", paste0(size, "_", size_ext), 
                     paste0(color, "_", color_ext))
    
    if (!all(cols_to_use %in% colnames(prioritization_table))) {
      stop(paste(paste0("`", cols_to_use %>% .[!. %in% colnames(prioritization_table)], 
                        "`", collapse = ", "), "column not in prioritization table"))
    }
    if (!is.logical(show_rankings) | length(show_rankings) != 1) 
      stop("show_rankings should be a TRUE or FALSE")
    if (!is.logical(show_all_datapoints) | length(show_all_datapoints) != 1) 
      stop("show_all_datapoints should be a TRUE or FALSE")
    if (!is.logical(true_color_range) | length(true_color_range) != 1) 
      stop("true_color_range should be a TRUE or FALSE")
    if (!is.logical(use_absolute_rank) | length(use_absolute_rank) != 1) 
      stop("use_absolute_rank should be a TRUE or FALSE")
    if (!is.numeric(top_n) | length(top_n) != 1) 
      stop("top_n should be a numeric vector of length 1")
    if (length(ligand_fill_colors) != 2) 
      stop("ligand_fill_colors should be a vector of length 2")
    if (length(receptor_fill_colors) != 2) 
      stop("receptor_fill_colors should be a vector of length 2")
    if (length(unranked_ligand_fill_colors) != 2) 
      stop("unranked_ligand_fill_colors should be a vector of length 2")
    if (length(unranked_receptor_fill_colors) != 2) 
      stop("unranked_receptor_fill_colors should be a vector of length 2")
    
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    requireNamespace("ggnewscale")
    requireNamespace("ggforce")
    requireNamespace("shadowtext")
    requireNamespace("cowplot")
    
    # ★★★ 핵심 수정: desc() 앞에 dplyr:: 추가 ★★★
    if (!"prioritization_rank" %in% colnames(prioritization_table)) {
      prioritization_table <- prioritization_table %>% 
        dplyr::mutate(prioritization_rank = rank(dplyr::desc(prioritization_score)))
    }
    prioritization_table <- prioritization_table %>% 
      dplyr::mutate(relative_rank = rank(dplyr::desc(prioritization_score)))
    
    rank_filter_col <- ifelse(use_absolute_rank, "prioritization_rank", "relative_rank")
    
    filtered_table <- prioritization_table %>% 
      dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = " - "))
    
    order_interactions <- unique(filtered_table %>% 
                                   dplyr::filter(.data[[rank_filter_col]] <= top_n) %>% 
                                   dplyr::pull(lr_interaction))
    
    filtered_table <- filtered_table %>% 
      dplyr::filter(lr_interaction %in% order_interactions) %>% 
      dplyr::mutate(lr_interaction = factor(lr_interaction, levels = rev(order_interactions)))
    
    if (nrow(filtered_table) == 0) {
      stop("No ligand-receptor interactions found in the top_n. Please try use_absolute_rank = FALSE or increase top_n.")
    }
    
    if (!is.factor(filtered_table$sender)) {
      filtered_table$sender <- as.factor(filtered_table$sender)
    } else {
      filtered_table$sender <- droplevels(filtered_table$sender)
    }
    
    lr_interaction_vec <- 1:length(order_interactions) %>% setNames(order_interactions)
    
    filtered_table <- filtered_table %>% 
      dplyr::select(c("lr_interaction", all_of(cols_to_use), "prioritization_rank", "relative_rank")) %>% 
      tidyr::pivot_longer(c(ligand, receptor), names_to = "type", values_to = "protein") %>% 
      dplyr::mutate(
        size = ifelse(type == "ligand", 
                      get(paste0(size, "_", size_ext[1])), 
                      get(paste0(size, "_", size_ext[2]))),
        color = ifelse(type == "ligand", 
                       get(paste0(color, "_", color_ext[1])), 
                       get(paste0(color, "_", color_ext[2])))
      ) %>% 
      dplyr::select(-contains(c("_ligand", "_receptor", "_sender", "_receiver"))) %>% 
      dplyr::mutate(start = rep(c(-pi, 0), nrow(filtered_table))) %>% 
      dplyr::mutate(x = as.numeric(sender), y = lr_interaction_vec[lr_interaction])
    
    if (any(filtered_table$size < 0) | any(filtered_table$size > 1.001)) {
      stop("Size column is not scaled between 0 and 1. Please use this column as the color instead.")
    }
    
    keywords_adj <- c("LFC", "pval", "", "product", "mean", "adjusted", "expression") %>% 
      setNames(c("lfc", "p", "val", "prod", "avg", "adj", "exprs"))
    
    size_title <- sapply(stringr::str_split(size, "_")[[1]], 
                         function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>% 
      paste0(., collapse = " ") %>% 
      stringr::str_replace("^\\w{1}", toupper)
    
    color_title <- sapply(stringr::str_split(color, "_")[[1]], 
                          function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>% 
      paste0(., collapse = " ") %>% 
      stringr::str_replace("^\\w{1}", toupper)
    
    color_lims <- c(0, 1)
    if (true_color_range) 
      color_lims <- NULL
    
    scale <- 0.5
    ncelltypes <- length(unique(filtered_table$sender))
    n_interactions <- length(lr_interaction_vec)
    
    legend2_df <- data.frame(
      values = c(0.25, 0.5, 0.75, 1), 
      x = (ncelltypes + 2.5):(ncelltypes + 5.5), 
      y = rep(floor(n_interactions/3), 4), 
      start = -pi
    )
    
    axis_rect <- data.frame(
      xmin = 0, xmax = ncelltypes + 1, 
      ymin = 0, ymax = n_interactions + 1
    )
    
    panel_grid_y <- data.frame(
      x = rep(seq(from = 0.5, to = ncelltypes + 0.5, by = 1), each = 2), 
      y = c(n_interactions + 1, 0), 
      group = rep(1:(ncelltypes + 1), each = 2)
    )
    
    panel_grid_x <- data.frame(
      y = rep(seq(from = 0.5, to = n_interactions + 0.5, by = 1), each = 2), 
      x = c(ncelltypes + 1, 0), 
      group = rep(1:(n_interactions + 1), each = 2)
    )
    
    theme_args <- list(
      panel.grid.major = element_blank(), 
      legend.box = "horizontal", 
      panel.background = element_blank()
    )
    theme_args[names(list(...))] <- list(...)
    
    scale_legend_title_size <- 3.57
    if ("legend.title" %in% names(theme_args)) {
      if (!is.null(theme_args$legend.title$size)) {
        scale_legend_title_size <- theme_args$legend.title$size * (5/14)
      }
    }
    
    scale_legend_text_size <- 3.57
    if ("legend.text" %in% names(theme_args)) {
      if (!is.null(theme_args$legend.text$size)) {
        scale_legend_text_size <- theme_args$legend.text$size * (5/14)
      }
    }
    
    if (!"legend.justification" %in% names(theme_args)) {
      theme_args$legend.justification <- c(1, 0.7)
    }
    if (!"legend.position" %in% names(theme_args)) {
      theme_args$legend.position <- c(1, 0.7)
    }
    
    p1 <- ggplot2::ggplot() + 
      ggforce::geom_arc_bar(
        data = filtered_table %>% dplyr::filter(type == "ligand", .data[[rank_filter_col]] <= top_n), 
        ggplot2::aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size) * scale, 
                     start = start, end = start + pi, fill = color), 
        color = "white"
      ) + 
      ggplot2::scale_fill_gradient(
        low = ligand_fill_colors[1], 
        high = ligand_fill_colors[2], 
        limits = color_lims, 
        oob = scales::squish, 
        n.breaks = 3, 
        guide = ggplot2::guide_colorbar(order = 1), 
        name = paste0(color_title, " (", color_ext[1], ")") %>% stringr::str_wrap(width = 15)
      ) + 
      ggnewscale::new_scale_fill() + 
      ggforce::geom_arc_bar(
        data = filtered_table %>% dplyr::filter(type == "receptor", .data[[rank_filter_col]] <= top_n), 
        ggplot2::aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size) * scale, 
                     start = start, end = start + pi, fill = color), 
        color = "white"
      ) + 
      ggforce::geom_arc_bar(
        data = legend2_df, 
        ggplot2::aes(x0 = x, y0 = y, r0 = 0, r = sqrt(values) * scale, 
                     start = start, end = start + pi), 
        fill = "black"
      ) + 
      ggplot2::geom_rect(
        data = legend2_df, 
        ggplot2::aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5), 
        color = "gray90", fill = NA
      ) + 
      ggplot2::geom_text(
        data = legend2_df, 
        ggplot2::aes(label = values, x = x, y = y - 0.6), 
        vjust = 1, size = scale_legend_text_size
      ) + 
      ggplot2::geom_text(
        data = data.frame(
          x = (ncelltypes + 4), 
          y = floor(n_interactions/3) + 1, 
          label = size_title %>% stringr::str_wrap(width = 15)
        ), 
        ggplot2::aes(x = x, y = y, label = label), 
        size = scale_legend_title_size, 
        vjust = 0, lineheight = 0.75
      ) + 
      ggplot2::geom_line(
        data = panel_grid_y, 
        ggplot2::aes(x = x, y = y, group = group), 
        color = "gray90"
      ) + 
      ggplot2::geom_line(
        data = panel_grid_x, 
        ggplot2::aes(x = x, y = y, group = group), 
        color = "gray90"
      ) + 
      ggplot2::geom_rect(
        data = axis_rect, 
        ggplot2::aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), 
        color = "black", fill = "transparent"
      ) + 
      ggplot2::scale_fill_gradient(
        low = receptor_fill_colors[1], 
        high = receptor_fill_colors[2], 
        limits = color_lims, 
        oob = scales::squish, 
        n.breaks = 3, 
        name = paste0(color_title, " (", color_ext[2], ")") %>% stringr::str_wrap(width = 15), 
        guide = ggplot2::guide_colorbar(order = 2)
      ) + 
      ggplot2::scale_y_continuous(
        breaks = n_interactions:1, 
        labels = names(lr_interaction_vec), 
        expand = ggplot2::expansion(add = c(0, 0))
      ) + 
      ggplot2::scale_x_continuous(
        breaks = 1:ncelltypes, 
        labels = levels(filtered_table$sender), 
        position = "top", 
        expand = ggplot2::expansion(add = c(0, 0))
      ) + 
      ggplot2::xlab("Sender cell types") + 
      ggplot2::ylab("Ligand-receptor interaction") + 
      ggplot2::coord_fixed() + 
      do.call(ggplot2::theme, theme_args)
    
    if (show_all_datapoints) {
      unranked_ligand_lims <- c(0, 1)
      unranked_receptor_lims <- c(0, 1)
      if (true_color_range) {
        unranked_ligand_lims <- filtered_table %>% 
          dplyr::filter(type == "ligand", .data[[rank_filter_col]] <= top_n) %>% 
          dplyr::select(color) %>% range
        unranked_receptor_lims <- filtered_table %>% 
          dplyr::filter(type == "receptor", .data[[rank_filter_col]] <= top_n) %>% 
          dplyr::select(color) %>% range
      }
      p1 <- p1 + 
        ggnewscale::new_scale_fill() + 
        ggforce::geom_arc_bar(
          data = filtered_table %>% dplyr::filter(type == "ligand", .data[[rank_filter_col]] > top_n), 
          ggplot2::aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size) * scale, 
                       start = start, end = start + pi, fill = color), 
          color = "white"
        ) + 
        ggplot2::scale_fill_gradient(
          low = unranked_ligand_fill_colors[1], 
          high = unranked_ligand_fill_colors[2], 
          limits = unranked_ligand_lims, 
          oob = scales::squish, 
          guide = "none"
        ) + 
        ggnewscale::new_scale_fill() + 
        ggforce::geom_arc_bar(
          data = filtered_table %>% dplyr::filter(type == "receptor", .data[[rank_filter_col]] > top_n), 
          ggplot2::aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size) * scale, 
                       start = start, end = start + pi, fill = color), 
          color = "white"
        ) + 
        ggplot2::scale_fill_gradient(
          low = unranked_receptor_fill_colors[1], 
          high = unranked_receptor_fill_colors[2], 
          limits = unranked_receptor_lims, 
          oob = scales::squish, 
          guide = "none"
        )
    }
    
    if (show_rankings) {
      p1 <- p1 + 
        shadowtext::geom_shadowtext(
          data = filtered_table %>% dplyr::filter(.data[[rank_filter_col]] <= top_n), 
          ggplot2::aes(x = x, y = y, label = prioritization_rank)
        )
    }
    
    p1
  }

}
