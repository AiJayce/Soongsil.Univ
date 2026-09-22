library(readr)


ligands_oi <- c('Ntm')
targets_oi <- c("Apoe","Trem2")

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 1.5,
                                                      minmax_scaling = TRUE) 


graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

Cairo(file = "/home/Data_Drive_8TB/kykim/6. Alzheimer/Nichnet/network_LR(Apoe).png", type = "png", width = 1500, height = 2000, units = "px", dpi = 1000)
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))
dev.off()

################### csv file output ####################
output_path <- "/home/Data_Drive_8TB/kykim/6. Alzheimer/Nichnet/"

edges_csv <- bind_rows( active_signaling_network$sig %>% dplyr::mutate(layer = "signaling"), active_signaling_network$gr %>% dplyr::mutate(layer = "regulatory")) %>% dplyr::select(from, to, weight, layer)
write_csv(edges_csv,  paste0(output_path, "active_signaling_network_together.csv"))

################### LE LR pair rank ####################

sender_cells <- subset(data, module_cluster == "CD8_T_cells")
receiver_cells <- subset(data, module_cluster == "module7")
expr_matrix_sender <- GetAssayData(sender_cells, assay = "RNA", slot = "data")
expr_matrix_receiver <- GetAssayData(receiver_cells, assay = "RNA", slot = "data")
expressed_ligands <- rownames(expr_matrix_sender)[Matrix::rowMeans(expr_matrix_sender > 0) > 0.1]
expressed_ligands <- intersect(expressed_ligands, unique(weighted_networks$lr_sig$from))
expressed_receptors <- rownames(expr_matrix_receiver)[Matrix::rowMeans(expr_matrix_receiver > 0) > 0.1]
expressed_receptors <- intersect(expressed_receptors, unique(weighted_networks$lr_sig$to))
ligand_oi <- c(
  "Inhba","Il24","Tdgf1","Gdf2","Col5a1","Gdf1","Hp","Fgf4","Il1b","Hpx",
  "Fgf1","Tnf","Mstn","Bmp4","Fgf10","Col3a1","Tgfb1","Bmp8a","Gdf7","Egf",
  "Col5a3","Eng","Orm1","Col5a2","Bmp1","Ccl21b","Ocln","Bmp7","Tgfb3","Csf3")
lr_pairs <- weighted_networks$lr_sig %>%
  dplyr::filter(from %in% expressed_ligands,to %in% expressed_receptors,from %in% ligand_oi ) %>%  dplyr::arrange(dplyr::desc(weight)) %>%  dplyr::slice_head(n = 5000)
lr_pairs <- lr_pairs %>%  rowwise() %>%  mutate(
    ligand_expr   = mean(expr_matrix_sender[from, ]),
    receptor_expr = mean(expr_matrix_receiver[to, ]),
    activity      = ligand_expr * receptor_expr) %>%  ungroup()
targets <- c("Tgfbr1","Smad4","Jun","Smad1","Trp53","Sp1","Stat3")
lr_pairs %>%
  dplyr::arrange(dplyr::desc(weight)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::filter(to %in% targets) %>%
  dplyr::select(from, to, rank) %>%
  print()

targets_df <- tibble::tribble(~from,  ~to, ~rank,
  "Tgfb1","Tgfbr1", 2,
  "Tgfb1","Smad4", 7,
  "Tgfb1","Trp53", 21,
  "Tgfb1","Stat3", 27,
  "Tgfb1","Smad1", 29,
  "Tgfb1","Sp1", 55,
  "Tgfb1","Jun", 56)

plot_df <- lr_pairs %>%
  dplyr::inner_join(targets_df, by = c("from","to")) %>%
  dplyr::mutate(to_label = paste0(to, " (", rank, ")")) %>%  dplyr::arrange(dplyr::desc(rank)) %>%  dplyr::mutate(to_label = factor(to_label, levels = to_label))

Cairo(file = "/home/Data_Drive_8TB/kykim/1. Oval/Nichenet/CD_8T_ligand.png",type = "png", width = 2000, height = 2500, units = "px", dpi = 350)
ggplot(plot_df, aes(x = from, y = to_label)) +
  geom_point(aes(size = activity, color = weight)) +
  scale_size_continuous(range = c(3, 12)) +
  scale_color_viridis_c(option = "D") +
  labs(x = "Ligand (CD8 T cell)", y = "Receptor (Luminal Epithelium)",size = "Activity",color = "LR weight" ) + theme_minimal(base_size = 16)
dev.off()







nodes <- unique(c(ligands_oi, targets_oi, unlist(active_signaling_network)))
nodes_type <- ifelse(nodes %in% ligands_oi, "ligand",
                     ifelse(nodes %in% targets_oi, "target",
                            ifelse(nodes %in% colnames(ligand_tf_matrix), "TF",
                                   ifelse(nodes %in% weighted_networks$lr_sig$to, "receptor",
                                          "unknown"))))
nodes_info <- data.frame(
  node = nodes,
  type = nodes_type
)