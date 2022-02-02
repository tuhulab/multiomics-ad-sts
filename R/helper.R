plot_heatmap_from_se <- function(
  se, 
  column_names = FALSE, 
  row_names = FALSE,
  sample_in_col = TRUE,
  assay_no = 2){ # plot heatmap from SummarizedExperiment
  
  se_to_plot <- se
  
  heatmap_data_m <-
    scale(assay(se_to_plot, assay_no) %>% t() %>% log1p(), 
          center = TRUE, scale = TRUE)
  
  heatmap_data_m <- if(sample_in_col == TRUE){
    heatmap_data_m %>% t()
    } else {heatmap_data_m}
  
  sample_site_color <- 
    colData(se_to_plot) %>% as_tibble() %>% pull(biopsy_area) %>% unique()
  
  color <- RColorBrewer::brewer.pal(length(sample_site_color), "Pastel2")
  names(color) <- sample_site_color
  
  heatmap_annotation <- 
    ComplexHeatmap::HeatmapAnnotation(
      `Tissue type` = colData(se_to_plot) %>% as_tibble() %>% pull(skin_type), 
      `Anatomic region` = colData(se_to_plot) %>% as_tibble() %>% pull(biopsy_area),
      col = list(`Tissue type` = pca_color, 
                 `Anatomic region` = color),
      annotation_legend_param = list(`Tissue type` = list(nrow = 1, 
                                                          title_gp = gpar(fontsize = 16),
                                                          labels_gp = gpar(fontsize = 16)),
                                     `Anatomic region` = list(nrow = 1,
                                                              title_gp = gpar(fontsize = 16),
                                                              labels_gp = gpar(fontsize = 16))
      ))
  
  heatmap <- ComplexHeatmap::Heatmap(heatmap_data_m,
                                     name = "gene expression",
                                     column_dend_reorder = TRUE,
                                     # clustering_method_columns = "complete",
                                     show_column_names = column_names,
                                     show_row_names = row_names,
                                     show_row_dend = FALSE,
                                     # row_names_gp = gpar(fontsize = min(10, 800 * 1.5 / dim(heatmap_data_m)[1])),
                                     # row_annotation = heatmap_annotation,
                                     heatmap_legend_param = list(
                                       title_gp = gpar(fontsize = 16),
                                       labels_gp = gpar(fontsize = 14)
                                     ))
}