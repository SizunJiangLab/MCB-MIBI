# Bunch of Functions written by Sizun Jiang for interfacing with MIBI data
# Plot clustering heatmaps from:
# CyTOF workflow, Mark Robinson lab

plot_clustering_heatmap_wrapper <- function(expr, expr01, cell_clustering, color_clusters, cluster_merging = NULL){
  # Function to plot pretty heatmaps, according to a expression df and cell clusters
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  colnames(expr_median) = c("cell_clustering", colnames(expr))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  colnames(expr01_median) = c("cell_clustering", colnames(expr01))
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)")
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering)) 
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)] 
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  #Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6,  legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
}

plot_clustering_heatmap_wrapper_mean <- function(expr, expr01, cell_clustering, color_clusters, cluster_merging = NULL){
  # Function to plot pretty heatmaps, according to a expression df and cell clusters
  # Calculate the mean expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  colnames(expr_median) = c("cell_clustering", colnames(expr))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  colnames(expr01_median) = c("cell_clustering", colnames(expr01))
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)")
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering)) 
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)] 
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  #Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6,  legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
}


plot_clustering_heatmap_wrapper2 <- function(expr, expr01,
                                             lineage_markers, functional_markers = NULL, sample_ids = NULL,
                                             cell_clustering, color_clusters, cluster_merging = NULL,
                                             plot_cluster_annotation = TRUE, title = NULL){
  # Calculate the median expression of lineage markers
  expr_median <- data.frame(expr[, lineage_markers],
                            cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  colnames(expr_median) = c("cell_clustering", colnames(expr[, lineage_markers]))
  expr01_median <- data.frame(expr01[, lineage_markers],
                              cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  colnames(expr01_median) = c("cell_clustering", colnames(expr01[, lineage_markers]))
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, lineage_markers], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, lineage_markers])
  # Median expression of functional markers in each sample per cluster
  expr_median_sample_cluster_tbl <- data.frame(expr01[, functional_markers,drop = FALSE], 
                                               sample_id = sample_ids, 
                                               cluster = cell_clustering) %>% group_by(sample_id, cluster) %>% summarize_all(funs(median))
  colnames(expr_median_sample_cluster_tbl) = c("sample_id", "cluster", functional_markers)
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , "%)")
  ### Annotation for the original clusters
  annotation_row1 <- data.frame(Cluster = factor(expr01_median$cell_clustering)) 
  color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)] 
  names(color_clusters1) <- levels(annotation_row1$Cluster)
  ### Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
    annotation_row2 <- data.frame(Cluster_merging =
                                    factor(cluster_merging$new_cluster[mm]))
    color_clusters2 <- color_clusters[1:nlevels(annotation_row2$Cluster_merging)]
    names(color_clusters2) <- levels(annotation_row2$Cluster_merging)
  }
  ### Heatmap annotation for the original clusters
  ha1 <- Heatmap(annotation_row1, name = "Cluster",
                 col = color_clusters1, cluster_columns = FALSE, cluster_rows = cluster_rows, row_dend_reorder = FALSE, 
                 show_row_names = FALSE, width = unit(0.5, "cm"), 
                 rect_gp = gpar(col = "grey"))
  ### Heatmap annotation for the merged clusters
  if(!is.null(cluster_merging)){
    ha2 <- Heatmap(annotation_row2, name = "Cluster \nmerging",
                   col = color_clusters2, cluster_columns = FALSE, 
                   cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_row_names = FALSE, width = unit(0.5, "cm"),
                   rect_gp = gpar(col = "grey"))
  }
  ### Cluster names and sizes - text
  ha_text <- rowAnnotation(text = row_anno_text(labels_row,
                                                gp = gpar(fontsize = 6)), width = max_text_width(labels_row))
  ### Cluster sizes - barplot
  ha_bar <- rowAnnotation("Frequency (%)" = row_anno_barplot(x = clustering_prop, border = FALSE, axis = TRUE,
                                                             axis_gp = gpar(fontsize = 5), gp = gpar(fill = "#696969", col = "#696969"), bar_width = 0.9), 
                          width = unit(0.7, "cm"), show_annotation_name = TRUE, annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"), 
                          annotation_name_gp = gpar(fontsize = 7))
  ### Heatmap for the lineage markers
  ht1 <- Heatmap(expr_heat, name = "Expr", column_title = "Lineage markers", col = color_heat, cluster_columns = FALSE, 
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks,labels = legend_breaks, color_bar = "continuous"),
                 show_row_names = FALSE, row_dend_width = unit(2, "cm"),
                 rect_gp = gpar(col = "grey"), column_names_gp = gpar(fontsize = 8))
  if(plot_cluster_annotation){
    draw_out <- ha1
  }else{
    draw_out <- NULL
  }
  if(!is.null(cluster_merging)){
    draw_out <- draw_out + ha2 + ht1 + ha_bar + ha_text
  }else{
    draw_out <- draw_out + ht1 + ha_bar + ha_text
  }
  ### Heatmaps for the signaling markers
  if(!is.null(functional_markers)){
    for(i in 1:length(functional_markers)){
      ## Rearange so the rows represent clusters
      expr_heat_fun <- as.matrix(dcast(expr_median_sample_cluster_tbl[, c("sample_id", "cluster", functional_markers[i])],
                                       cluster ~ sample_id, value.var = functional_markers[i])[, -1])
      draw_out <- draw_out + Heatmap(expr_heat_fun, column_title = functional_markers[i], col = color_heat, cluster_columns = FALSE, 
                                     cluster_rows = cluster_rows, row_dend_reorder = FALSE, show_heatmap_legend = FALSE, 
                                     show_row_names = FALSE, rect_gp = gpar(col = "grey"), column_names_gp = gpar(fontsize = 8))
    }
  }
  draw(draw_out, row_dend_side = "left", column_title = title)
}

## Function to color per location of newLmod
color_overlay = function(newLmodlist, expr, cell_clustering, color_clusters){
  # newLmodlist is a list of newLmods from segParams (segmentationParams.mat)
  # List contains each FOV (or Point number)
  # expr = expr2_df
  # cell_clustering = cell_clustering1
  # color_clusters = col_vector
  
  ## Subset the data into a list of points
  expr$CellType = cell_clustering
  expr_list = split(expr, expr$PointNum)
  
  unique_FOVs = sort(as.numeric(names(expr_list)), decreasing = F)
  all_ggplots = lapply(unique_FOVs, function(fov){
    # print(fov)
    newLmod = newLmodlist[[as.character(fov)]]
    # ## Try for Point1
    # expr_point = expr_list[[1]]
    expr_point = expr_list[[as.character(fov)]]
    unique_clusters = sort(unique(expr_point[["CellType"]]), decreasing = F)
    # GGplot compatible df with row, col, cluster
    cluster_df = do.call(rbind, lapply(unique_clusters, function(c){
      # print(c)
      # w1 is getting cell label in image from cluster
      w1 = expr_point[["cellLabelInImage"]][which(expr_point[["CellType"]] == c)]
      # w2 is getting xy coordinates from cell label in image, from newLmod
      w2 = data.frame(which(newLmod == w1, arr.ind = T))
      w2$CellType = as.factor(c)
      w2$row = dim(newLmod)[1] - w2$row
      return(w2)
    }))
    
    ggplot(cluster_df, aes_string(x="col", y="row", color = "CellType")) +
      geom_point(size=2, shape = 15) +
      # geom_point(shape = 1, size=3, colour = "black") +
      #guides(colour=guide_legend(override.aes=list(size=6))) +
      xlab("") + ylab("") +
      ggtitle(paste("Cell Type Plot for FOV number:", fov)) + 
      theme_light()+
      # theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            legend.box = "horizontal") +
      ylim(c(0, dim(newLmod)[2])) +
      xlim(c(0, dim(newLmod)[1])) +
      scale_colour_manual(values = mycolor)
  })
}

plot_gradient_mean=function(data, var_cluster, midscale = F)
{
  w = which(colnames(data) == var_cluster)
  mid = mean(data[,w])
  df = data.frame(x = data[["UMAP1"]],
                  y = data[["UMAP2"]],
                  var = data[,w])
  plt = ggplot(df) + coord_fixed(ratio=graphical.ratio) + 
    geom_point(aes(x=x, y=y, color=var),cex = 1.5) + 
    labs(x = "UMAP 1", y = "UMAP 2") +
    # ggplot(df, aes_string(x="UMAP1", y="UMAP2", color=var))+
    #   geom_point(size=0.25) +
    #   # guides(colour=guide_legend(override.aes=list(size=6))) +
    #   xlab("") + ylab("") +
    #ggtitle(as.character(var_cluster)) +
    theme_light(base_size=20) +
    # theme(axis.text.x=element_blank(),
    #       axis.text.y=element_blank(),
    #       # legend.direction = "horizontal", 
    #       # legend.position = "bottom",
    #       legend.box = "horizontal",
    #       legend.key.size = unit(1,"line")) + 
    ggtitle(paste0(var_cluster, " Expression Level"))
  if (midscale == F){
    plt + scale_color_gradient2(low="blue", mid="white",
                                high="red", space ="Lab")}
  else {plt+scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                                  high="red", space ="Lab")}
  
  # No midpoint scaling
  
  
  #scale_color_gradient2(palette) 
}

# plot_discrete=function(data, var_cluster)
#   # Useless and old now
# {
#   mycolor = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))
#   df = data.frame(data["UMAP1"],
#                   data["UMAP2"],
#                   data[var_cluster])
#   ggplot(df) + coord_fixed(ratio=graphical.ratio) + 
#     # geom_point(aes(x=UMAP1, y=UMAP2, color=var_cluster),cex = 1.5) + 
#     # labs(x = "UMAP1", y = "UMAP2") +
#     # # ggplot(df, aes_string(x="UMAP1", y="UMAP2", color=var))+
#     # #   geom_point(size=0.25) +
#     # #   # guides(colour=guide_legend(override.aes=list(size=6))) +
#     # #   xlab("") + ylab("") +
#     # #ggtitle(as.character(var_cluster)) +
#     # theme_light(base_size=20) +
#     # scale_color_manual(values=mycolor)
#     # # ggtitle(title = "Counts of Cell Types")
#     # 
#     # 
#   geom_point(aes(x=UMAP1, y=UMAP2, color=as.factor(var_cluster)),cex = 1.5) + 
#     guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) +
#     labs(x = "UMAP 1", y = "UMAP 2",title = paste0("UMAP representation colored by ",var_cluster), 
#          color = "var_cluster") + theme_bw() + 
#     scale_color_manual(values=mycolor) + geom_label_repel(data = umap.cent.celltype, aes(x = x, y = y, label = var_cluster))
# }

make_cellcluster_layers = function(expr_df_anon, newLmod, plot_cluster_name, mycolor, cell_types){
  require(EBImage)
  # expr_df_anon = FOV specific df containing the cellLabelInImage, channels and more importantly, a column for "cell_type_label" (plot_cluster_name)
  # newLmod is from the .mat output segmentationParams.mat
  # my color is a vector of colors, eg "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"
  # outputs an EBImage class, essentially an array of x by y by z, where z is a matrix consisting of the assigned color for each cluster, and white for cell boundaries and black for no cells
  # Get cell types
  if(missing(plot_cluster_name)){
    plot_cluster_name = "cell_type_label"
  }
  if(missing(cell_types)){
    cell_types = unique(expr_df_anon[["cell_type_label"]])[order(unique(expr_df_anon[[plot_cluster_name]]), decreasing = F)]
  }
  # Map cell type to color
  color_matrix = matrix(0, nrow = dim(newLmod)[1], ncol = dim(newLmod)[2])
  for(i in 1:length(cell_types)){
    ilayer = expr_df_anon[["cellLabelInImage"]][which(expr_df_anon[[plot_cluster_name]] == cell_types[i])]
    color_matrix[newLmod %in% ilayer] = mycolor[i]
  }
  color_matrix[newLmod == 0] = "#000000" # Black for boundaries
  color_matrix[color_matrix == 0] = "#FFFFFF" # White for no cells
  return(Image(t(color_matrix)))
}

euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

get_expt_interaction = function(expr_df_anon, colname_celltype, radius, mc.cores){
  # Takes a df of many cells and their expressions, expr_df_anon, containing 
  # 1. colname_celltype, the column name for the cell type label. Default i cell_type_label
  # 2. x_cent and y_cent for the coordinates (in pixels)
  # radius is the number of pixels apart the center of cells have to be to consider as "interacting"
  # mc.cores is the number of cores to use in mclapply
  # Returns a list of each FOV, as defined by PointNum, and a 2D array of cell_type x cell_type
  require(parallel)
  require(spatstat)
  if(missing(colname_celltype)){
    colname_celltype = "cell_type_label"}
  if(missing(radius)){
    radius = 20}
  # if(missing(cols_to_average)){
  #   cols_to_average = NULL}
  # if(missing(pn)){
  #   pn = 1000}
  # if(missing(mc.cores)){
  #   mc.cores = 4}
  # Split into list by pointnum
  expr_list = split(expr_df_anon, expr_df_anon$PointNum)
  celltypes = unique(expr_df_anon[[colname_celltype]])[order(unique(expr_df_anon[[colname_celltype]]), decreasing = F)]
  # Find expected interactions for each AB pair
  expt_array = mclapply(expr_list, function(pointcells){
    distances = crossdist(pointcells$x_cent, pointcells$y_cent, pointcells$x_cent, pointcells$y_cent, squared = F)
    mat = matrix(0, nrow = length(celltypes), ncol = length(celltypes))
    for(i in 1:length(celltypes)){
      for(j in 1:length(celltypes)){
        A = celltypes[i]
        B = celltypes[j]
        distances_AB = distances
        w1 = which(pointcells[[colname_celltype]] == A)
        w2 = which(pointcells[[colname_celltype]] == B)
        if(length(w1) > 0 && length(w2) > 0){
          kick = setdiff(1:dim(pointcells)[1], c(w1,w2))
          distances_AB[kick,] = 0
          distances_AB[,kick] = 0
          exp_interaction = sum(as.vector(distances_AB) < radius & as.vector(distances_AB) > 0)/2  
        }else{
          exp_interaction = NA}
        mat[i,j] = exp_interaction
      }
    }
    return(mat)
  }, mc.cores = mc.cores)
  return(expt_array)
}

# find_interactions = function(expr_list, colname_celltype, mc.cores){
#   # Takes in a list, each containing a df of interactions
#   if(missing(mc.cores)){
#     mc.cores = 4}
#   expected_interactions_FOV = mclapply(expr_list, function(pointcells){
#     return(df_interaction(pointcells, colname_celltype, radius))
#   }, mc.cores = mc.cores)
# }

# df_interaction = function(pointcells, colname_celltype, radius){
#   if(missing(radius)){
#     radius = 20}
#   if(missing(colname_celltype)){
#     colname_celltype = "cell_type_label"}
#   distances = crossdist(pointcells$x_cent, pointcells$y_cent, pointcells$x_cent, pointcells$y_cent, squared = F)
#   # Retrieval of crossdist between 2 cell types of interest
#   # Keep the crossdist matrix layout/dimension, just make everything else = 0
#   # No need to recalculate the distance matrix, its going to be the same however labels are shuffled
#   A_vec = vector()
#   B_vec = vector()
#   exp_freq_vec = vector()
#   for(i in 1:length(celltypes)){
#     for(j in 1:length(celltypes)){
#       A = celltypes[i]
#       B = celltypes[j]
#       distances_AB = distances
#       w1 = which(pointcells[[colname_celltype]] == A)
#       w2 = which(pointcells[[colname_celltype]] == B)
#       if(length(w1) > 0 && length(w2) > 0){
#         kick = setdiff(1:dim(pointcells)[1], c(w1,w2))
#         distances_AB[kick,] = 0
#         distances_AB[,kick] = 0
#         # distances_AB[!(pointcells[[colname_celltype]] %in% c(A, B)), ] = 0
#         # distances_AB[,!(pointcells[[colname_celltype]] %in% c(A, B))] = 0
#         # Double counting from symmetrical matrix
#         exp_interaction = sum(as.vector(distances_AB) < radius & as.vector(distances_AB) > 0)/2  
#       }else{
#         exp_interaction = NA
#       }
#       A_vec = c(A_vec, A)
#       B_vec = c(B_vec, B)
#       exp_freq_vec = c(exp_freq_vec, exp_interaction)
#     }
#   }
#   df = data.frame(A = A_vec, B = B_vec, exp_freq = exp_freq_vec)
# }

df_interaction_perm = function(pointcells, colname_celltype, radius, pn, mc.cores){
  require(bigstatsr)
  if(missing(radius)){
    radius = 20}
  if(missing(colname_celltype)){
    colname_celltype = "cell_type_label"}
  distances = crossdist(pointcells$x_cent, pointcells$y_cent, pointcells$x_cent, pointcells$y_cent, squared = F)
  # Retrieval of crossdist between 2 cell types of interest
  # Keep the crossdist matrix layout/dimension, just make everything else = 0
  # No need to recalculate the distance matrix, its going to be the same however labels are shuffled
  # Make a 3D Array for celltype A x celltype Bx pn
  # Order follows celltypes
  # nd_array = array(0, c(length(celltypes), length(celltypes), pn))
  shuff_df = pointcells
  nd_array = abind(mclapply(1:pn, function(n){
    mat = matrix(0, nrow = length(celltypes), ncol = length(celltypes))
    for(i in 1:length(celltypes)){
      for(j in 1:length(celltypes)){
        A = celltypes[i]
        B = celltypes[j]
        distances_AB = distances
        shuff_df[[colname_celltype]] = sample(pointcells[[colname_celltype]])
        w1 = which(shuff_df[[colname_celltype]] == A)
        w2 = which(shuff_df[[colname_celltype]] == B)
        if(length(w1) > 0 && length(w2) > 0){
          kick = setdiff(1:dim(shuff_df)[1], c(w1,w2))
          distances_AB[kick,] = 0
          distances_AB[,kick] = 0
          exp_interaction = sum(as.vector(distances_AB) < radius & as.vector(distances_AB) > 0)/2  
        }else{
          exp_interaction = NA}
        mat[i,j] = exp_interaction
      }
    }
    return(mat)
  }, mc.cores = mc.cores), along = 3)
  return(nd_array)
}

# get_perm_interaction = function(expr_df_anon, colname_celltype, radius, pn, mc.cores){
#   # Takes a df of many cells and their expressions, expr_df_anon, containing 
#   # 1. colname_celltype, the column name for the cell type label. Default is cell_type_label
#   # 2. x_cent and y_cent for the coordinates (in pixels)
#   # radius is the number of pixels apart the center of cells have to be to consider as "interacting", default is 20
#   # mc.cores is the number of cores to use in mclapply
#   # Returns a list of each FOV, as defined by PointNum, and a df of the permutated interactions between each pairwise cell type
#   require(parallel)
#   require(spatstat)
#   require(parallel)
#   require(spatstat)
#   if(missing(colname_celltype)){
#     colname_celltype = "cell_type_label"}
#   if(missing(radius)){
#     radius = 20}
#   if(missing(mc.cores)){
#     mc.cores = 4}
#   if(missing(pn)){
#     pn = 1000}
#   
#   # Split into list by pointnum
#   expr_list = split(expr_df_anon, expr_df_anon$PointNum)
#   celltypes = unique(expr_df_anon[[colname_celltype]])[order(unique(expr_df_anon[[colname_celltype]]), decreasing = F)]
#   
#   # Permutation test
#   perm_interactions_list = lapply(expr_list, function(pointcells){
#     shuff_df = pointcells
#     perm_interactions_FOV = mclapply(1:pn, function(i){
#       shuff_df[[colname_celltype]] = sample(pointcells[[colname_celltype]])
#       perm_interaction = df_interaction(shuff_df, colname_celltype, radius)
#     }, mc.cores = mc.cores)
#     # Breakdown list into dataframe of pn + 2 number of columns
#     perm_interactions_FOVdf = perm_interactions_FOV %>% reduce(left_join, by = c("A", "B"))
#     # perm_interactions_FOVmean = data.frame(perm_interactions_FOVdf[,1:2],
#     #                                       mean =rowMeans(perm_interactions_FOVdf[3:dim(perm_interactions_FOVdf)[2]]))
#     return(perm_interactions_FOVdf)
#   })
#   return(perm_interactions_list)
# }


get_perm_interaction_array = function(expr_df_anon, colname_celltype, radius, pn, mc.cores){
  # Takes a df of many cells and their expressions, expr_df_anon, containing
  # 1. colname_celltype, the column name for the cell type label. Default is cell_type_label
  # 2. x_cent and y_cent for the coordinates (in pixels)
  # radius is the number of pixels apart the center of cells have to be to consider as "interacting", default is 20
  # mc.cores is the number of cores to use in mclapply
  # Returns a list of containing a DF corresponding to the pairwise interactions, and
  # A list of each FOV, as defined by PointNum, containing a celltype * cell type * pn array of the permutated interactions between each pairwise cell type
  require(parallel)
  require(spatstat)
  
  # if(missing(colname_celltype)){
  #   colname_celltype = "cell_type_label"}
  # if(missing(radius)){
  #   radius = 20}
  # if(missing(mc.cores)){
  #   mc.cores = 4}
  # if(missing(pn)){
  #   pn = 1000}
  
  # Split into list by pointnum
  expr_list = split(expr_df_anon, expr_df_anon$PointNum)
  celltypes = unique(expr_df_anon[[colname_celltype]])[order(unique(expr_df_anon[[colname_celltype]]), decreasing = F)]
  
  # Permutation test
  perm_interactions_list = mclapply(expr_list, function(pointcells){
    perm_interaction = df_interaction_perm(pointcells, colname_celltype, radius, pn, mc.cores)
  }, mc.cores = mc.cores)
  return(perm_interactions_list)
}

# Pair wise correlation plots
# https://www.r-bloggers.com/2011/03/five-ways-to-visualize-your-pairwise-comparisons/
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}


#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
  for(i in rev(new_order)){
    x=relevel(x,ref=i)
  }
  return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
  # normalize data
  for(i in 1:length(data)){
    data[,i]=(data[,i]-mean(data[,i]))/sd(data[,i])
  }
  # obtain new data frame
  z=data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      temp=as.data.frame(cbind(x,y))
      temp=cbind(temp,names(data)[i],names(data)[j])
      z=rbind(z,temp)
      j=j+1
    }
  }
  names(z)=c('x','y','x_lab','y_lab')
  z$x_lab = ezLev(factor(z$x_lab),names(data))
  z$y_lab = ezLev(factor(z$y_lab),names(data))
  z=z[z$x_lab!=z$y_lab,]
  #obtain correlation values
  z_cor = data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      x_mid = min(x)+diff(range(x))/2
      y_mid = min(y)+diff(range(y))/2
      this_cor = cor(x,y)
      this_cor.test = cor.test(x,y)
      this_col = ifelse(this_cor.test$p.value<.05,'<.05','>.05')
      this_size = (this_cor)^2
      cor_text = ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b=as.data.frame(cor_text)
      b=cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor=rbind(z_cor,b)
      j=j+1
    }
  }
  names(z_cor)=c('cor','x_mid','y_mid','p','rsq','x_lab','y_lab')
  z_cor$x_lab = ezLev(factor(z_cor$x_lab),names(data))
  z_cor$y_lab = ezLev(factor(z_cor$y_lab),names(data))
  diag = z_cor[z_cor$x_lab==z_cor$y_lab,]
  z_cor=z_cor[z_cor$x_lab!=z_cor$y_lab,]
  #start creating layers
  points_layer = layer(
    geom = 'point'
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_line_layer = layer(
    geom = 'line'
    , geom_params = list(colour = 'red')
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_ribbon_layer = layer(
    geom = 'ribbon'
    , geom_params = list(fill = 'green', alpha = .5)
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  cor_text = layer(
    geom = 'text'
    , data = z_cor
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=cor
      , size = rsq
      , colour = p
    )
  )
  var_text = layer(
    geom = 'text'
    , geom_params = list(size=var_text_size)
    , data = diag
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=x_lab
    )
  )
  f = facet_grid(y_lab~x_lab,scales='free')
  o = opts(
    panel.grid.minor = theme_blank()
    ,panel.grid.major = theme_blank()
    ,axis.ticks = theme_blank()
    ,axis.text.y = theme_blank()
    ,axis.text.x = theme_blank()
    ,axis.title.y = theme_blank()
    ,axis.title.x = theme_blank()
    ,legend.position='none'
  )
  size_scale = scale_size(limits = c(0,1),to=cor_text_limits)
  return(
    ggplot()+
      points_layer+
      lm_ribbon_layer+
      lm_line_layer+
      var_text+
      cor_text+
      f+
      o+
      size_scale
  )
}

# For a DF of interest, given the anchor cells, find the sum and mean of the markers/celltype counts across the columns_to_use for finding the distance of markers from anchor cells
get_distance_matrix <- function(plot_df, anchor_cells, split_col1, split_col2, columns_to_use, micron_dist = 100, pixel_to_microns = 400/512, bins = 11, cores = 6){
  # plot_df: Your big DF to use
  # anchor_cells: boolean for which cells are teh anchor cells
  # columns_to_use: which columns to use for calculating the distance metric_types
  # split_col1/2 = which column to split by
  
  dist = ceiling(micron_dist * pixel_to_microns)
  
  dist_df = do.call(rbind,lapply(unique(plot_df[["PointNum"]]), function(fov){
    print(paste("Calculating FOV:",fov))
    pointdf = plot_df[plot_df$PointNum == fov,]
    # Go through anchorcells
    w = anchor_cells[plot_df$PointNum == fov]
    w = which(w == TRUE)
    # Ensure anchor cells present in FOV
    if(!is.empty(w)){
      # Loop through these cells and calculate distance of cells in a radius around
      df_keep = do.call(dplyr::bind_rows, mclapply(1:length(w), function(i){
        
        # df_keep = data.frame()
        # for(i in 1:length(w)){
        # print(i)
        # df_keep = lapply(1:length(w), function(i){
        # print(i)
        # Zero around cell of interest
        new_x = pointdf$x_cent - pointdf[w[i],][["x_cent"]]
        new_y = pointdf$y_cent - pointdf[w[i],][["y_cent"]]
        # keep = which(abs(new_x) <= radius & abs(new_y) <= radius)
        # This version uses a square
        # keep = which(abs(new_x) <= radius & abs(new_y) <= radius & pointdf[[marker_of_interest]] > 0)
        # This version uses radius. Includes self!
        dist_from_anchor = crossdist(new_x,new_y, 0,0, squared = F)
        keep = which(dist_from_anchor <= radius)
        # Keep dist_df by keep, rearranged by dist_from_anchor increasingly
        ordered = sort(dist_from_anchor[keep], decreasing = F, index.return = T) # This is the order of keep!
        
        # Bin data for each marker from these ordered cells
        # Loop through columns_to_use
        binned_cols = do.call(rbind, lapply(columns_to_use, function(c){
          # Make marker df of x = distance, y = expression/count
          marker_df = data.frame(dist = dist_from_anchor[keep[ordered$ix]],
                                 expr = pointdf[,c][keep[ordered$ix]])
          marker_df = transform(marker_df, bin = cut(dist, seq(0, micron_dist, length.out = bins), include.lowest = T))
          marker_sum = plyr::ddply(marker_df, "bin", summarize, totVal = sum(expr))
          marker_mean = plyr::ddply(marker_df, "bin", summarize, totVal = mean(expr))
          marker_positive = plyr::ddply(marker_df, "bin", summarize, totVal = mean(expr>0))
          
          df_ret = rbind(data.frame(t(marker_sum$totVal)),
                         data.frame(t(marker_mean$totVal)),
                         data.frame(t(marker_positive$totVal)))
          colnames(df_ret) = marker_sum$bin
          df_ret$Marker = colnames(pointdf)[c]
          df_ret$PointNum = fov
          df_ret$cellLabelInImage = pointdf$cellLabelInImage[w[i]]
          df_ret$Type = c("Sum", "Mean", "Positive")
          df_ret$Split1 = pointdf[w[i],][[split_col1]]
          df_ret$Split2 = pointdf[w[i],][[split_col2]]
          return(df_ret)
        }))
        # print(dim(binned_cols))
        # df_keep = rbind(df_keep,binned_cols)
        # }
        return(binned_cols)
      }, mc.cores = cores))
    }
    else{
      df_keep = data.frame()
    }
    return(df_keep)
  }))
  return(dist_df)
}


