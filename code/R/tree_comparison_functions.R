################################################################
## Tree comparison functions
##
## Author: Jana S. Huisman
## Last Edit: March 2020
################################################################

plot_by_method <- function(clade_df, name_y = 'Illumina_SPAdes', 
                           color_var = 'plasmid', facet_var = 'assembly_method', method_labels){
  nfacets <- dim(unique(clade_df[, facet_var]))[1]
  
  plot_by_method <- ggplot(clade_df, aes(x=clade_prob, y=get(name_y), color=get(color_var))) + 
    geom_point(size = 4) +
    #position = position_jitter(width=0.01), 
    lims(x = c(-0.1,1.1), y = c(-0.1,1.1)) +
    labs(y = paste0("Clade Probabilities in ", name_y, " tree posterior"),
         x = "Clade Probabilities in 2nd tree posterior", color = "2nd Assembly \nMethod") +
    geom_abline(slope=1, intercept=0, alpha=0.3) + 
    facet_wrap(facets = facet_var, nrow = 2, ncol = ceiling(nfacets/2)) +
    theme(legend.position = 'bottom')
  return(plot_by_method)
}

# concat_plasmid_clade_probs_2 <- function(input_folder, methods){
#   method_clade_df <- data.frame()
#   for (method_y in methods){
#     
#     for (method_x in setdiff(methods, method_y)){
#       
#       for (alignment_method in c('roary')){
#         #y vs x is stored in methodx/methody for historical reasons
#         plasmid_path <- paste0(input_folder, '/', method_x, '/', method_y, '_', alignment_method)
#         plasmid_files <- list.files(path = plasmid_path)
#         plasmids <- gsub('.txt', '', plasmid_files)
#         for (plasmid in plasmids){
#           clades <- try(read.csv(paste0(plasmid_path, '/', plasmid, '.txt'), 
#                                  sep=' ', stringsAsFactors = FALSE))
#           
#           if(class(clades) != 'try-error'){
#             colnames(clades) <- c('Clade', 'clade_prob_y', 'clade_prob_x')  
#             
#             clades['method_y'] <- method_y
#             clades['method_x'] <- method_x
#             clades['alignment_method'] <- alignment_method
#             clades['plasmid'] <- plasmid
#             method_clade_df <- rbind(method_clade_df, clades)
#           }
#         }
#       }
#     }
#   }
#   return(method_clade_df)
# }

################################################################

clade_heatmap<- function(clade_df, method_labels, col_order = NULL ){
  subset_df = clade_df
  subset_df$Clade <- gsub("ESBL|.fasta|\\{|\\}", "", subset_df$Clade)
  clade_size_order = order(unlist(lapply(unique(subset_df$Clade), function(x){length(unlist(strsplit(x, ',')))})))
  subset_df$Clade <- factor(subset_df$Clade, levels = unique(subset_df$Clade)[clade_size_order])
  
  # column, i.e. method, order
  if (!is.null(col_order)){
    subset_df$method <- factor(subset_df$method, levels = col_order)
  }
  
  ggplot(subset_df, 
         #aes(y = Clade, fill = clade_prob_y, x = method_y)) +
         aes(y = Clade, fill = clade_prob, x = method)) +
    geom_tile() +
    labs(x = "Assembly method", y = "Clade", fill = "Clade \nprobability") +
    #scale_fill_gradient2(midpoint = 0.5) +
    scale_fill_viridis(option="viridis", discrete=FALSE) +
    scale_x_discrete(labels = method_labels) +
    theme(axis.text.y = element_blank(), #element_text(size = 15),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
          axis.title=element_text(size=24, face="bold"))
}


clade_heatmap_by_plasmid <- function(clade_df, plasmid, alignment, method_labels){
  subset_df = clade_df[clade_df$plasmid == plasmid & clade_df$alignment_method == alignment,]
  
  # compute column order
  subset_by_row = pivot_wider(subset_df, names_from = Clade, values_from = clade_prob, id_cols = method)
  dist_m <- dist(subset_by_row[,2:dim(subset_by_row)[2]], method = "euclidean")
  hierarch_clust <- hclust(dist_m)
  col_order <- subset_by_row[hierarch_clust$order,1]
  
  clade_heatmap(subset_df, method_labels, col_order$method) 
}

################################################################
# 
# concat_clade_probs <- function(input_folder, methods){
#   method_clade_df <- data.frame()
#   for (method_y in methods){
#     
#     for (method_x in setdiff(methods, method_y)){
#       # method y vs x is stored in method_x/method_y_cgmlst.txt for historical reasons
#       clades <- try(read.csv(paste0(input_folder, '/', method_x, '/', method_y, '_cgmlst.txt'), 
#                              sep=' ', stringsAsFactors = FALSE))
#       
#       if(class(clades) != 'try-error'){
#         colnames(clades) <- c('Clade', 'clade_prob_y', 'clade_prob_x')  
#         
#         clades['method_y'] <- method_y
#         clades['method_x'] <- method_x
#         
#         method_clade_df <- rbind(method_clade_df, clades)
#       }
#     }
#   }
#   return(method_clade_df)
# }

plot_cgmlst_clades <- function(clade_df, methods){
  
  plot_by_method <- ggplot(clade_df, aes(x=clade_prob_x, y=clade_prob_y, color=method_x)) + 
    geom_point(size = 4, alpha = 0.75) +
    #position = position_jitter(width=0.01),
    lims(x = c(-0.1,1.1), y = c(-0.1,1.1)) +
    labs(y = paste0("Clade Probabilities in 1st tree posterior"),
         x = "Clade Probabilities in 2nd tree posterior", color = "2nd Assembly \nMethod") +
    scale_color_manual(labels = unname(methods),
                       breaks = names(methods),
                       values = brewer.pal(length(methods), 'Paired')) + 
    # inferno(length(methods)) # 'Set1' also nice
    geom_abline(slope=1, intercept=0, alpha=0.3) + 
    facet_wrap(facets = vars(method_y), nrow = 3, ncol = 3,
               labeller = labeller(method_y = methods)) +
    theme(legend.position = 'bottom')
  return(plot_by_method)
}
################################################################

.base_breaks <- function(n = 5){
  function(x) {
    grDevices::axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    #pretty(gamma_max_growth, n = 5)
  }
}


