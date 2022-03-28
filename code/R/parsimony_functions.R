###### Gene Presence #####################################

check_gene_presence <- function(x){as.numeric(any(x != ''))}

read_tip_data <- function(method, base_dir = '../../data'){
  file_path <- paste0(base_dir, '/gene_summaries/', method, '.csv')
  tip_data <- read.csv(file_path, stringsAsFactors = F)

  compressed_tip_data <- tip_data %>%
    group_by(sample) %>%
    summarise_at(vars(-group_cols()), check_gene_presence)
  return(compressed_tip_data)
}

####### Get parsimony permutation #########################
 
# Used in compute_pars_for_gene()
get_parsimony_permutations <- function(tree, tip_states, n_perm = 1000){

  # this only finishes in appreciable time if the tree is small (otherwise just use n_perm)
  n_total_perm <- multinom(tip_states, counts = FALSE, useDouble = TRUE)
  n_true_perm <- min(n_total_perm, n_perm)

  n_pres = sum(tip_states)

  permuted_tip_states <- matrix(data = 0, nrow = n_true_perm, ncol = length(tip_states))
  for(i in 1:n_true_perm){
    permuted_tip_states[i,sample(1:length(tip_states), n_pres)] <- 1
  }

  colnames(permuted_tip_states) <- tree$tip.label

  parsimony.v <- vector(mode='numeric', length=dim(permuted_tip_states)[1])
  for (i_row in 1:dim(permuted_tip_states)[1]){
    data <- as.phyDat(as.factor(permuted_tip_states[i_row,]))
    parsimony.v[i_row]<- parsimony(tree, data)
  }

  return(parsimony.v)
}

# ##### Compute Parsimony Dataframes ###############################
compute_pars_for_gene <- function(method, tree, gene){
  tip_data <- read_tip_data(method)

  no_data_tips <- setdiff(tree$tip.label, tip_data$sample)
  if(length(no_data_tips)>0){
    extra_df <- as.data.frame(matrix(data = 0, nrow = length(no_data_tips), ncol = dim(tip_data)[2],
                                     dimnames = list(samples = no_data_tips, attributes = colnames(tip_data))))
    extra_df['sample'] <- no_data_tips
    tip_data <- rbind(tip_data, extra_df)
  }

  state_of_interest_id <- grep(gene, sub('\\.', '_', colnames(tip_data))) #paste0(gene, '$')
  if (length(state_of_interest_id)>1){
    print(paste0('Gene: ', gene, ' Method: ', method))
  }

  tip_states <- unlist(sapply(tree$tip.label,
                              function(x){unname(tip_data[tip_data$sample == x, state_of_interest_id[1]]) }))

  max_pars = min(sum(tip_states), length(tip_states) - sum(tip_states))

  if (any(is.na(tip_states))){
    return(list(perm = data.frame(), real = data.frame()))
  } else {
    parsimony.v <- get_parsimony_permutations(tree, tip_states)
    method_pars_df <- data.frame('parsimony' = parsimony.v, 'norm_pars' = parsimony.v/max_pars,
                                 'method' = method, 'gene' = gene)

    real_pars <- parsimony(tree, as.phyDat(as.factor(tip_states)))
    real_method_pars_df <- data.frame('parsimony' = real_pars, 'norm_pars' = real_pars/max_pars,
                                      'method' = method, 'gene' = gene)

    return(list(perm = method_pars_df, real = real_method_pars_df))
  }
}

compute_pars_dfs <- function(method_list, gene_list, plasmid = NULL){
  parsimony_df <- data.frame()
  real_parsimony_df <- data.frame()
  for (gene in gene_list){
    for (method in method_list){

      if (is.null(plasmid)){
        # from tree_reading_functions
        tree <- read_cgmlst_tree(method)
        tree_type = 'cgMLST'
      } else {
        # from tree_reading_functions
        tree <- read_plasmid_mcc_tree(method, plasmid)
        tree_type = plasmid
      }

      new_pars = compute_pars_for_gene(method, tree, gene)
      new_pars$perm['tree'] <- tree_type
      new_pars$real['tree'] <- tree_type

      parsimony_df <- rbind(parsimony_df, new_pars$perm)
      real_parsimony_df <- rbind(real_parsimony_df, new_pars$real)
    }
  }
  return(list(perm_parsimony = parsimony_df, real_parsimony = real_parsimony_df))
}

compute_multitree_pars_dfs <- function(method_list, gene_list, plasmid_list, cgmlst = TRUE){

  if (cgmlst){
    pars_dfs <- compute_pars_dfs(methods, gene_list, plasmid = NULL)
  } else {
    pars_dfs <- data.frame()
  }

  for (plasmid in plasmid_list){
    new_pars_dfs <- compute_pars_dfs(methods, gene_list, plasmid = plasmid)

    pars_dfs$perm_parsimony <- rbind(pars_dfs$perm_parsimony, new_pars_dfs$perm_parsimony)
    pars_dfs$real_parsimony <- rbind(pars_dfs$real_parsimony, new_pars_dfs$real_parsimony)
  }

  return(pars_dfs)
}


###########################################################
filter_pars_dfs <- function(pars_dfs, filter_expr){
  subset_perm_pars_dfs <- pars_dfs$perm_parsimony %>%
    filter(eval(parse(text=filter_expr)))
  subset_real_pars_dfs <- pars_dfs$real_parsimony %>%
    filter(eval(parse(text=filter_expr)))
  subset_pars_dfs <- list(perm_parsimony = subset_perm_pars_dfs,
                          real_parsimony = subset_real_pars_dfs)
  return(subset_pars_dfs)
}

# #### Compute p value dfs ##################################

pars_pvalue_from_df <- function(pars_dfs, tree_type, sel_gene, sel_method){

  perm_pars_vec <- pars_dfs$perm_parsimony %>%
    filter(tree == tree_type, gene == sel_gene, method == sel_method) %>%
    pull(parsimony)

  real_pars_val <- pars_dfs$real_parsimony %>%
    filter(tree == tree_type, gene == sel_gene, method == sel_method) %>%
    pull(parsimony)

  pars_density <- ecdf(perm_pars_vec)

  return(pars_density(real_pars_val))
}

compute_pars_pvalue_df <- function(pars_dfs, genes_list, method_list){

  trees <- unique(pars_dfs$real_pars$tree)
  pars_pvalue_df <- data.frame()
  for (tree_type in trees){
    for (sel_gene in genes_list){
      for (sel_method in method_list){
        pval <- pars_pvalue_from_df(pars_dfs, tree_type, sel_gene, sel_method)
        pval_df <- data.frame('p_value' = pval, 'tree' = tree_type,
                              'method' = sel_method, 'gene' = sel_gene)
        pars_pvalue_df <- rbind(pars_pvalue_df, pval_df)
      }
    }
  }
  return(pars_pvalue_df)
}

### Plot parsimony plots ##################################

plot_pars_densities <- function(pars_dfs, set_ncol, max_pars = 12, 
                                adjust = 3, dot_dist = 0.05){
  set_nrow = ceiling(length(unique(pars_dfs$perm_parsimony$gene))/set_ncol)
  
  ggplot(data=pars_dfs$perm_parsimony, aes(x=parsimony)) + 
    geom_density(aes(x=parsimony, colour = method), fill = "white", alpha = .3, adjust = adjust) +
    theme(axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold")) +
    scale_x_continuous(breaks = pretty(c(0.5, max_pars+0.5), n = 5), limits = c(0.5, max_pars+0.5)) +
    scale_colour_discrete(labels = method_labels) +  
    geom_point(data = pars_dfs$real_parsimony, position = position_stack(vjust = 1), show.legend = FALSE,
               alpha = .8, size = 3,
               mapping = aes(x = parsimony, colour = method, y = rep(dot_dist, dim(pars_dfs$real_parsimony)[1]) )) +
    facet_wrap(facets = vars(gene), nrow = set_nrow, ncol = set_ncol) + 
    labs(x = 'Parsimony score', y = 'Density', colour = 'Assembly method') +
    theme(legend.position = 'bottom') + guides(colour = guide_legend(nrow = 3))
}

plot_resgene_scan <- function(pars_dfs, sorted_spec_res_genes, x_var = 'norm_pars',
                              facet_row_var = 'method', facet_col_var = 'gene',
                              colour_var = 'tree', colour_name = 'Tree type',
                              max_pars = 12, method_labels = method_labels){
  pars_dfs$perm_parsimony$gene <- factor(pars_dfs$perm_parsimony$gene, levels = sorted_spec_res_genes)
  pars_dfs$real_parsimony$gene <- factor(pars_dfs$real_parsimony$gene, levels = sorted_spec_res_genes)
  
  pars_plot <- ggplot(data=pars_dfs$perm_parsimony, aes(x= get(x_var) ) ) + 
    geom_density(aes(x = get(x_var), colour = tree), alpha = .3, fill = "white", adjust = 3, show.legend = FALSE) + 
    geom_point(data = pars_dfs$real_parsimony, 
               mapping = aes(x = get(x_var), colour = get(colour_var), 
                             y = rep(0.2, dim(pars_dfs$real_parsimony)[1]) ),
               position = position_stack(vjust = 1), 
               alpha = .8, size = 3) +
    geom_vline(data = pars_dfs$real_parsimony, show.legend = FALSE, 
               mapping = aes(xintercept = get(x_var), colour = get(colour_var) ) ) +
    facet_grid(rows = vars(get(facet_row_var)), cols = vars(get(facet_col_var)),
               labeller = labeller( .rows = method_labels)) + 
    theme(axis.text=element_text(size=24), axis.title=element_text(size=24,face="bold"),
          legend.position = 'bottom') +
    labs(x = 'Parsimony score', y = 'Density', colour = colour_name) +
    theme_classic() #+ guides(colour = guide_legend(nrow = 1))
  
  if (x_var == 'parsimony'){
    pars_plot <- pars_plot +
      scale_x_continuous(breaks = pretty(c(0.5, max_pars+0.5), n = 5), limits = c(0.5, max_pars+0.5))
  }
  
  return(pars_plot)
}
