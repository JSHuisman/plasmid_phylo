##### Reading Trees #######################################

rename_tip <- function(tree){
  tree$tip.label <- gsub(".fasta", '', tree$tip.label)
  tree$tip.label <- gsub("_[[:digit:]]*$", '', tree$tip.label)
  return(tree)
}

read_cgmlst_tree <- function(method, base_dir = '../../data'){
  tree_file <- paste0(base_dir, '/chrom_trees/', method, '_MCC.tree')
  tree <- read.nexus(file = tree_file, tree.names = NULL) 
  tree <- rename_tip(tree)
  return(tree)  
}

read_plasmid_mcc_tree <- function(method, plasmid, 
                                  base_dir = '../../data'){
  tree_file <- paste0(base_dir, '/plasmid_trees/', method, '_', plasmid, '_MCC.trees')
  tree <- read.nexus(file = tree_file, tree.names = NULL)
  tree <- rename_tip(tree)
  return(tree)  
}

###########################################################
get_plasmid_tips <- function(base_dir, plasmid){
  tip_file <- paste0(base_dir, '/plasmid_trees/subsetted_tips/', plasmid, '.txt')
  tips = read.csv(tip_file, header = FALSE, col.names = 'tips')
  tips['taxa'] = gsub("_[[:digit:]]$", '', tips$tips)
  if (anyDuplicated(tips['taxa']) > 0){
    tips <- tips %>% distinct(taxa, .keep_all	= TRUE)
    #  stop("Plasmid tree contains duplicated taxa")
  }
  return(tips)
}

drop_tips_from_tree <- function(tree, tips_to_drop){
  if(class(tree) == 'phylo'){
    new_tree <- drop.tip(tree, tips_to_drop)
  } else if (class(tree) == 'multiPhylo'){
    new_tree <- lapply(tree, drop.tip, tips_to_drop)
    new_tree <- do.call(c, new_tree)
  }
  return(new_tree)
}
