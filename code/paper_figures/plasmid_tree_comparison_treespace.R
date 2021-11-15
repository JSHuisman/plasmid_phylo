################################################################
## Plasmid tree comparison between methods
##
## Author: Jana S. Huisman
## Last Edit: March 2020
################################################################

library(ggplot2)
#library(reshape2)
library(dplyr)
library(tidyverse)

library("viridis")
library(RColorBrewer)

source('../R/Log_file_analysis_functions.R')

theme_set(theme_light(base_size = 24, base_family = 'sans'))

################################################################
base_dir='../../data'
###########################################################
library('treespace')

plasmids_to_keep <- read_csv('../../data/plasmids_to_keep.csv')
methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler',
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid',
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods

###########################################################

plasmid_path <- function(method, plasmid = 'Col_MG828'){
  path <- paste0(base_dir, '/plasmid_trees/', method, '_', plasmid, '_subset.trees')
  # by selecting subset, we exclude the fixed set
  return(path)
  }

read_trees <- function(x){
  trees <- try(read.nexus(plasmid_path(x, plasmid = plasmid))[seq(900, 9001, 100)])
  if (class(trees) == 'try-error'){
   trees <- NULL
  }
  return(trees)
}

ntrees = length(seq(900, 9001, 100))

for (plasmid in plasmids_to_keep$name){

  plasmid_trees <- lapply(methods, FUN = read_trees)
  plasmid_ntips <- sapply(1:length(plasmid_trees), function(i){length(plasmid_trees[[i]]$STATE_0$tip.label)})
  methods_included <- rep(1:9,1)[plasmid_ntips == median(plasmid_ntips)]

  if (median(plasmid_ntips) > 0){

    all_plasmid_trees <- plasmid_trees[[methods_included[1]]]
    for (i in methods_included[2:length(methods_included)]){
      all_plasmid_trees <- c(all_plasmid_trees, plasmid_trees[[i]])
    }

    if (all_plasmid_trees[[1]]$Nnode > 2){

      test <- findGroves(all_plasmid_trees, nf = 2, nclust = 2)

      #scree.size determines size of histogram inset
      # d in upper corner indicates distance between grid lines; can add ylim/xlim to change
      p <- plotGroves(test$treespace$pco, groups = sort(rep(methods[methods_included], ntrees)),
                      type="ellipse", scree.size = 0)#, lab.show=TRUE)
      print(p)
      #plotGroves(test, type="ellipse", lab.show=TRUE)
      quartz.save(paste0("../../figures/", plasmid,
                         "_treespace.pdf"), type = 'pdf', dpi = 100)
      rm(test)
    }
  }
}

