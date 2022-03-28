####################################
# Plotting Fig 1
#
# Author: J.S.Huisman, March 2022
####################################

library(ape)
library("ggplot2")
#install_github("YuLab-SMU/ggtree")
library("ggtree")
library(tidyverse)

####################################

theme_set(theme_light(base_size = 24, base_family = 'sans'))

methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler', 
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid', 
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods

base_dir = '../../data'
fig_folder <- '../../figures'

source('../R/tree_reading_functions.R')

####################################
# Reading in the cgmlst trees

tree_list = lapply(methods, function(method){read_cgmlst_tree(method, base_dir)})
class(tree_list) <- "multiPhylo"
names(tree_list) <- methods

colour_map <- function(x){
  if (x == 'Illumina_SPAdes'){
    return(factor("Illumina_SPAdes", levels = c("Illumina_SPAdes", "NP_Canu", "Other")))
  } else if (x == "NP_Canu"){
    return(factor("NP_Canu", levels = c("Illumina_SPAdes", "NP_Canu", "Other")))
  } else {
    return(factor("Other", levels = c("Illumina_SPAdes", "NP_Canu", "Other")))
    }
}

###########
# Joint cgMLST tree - correct tip labels
p3 <- ggtree(tree_list[c(1, 3:9)], layout='roundrect', 
             aes(colour = sapply(.id, colour_map)), ladderize = FALSE, size = 1 ) 
d3 <- p3$data

p_NP_Canu <- ggtree(tree_list[[2]], layout='roundrect', ladderize = FALSE )
dCanu <- p_NP_Canu$data %>%
  mutate(.id = "NP_Canu")

d3_mod_tips <- bind_rows(d3, dCanu) %>%
  filter(!is.na(label))

# cgMLST_plot <- p3 + 
#   geom_tree(data=d3, layout='roundrect', aes(colour = sapply(.id, colour_map)) ) +
#   geom_tree(data=dCanu, layout='roundrect', aes(colour = sapply(.id, colour_map)) ) +
#   geom_tiplab(size=5, linesize=.5, colour = 'black', offset = 20) + 
#   geom_line(aes(x, y, group=label), data=d3_mod_tips, color='grey', linetype = 'dotted') +
#   scale_color_manual(values = c('red', 'black', 'grey'))
# 
# cgMLST_plot

########################################################################
# Reading in the plasmid trees

plasmid_tree_list = lapply(methods, function(method){read_plasmid_mcc_tree(method, 'IncFIB_AP001918')})
class(plasmid_tree_list) <- "multiPhylo"
names(plasmid_tree_list) <- methods

IncI1_tree_list = lapply(methods, function(method){read_plasmid_mcc_tree(method, 'IncI1')})
class(IncI1_tree_list) <- "multiPhylo"
names(IncI1_tree_list) <- methods

###########
# Joint plasmid tree - correct tip labels
plas <- ggtree(plasmid_tree_list[c(1)], layout='roundrect', aes(colour = sapply(.id, colour_map)), ladderize = FALSE ) 
plas3 <- plas$data

p_NP_Canuplas <- ggtree(plasmid_tree_list[[2]], layout='roundrect', ladderize = FALSE )
dCanuplas <- p_NP_Canuplas$data %>%
  mutate(.id = "NP_Canu")

p_PBplas <- ggtree(plasmid_tree_list[[9]], layout='roundrect', ladderize = FALSE )
dPBplas <- p_PBplas$data %>%
  mutate(.id = "PB_Unicycler")

d3_plas_mod_tips <- bind_rows(plas3, dPBplas) %>%
  filter(!is.na(label))

# plasmid_plot <- plas + 
#   geom_tree(data=plas3, layout='roundrect', aes(colour = sapply(.id, colour_map)) ) +
#   #geom_tree(data=dCanu, layout='roundrect', aes(colour = sapply(.id, colour_map)) ) +
#   geom_tree(data=dPBplas, layout='roundrect', aes(colour = sapply(.id, colour_map)) ) +
#   geom_tiplab(size=5, linesize=.5) + 
#   geom_line(aes(x, y, group=label), data=d3_plas_mod_tips, color='black', alpha = 0.7, linetype = 'dotted') +
#   scale_color_manual(values = c('red', 'black', 'grey'))
# 
# plasmid_plot

########################################################################
## reverse x-axis and 
## set offset to create tree on the right hand side of the first tree
cgMLST_plot <- p3 + 
  geom_tree(data=d3, layout='roundrect', aes(colour = sapply(.id, colour_map)), size = 1 ) +
  geom_tree(data=dCanu, layout='roundrect', aes(colour = sapply(.id, colour_map)), size = 1 ) +
  geom_line(aes(x, y, group=label), data=d3_mod_tips, color='black', alpha = 0.7, linetype = 'dotted') +
  scale_color_manual(values = c('red', 'black', 'grey')) +
  geom_treescale(width = 50)

cgMLST_plot

## IncF
Illumina_inv <- plas3
PB_inv <- dPBplas
Canu_inv <- dCanuplas
Illumina_inv$x <- max(plas3$x) - plas3$x + max(cgMLST_plot$data$x) + 150
PB_inv$x <- max(dPBplas$x) - dPBplas$x + max(cgMLST_plot$data$x) + 150
Canu_inv$x <- max(dCanuplas$x) - dCanuplas$x + max(cgMLST_plot$data$x) + 150

# d3_inv_tips <- bind_rows(Illumina_inv, PB_inv, Canu_inv) %>%
#   filter(!is.na(label))

dI_cross <- bind_rows(cgMLST_plot$data, Illumina_inv) %>% 
   filter(!is.na(label))
dP_cross <- bind_rows(cgMLST_plot$data, PB_inv) %>% 
  filter(!is.na(label))
# dC_cross <- bind_rows(cgMLST_plot$data, Canu_inv) %>% 
#   filter(!is.na(label))

## IncI1
IncI1_I <- ggtree(IncI1_tree_list[c(1)], layout='roundrect', ladderize = FALSE ) 
IncI1_C <- ggtree(IncI1_tree_list[[2]], layout='roundrect', ladderize = FALSE )
IncI1_P <- ggtree(IncI1_tree_list[[9]], layout='roundrect', ladderize = FALSE )
IncI1_I$data$x <- max(plas3$x) + max(IncI1_I$data$x) - IncI1_I$data$x + max(cgMLST_plot$data$x) + 300
IncI1_P$data$x <- max(dPBplas$x) + max(IncI1_P$data$x) - IncI1_P$data$x + max(cgMLST_plot$data$x) + 300
IncI1_C$data$x <- max(dCanuplas$x) + max(IncI1_C$data$x) - IncI1_C$data$x + max(cgMLST_plot$data$x) + 300
dP_cross_IncI <- bind_rows(cgMLST_plot$data, IncI1_P$data) %>% 
  filter(!is.na(label))

cgMLST_plot + 
  geom_tree(data=Illumina_inv, layout='roundrect', aes(colour = sapply(.id, colour_map)), size = 1 ) +
  geom_tree(data=PB_inv, layout='roundrect', aes(colour = sapply(.id, colour_map)), size = 1 ) +
  geom_tree(data=Canu_inv, layout='roundrect', aes(colour = sapply(.id, colour_map)), size = 1 ) +
  #geom_line(aes(x, y, group=label), data=d3_inv_tips, color='grey', linetype = 'dotted') +
  #geom_tiplab(size=5, align = TRUE, offset = 20, linesize=.5) + 
  geom_line(aes(x, y, group=label), data=dI_cross, color='red', alpha = 0.7, linetype = 'dashed') +
  #geom_line(aes(x, y, group=label), data=dP_cross, color='grey', linetype = 'longdash') +
  #geom_line(aes(x, y, group=label), data=dC_cross, color='black', linetype = 'dashed') +
  geom_tree(data=IncI1_I$data, layout='roundrect', aes(colour = "Illumina_SPAdes"), size = 1 ) +
  geom_tree(data=IncI1_P$data, layout='roundrect', aes(colour = "Other"), size = 1 ) +
  geom_tree(data=IncI1_C$data, layout='roundrect', aes(colour = "NP_Canu"), size = 1 ) +
  geom_line(aes(x, y, group=label), data=dP_cross_IncI, color='black', alpha = 0.9, linetype = 'dotted') +
  labs(colour = 'Method') +
  theme(legend.position = 'bottom',
        text = element_text(size = 20))

ggsave(paste0(fig_folder, '/Fig1.pdf'), width = 15, height = 10, device="pdf")


###########################################################################
#Add supplementary figure with other IncFIB trees
new_colour_map <- function(x){
  if (x == 'Illumina_SPAdes'){
    return(factor("Illumina", levels = c("Illumina", "NP", "PB")))
  } else if (grepl('NP', x, fixed = TRUE)){
    return(factor("NP", levels = c("Illumina", "NP", "PB")))
  } else {
    return(factor("PB", levels = c("Illumina", "NP", "PB")))
  }
}

plasmid_tree_list = lapply(methods, function(method){read_plasmid_mcc_tree(method, 'IncFIB_AP001918')})
class(plasmid_tree_list) <- "multiPhylo"
names(plasmid_tree_list) <- methods

ggtree(plasmid_tree_list, layout='roundrect', aes(colour = sapply(.id, new_colour_map)), ladderize = FALSE ) +
  facet_wrap(vars(.id), labeller = labeller(.id = method_labels)) +
  geom_treescale(width = 10, x = 0.1) +
  geom_tiplab(size=4, linesize=.5, align = FALSE, show.legend = F) + 
  xlim(0, 60) + 
  scale_color_manual(values = c('red', 'black', 'grey')) +
  labs(colour = 'Method') +
  theme_void() +
  theme(legend.position = 'bottom',
        text = element_text(size = 20))

ggsave(paste0(fig_folder, '/FigS_IncFIB_trees.pdf'), width = 10, height = 15, device="pdf")

##### IncI1

plasmid_tree_list = lapply(methods, function(method){read_plasmid_mcc_tree(method, 'IncI1')})
class(plasmid_tree_list) <- "multiPhylo"
names(plasmid_tree_list) <- methods

ggtree(plasmid_tree_list, layout='roundrect', aes(colour = sapply(.id, new_colour_map)), ladderize = FALSE ) +
  facet_wrap(vars(.id), labeller = labeller(.id = method_labels)) +
  geom_treescale(width = 10, x = 0.1) +
  geom_tiplab(size=4, linesize=.5, align = FALSE, show.legend = F) + 
  xlim(0, 30) + 
  scale_color_manual(values = c('red', 'black', 'grey')) +
  labs(colour = 'Method') +
  theme_void() +
  theme(legend.position = 'bottom',
        text = element_text(size = 20))

ggsave(paste0(fig_folder, '/FigS_IncI1_trees.pdf'), width = 10, height = 15, device="pdf")

