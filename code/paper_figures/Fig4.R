### Parsimony comparisons on trees ########################
#   Author: J.S. Huisman
#   Date: June 2019
#   Last Edit: March 2022

### Packages + General Settings ###########################

library(ape)
library("multicool") #allPerm
library(phytools) # needed for random tree - force ultrametric
library(phangorn) #as.phyDat

library(tidyverse)
library(ggplot2)

###########################################################
theme_set(theme_light(base_size = 24, base_family = 'sans'))

methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler', 
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid', 
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods
###########################################################
base_dir = '../../data'

source('../R/tree_reading_functions.R')
source('../R/parsimony_functions.R')

### Plasmids on cgMLST ####################################

plasmids_to_keep<- read_csv(paste0(base_dir,'/plasmids_to_keep.csv'))
plasmids_removed <- c('Col_MGD2', 'Col8282', 'IncB', 'IncX1', 'IncX3', 'p0111')
plasmids_to_keep <- plasmids_to_keep %>%
  filter(! (name %in% plasmids_removed))

#pars_dfs <- compute_pars_dfs(methods, plasmids_to_keep$name)
# write.csv(pars_dfs$perm_parsimony, file = paste0(base_dir, '/Parsimony_plasmids_on_cgMLST_perm.csv'), row.names = FALSE)
# write.csv(pars_dfs$real_parsimony, file = paste0(base_dir, '/Parsimony_plasmids_on_cgMLST_real.csv'), row.names = FALSE)

perm_parsimony <- read.csv(paste0(base_dir, '/Parsimony_plasmids_on_cgMLST_perm.csv'), 
                     stringsAsFactors = FALSE)
real_parsimony <- read.csv(paste0(base_dir, '/Parsimony_plasmids_on_cgMLST_real.csv'), 
                     stringsAsFactors = FALSE)
pars_dfs <- list('perm_parsimony' = perm_parsimony, 'real_parsimony' = real_parsimony)

# Plasmids on cgMLST: Distribution Plots ------------------

plot_pars_densities(pars_dfs, set_ncol = 4)
ggsave(paste0("../../figures/Parsimony_plasmids_on_cgMLST.pdf"),
       plot=last_plot(), width = 19, height = 14, device="pdf")

# Plasmids on cgMLST: P value -----------------------------

plas_pars_pvalue_df <- compute_pars_pvalue_df(pars_dfs, plasmids_to_keep$name, methods)  %>%
  separate(method, into = c('seq_method', 'align_method'), sep = '_', extra = 'merge', remove = F)

plasmid_plot <- ggplot(data=plas_pars_pvalue_df, aes(y = p_value, x = gene) ) + 
  geom_boxplot() +
  geom_point(aes(y = p_value, colour = seq_method, fill = align_method), shape = 21, 
             show.legend = TRUE, alpha = 1, size = 6, stroke = 3) +
  geom_hline(mapping = aes(yintercept = 0.05), size = 1) +
  scale_colour_manual(values = c('red', 'black', 'grey', 'white')) +  
  scale_fill_manual(values = c('red', 'black', 'grey', 'white')) +  
  facet_wrap(vars(gene), nrow = 1, scales = 'free_x', 
             labeller = labeller(gene = c("Col_MG828" = "ColMG828", 
                                          "Col156"="Col156", "ColRNAI"="ColRNAI", 
                                          "IncFIA" ="IncFIA", "IncFIB_AP001918" = "IncFIB",
                                          "IncFIC_FII"="IncFIC_FII",
                                          "IncFII_pRSB107"= "IncFII", "IncI1"="IncI1") )) +
  labs(x = 'Plasmid', y = 'P-value', colour = 'Sequencing method', 
       fill = 'Assembly method') +
  theme_classic() + 
  theme(axis.text= element_text(size=20), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=24, face="bold"),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 25),
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) +
  guides(colour = guide_legend(nrow = 2, override.aes=list(size = 5) ),
         fill = guide_legend(nrow = 2, override.aes=list(colour = NA, stroke = 2) ))

plasmid_plot

ggsave(paste0("../../figures/Pval_plasmid_parsimony.pdf"),
       plot=last_plot(), width = 17, height = 10, device="pdf")  

################## new depiction

plasmid_plot <- ggplot(data=plas_pars_pvalue_df, aes(y = p_value, x = method) ) + 
  #geom_boxplot(aes(group = gene)) +
  geom_point(aes(y = p_value, colour = seq_method, fill = align_method), shape = 21, 
             show.legend = TRUE, alpha = 1, size = 6, stroke = 3) +
  geom_hline(mapping = aes(yintercept = 0.05), size = 1) +
  scale_colour_manual(values = c('red', 'black', 'grey', 'white')) +  
  scale_fill_manual(values = c('red', 'black', 'grey', 'white')) +  
  facet_wrap(vars(gene), nrow = 2, scales = 'free_x', 
             labeller = labeller(gene = c("Col_MG828" = "ColMG828", 
                                          "Col156"="Col156", "ColRNAI"="ColRNAI", 
                                          "IncFIA" ="IncFIA", "IncFIB_AP001918" = "IncFIB",
                                          "IncFIC_FII"="IncFIC_FII",
                                          "IncFII_pRSB107"= "IncFII", "IncI1"="IncI1") )) +
  labs(x = 'Plasmid', y = 'P-value', colour = 'Sequencing method', 
       fill = 'Assembly method') +
  theme_classic() + 
  theme(axis.text= element_text(size=20), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=24, face="bold"),
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 25),
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) +
  guides(colour = guide_legend(nrow = 2, override.aes=list(size = 5) ),
         fill = guide_legend(nrow = 2, override.aes=list(colour = NA, stroke = 2) ))

plasmid_plot

ggsave(paste0("../../figures/Pval_plasmid_parsimony.pdf"),
       plot=last_plot(), width = 12, height = 12, device="pdf")  



signif_plasmid_plot <- ggplot(plas_pars_pvalue_df, aes(y = gene, x = method)) + 
  geom_tile(aes(fill = (p_value < 0.05))) +
  scale_fill_manual(values = c('grey', 'black')) + 
  scale_x_discrete(labels = method_labels) +
  scale_y_discrete(labels = c("Col_MG828" = "ColMG828", 
                                     "Col156"="Col156", "ColRNAI"="ColRNAI", 
                                     "IncFIA" ="IncFIA", "IncFIB_AP001918" = "IncFIB",
                                     "IncFIC_FII"="IncFIC_FII",
                                     "IncFII_pRSB107"= "IncFII", "IncI1"="IncI1") ) +
  theme_classic() + 
  labs(x = 'Method', y = 'Plasmid', 
       fill = 'Significant association') +
  theme(axis.text = element_text(size=20), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_text(size=24, face="bold"),
        #panel.spacing = unit(2, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) 

signif_plasmid_plot

### Resgenes on cgMLST ####################################
tip_data <- read_tip_data('PB_Unicycler')

genes <- sort(colSums(tip_data[,12:74]))
res_genes <- genes[setdiff(seq_along(genes), grep('Col|Inc|p0111', names(genes)))]

# sel_res_genes <- c("aadA5", "aph.3....Ib", "dfrA17", "aac.6...Ib.cr","blaOXA.1",
#                    "blaTEM.1B","catB3", "sul2", "tet.34.", "blaCTX.M.15",
#                    "mph.A.", "ant.3....Ia", "blaCTX.M.1", "blaCTX.M.27", "sul1")
spec_res_genes <- c("aadA5", "blaCTX.M.1", "blaCTX.M.27","blaCTX.M.15","blaTEM.1B", "sul1", "sul2")
sorted_spec_res_genes <- names(sort(res_genes[spec_res_genes]))

#pars_dfs <- compute_multitree_pars_dfs(methods, sel_res_genes, c("IncI1", "IncFIA", "IncFII_pRSB107"), cgmlst = TRUE)
# write.csv(pars_dfs$perm_parsimony, file = paste0(base_dir, '/Parsimony_resgenes_on_cgMLST_perm.csv'), row.names = FALSE)
# write.csv(pars_dfs$real_parsimony, file = paste0(base_dir, '/Parsimony_resgenes_on_cgMLST_real.csv'), row.names = FALSE)

perm_parsimony <- read.csv(paste0(base_dir, '/Parsimony_resgenes_on_cgMLST_perm.csv'), 
                     stringsAsFactors = FALSE)
real_parsimony <- read.csv(paste0(base_dir, '/Parsimony_resgenes_on_cgMLST_real.csv'), 
                     stringsAsFactors = FALSE)
pars_dfs <- list('perm_parsimony' = perm_parsimony, 'real_parsimony' = real_parsimony)


# Resgenes on cgMLST: Distribution Plots ------------------

subset_pars_dfs <- filter_pars_dfs(pars_dfs, "tree == 'cgMLST'")

plot_pars_densities(subset_pars_dfs, set_ncol = 5, 
                    max_pars = floor(dim(tip_data)[1]/2),
                    adjust = 3, dot_dist = 0.2)

ggsave(paste0("../../figures/Parsimony_resgenes_on_cgMLST.pdf"),
       plot=last_plot(), width = 19, height = 14, device="pdf")

# Resgenes on plasmids / cgMLST - large overview ----------

subset_pars_dfs <- filter_pars_dfs(pars_dfs, 'gene %in% sorted_spec_res_genes')
plot_resgene_scan(subset_pars_dfs, sorted_spec_res_genes, x_var = 'norm_pars',
                  method_labels = method_labels)

ggsave(paste0("../../figures/Parsimony_resgenes_dist_overview.pdf"),
       plot=last_plot(), width = 15, height = 12, device="pdf")

# P value - Resgenes on plasmids / cgMLST -----------------

pars_pvalue_df <- compute_pars_pvalue_df(pars_dfs, sorted_spec_res_genes, methods)

# Violin-PLOT
res_pars_pvalue_df <- pars_pvalue_df %>%
  separate(method, into = c('seq_method', 'align_method'), sep = '_', extra = 'merge', remove = F) %>%
  mutate(tree_recode = recode(tree,
                              "cgMLST" = "Canu", #"Illumina",
                              "IncI1" = "SPAdes", #"NP",
                              "IncFIA" = "Unicycler", #"PB",
                              "IncFII_pRSB107" = "Canu_Hybrid_Polished"))


res_plot <- ggplot(data=res_pars_pvalue_df, aes(y = p_value, x = tree, 
                                group = interaction(gene, tree))) + 
  #geom_violin(aes(fill = tree_recode), scale = 'width', width = 0.5, show.legend = F) +
  geom_point(aes(colour = seq_method, fill = align_method), shape = 21, 
             show.legend = TRUE, alpha = 1, size = 6, stroke = 3) +
  facet_wrap(vars(gene), nrow = 1) +
  scale_colour_manual(values = c('red', 'black', 'grey', 'white')) +  
  scale_fill_manual(values = c('red', 'black', 'grey', 'white')) +  
  geom_hline(mapping = aes(yintercept = 0.05), size = 1) +
  theme_classic() + 
  labs(x = 'Resistance gene', y = 'P-value', colour = 'Sequencing method', 
       fill = 'Assembly method', shape = 'Tree') +
  theme(axis.text= element_text(size=20), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=24, face="bold"),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) 

res_plot

ggsave(paste0("../../figures/Pval_resgene_parsimony_boxplot.pdf"),
       plot=last_plot(), width = 17, height = 10, device="pdf")  

### New depictions ####################################

#spec_res_genes <- c("aadA5", "blaCTX.M.1", "blaCTX.M.27","blaCTX.M.15","blaTEM.1B", "sul1", "sul2")

res_plot <- ggplot(data=res_pars_pvalue_df %>% filter(gene %in% c("blaCTX.M.1", "blaCTX.M.27","sul1", "sul2")), 
       aes(y = p_value, x = method)) + 
  geom_point(aes(colour = seq_method, fill = align_method), shape = 21, 
             show.legend = FALSE, alpha = 1, size = 6, stroke = 3) +
  facet_grid(rows = vars(gene), cols = vars(tree)) +
  scale_colour_manual(values = c('red', 'black', 'grey', 'white')) +  
  scale_fill_manual(values = c('red', 'black', 'grey', 'white')) +  
  geom_hline(mapping = aes(yintercept = 0.05), size = 1) +
  theme_classic() + 
  labs(x = 'Resistance gene', y = 'P-value', colour = 'Sequencing method', 
       fill = 'Assembly method', shape = 'Tree') +
  theme(axis.text= element_text(size=20), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title = element_text(size=24, face="bold"),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) 

res_plot

signif_res_plot <- ggplot(res_pars_pvalue_df, aes(y = gene, x = method)) + 
  geom_tile(aes(fill = (p_value < 0.05)), show.legend = F) +
  facet_wrap(vars(tree), nrow = 1) +
  scale_fill_manual(values = c('grey', 'black')) + 
  scale_x_discrete(labels = method_labels) +
  theme_classic() + 
  labs(x = 'Method', y = 'Resistance gene', 
       fill = 'Significant association') +
  theme(axis.text = element_text(size=20), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_text(size=24, face="bold"),
        panel.spacing = unit(2, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=20)) 

signif_res_plot

### Combine both ####################################
library(patchwork)

plasmid_plot + res_plot +
plot_layout(ncol = 1, heights =c(1, 2), guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30),
        legend.position = 'bottom')

ggsave("../../figures/Pval_plot.pdf", 
       height = 20, width = 16)

##
signif_plasmid_plot + signif_res_plot +
  plot_layout(nrow = 1, widths =c(1, 4), guides = "collect") + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 30),
        legend.position = 'bottom',
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

ggsave("../../figures/Pval_signif_plot.pdf", 
       height = 13, width = 20)
  

