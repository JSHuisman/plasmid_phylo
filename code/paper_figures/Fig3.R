###########################################################
# Compare plasmid trees to chromosomal version.
###########################################################
library('ape')
library('phangorn')
library("HDInterval")
library(ggplot2)
library(tidyverse)

theme_set(theme_light(base_size = 24, base_family = 'sans'))

source('../R/tree_reading_functions.R')

base_dir = '../../data'

###########################################################
methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler', 
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid', 
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods

##########################
# New normed RF plot ####
# These csv files were created similar to the random tree comparison below, 
# but comparing the plasmid tree posterior to the chromosomal MCC tree

RF_df <- read.csv(paste0(base_dir, '/RF_norm_cgMLST_MCC.csv'), stringsAsFactors = FALSE)
#RF_df <- read.csv(paste0(base_dir, '/RF_norm_cgMLST_MCC_PB_Unicycler.csv'), stringsAsFactors = FALSE)
RF_df <- RF_df %>%
  separate(method, into = c('sequence', 'other'), extra = "drop", remove = FALSE) %>%
  mutate(method = factor(RF_df$method, levels = methods)) %>%
  filter(alignment_method == 'roary',
         plasmid %in% c('ColRNAI', "IncI1", "IncFIA", "IncFIB_AP001918"))

RF_plot <- ggplot(RF_df, aes(y=upper.HPD, x=method)) +
  geom_errorbar(aes(ymin=lower.HPD, ymax=upper.HPD, colour = sequence), width=0.2, size=1.2,
                position = position_dodge(width = 0.2), show.legend = F) +
  geom_point(aes(y=mean.RF, colour = sequence), size=3, position = position_dodge(width = 0.2), 
             show.legend = F)  +
  labs(y = 'Robinson-Foulds dist. to cgMLST tree',
       x = "Assembly method") +
  facet_wrap(facets = vars(plasmid), nrow = 2, ncol = 2) + 
  scale_x_discrete(labels= method_labels) +
  scale_color_manual(values = c('red', 'black', 'grey')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15), 
        axis.title=element_text(size=20, face="bold"),
        strip.text = element_text(size=15))

RF_plot
ggsave(paste0("../../figures/RF_norm_cgMLST_MCC.pdf"),
       plot=RF_plot, width = 12, height = 10, device="pdf")


#ggsave(paste0("../../figures/RF_norm_cgMLST_MCC_PB_Unicycler.pdf"),
#        plot=RF_plot, width = 12, height = 10, device="pdf")


###########################################################
# Compare random trees with same tips

RF_df <- data.frame()

included_plasmids <- c('IncI1', 'IncFII_pRSB107', 'IncFIC_FII',
                       'IncFIB_AP001918', 'IncFIA', 'ColRNAI',
                       'Col156', 'Col_MG828')

for (plasmid in included_plasmids){
  # Plasmid tips
  tips = try(get_plasmid_tips(base_dir, plasmid))
  
  for (method in methods){
    # cgMLST files / trees
    cgMLST_MCC_tree <- read_cgmlst_tree(method, base_dir)
    
    if (class(tips)!= 'try-error'){
      # Generate random trees
      rand_trees <- rmtree(100, length(tips$taxa), rooted = TRUE, tip.label = tips$taxa, br = runif,
                           equiprob = FALSE)
      
      # Drop tips from chromosomal
      tips_to_drop = setdiff(cgMLST_MCC_tree$tip.label, tips$taxa)
      red_MCC_tree <- drop_tips_from_tree(cgMLST_MCC_tree, tips_to_drop)
      
      # Calculate RF between plasmid tree and MCC
      RF_to_MCC <- suppressWarnings(RF.dist(rand_trees, tree2 = red_MCC_tree, normalize = T, check.labels = TRUE,
                                            rooted = FALSE))
      HPD = hdi(RF_to_MCC, credMass=0.95)
      
      new_row <- data.frame("min RF" = min(RF_to_MCC),
                            'mean RF' = mean(RF_to_MCC),
                            'max RF' = max(RF_to_MCC),
                            'lower HPD' = HPD[["lower"]],
                            'upper HPD' = HPD[["upper"]],
                            'plasmid' = plasmid,
                            'method' = method,
                            'alignment_method' = 'roary')
      
      RF_df <- rbind(RF_df, new_row)
      
    } 
  }
  
}

write.csv(RF_df, file = paste0(base_dir, '/RF_norm_MCC_rand.csv'))

RF_df <- RF_df %>%
  separate(method, into = c('sequence', 'other'), extra = "drop", remove = FALSE) %>%
  mutate(method = factor(RF_df$method, levels = methods)) %>%
  filter(alignment_method == 'roary',
         plasmid %in% c('ColRNAI', "IncI1", "IncFIA", "IncFIB_AP001918"))

##########################
# New normed RF plot  with random comparison ####

ggplot(NULL, aes(y=upper.HPD, x=method)) +
  geom_errorbar(data = RF_df, aes(ymin=lower.HPD, ymax=upper.HPD, colour = sequence), width=0.2, size=1.2,
                position = position_dodge(width = 0.2), show.legend = F) +
  geom_point(data = RF_df, aes(y=mean.RF, colour = sequence), size=3, position = position_dodge(width = 0.2), 
             show.legend = F)  +
  labs(y = 'Robinson-Foulds dist. to cgMLST tree',
       x = "Assembly method") +
  facet_wrap(facets = vars(plasmid), nrow = 2, ncol = 2) + 
  scale_x_discrete(labels= method_labels) +
  scale_color_manual(values = c('red', 'black', 'grey')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=15), 
        axis.title=element_text(size=20, face="bold"),
        strip.text = element_text(size=15))

ggsave(paste0("../../figures/RF_norm_MCC_rand.pdf"),
       width = 12, height = 10, device="pdf")






