####################################
# Fig2
#
# Author: J.S. Huisman 
# Last edit: November 2021
####################################

library("ggplot2")
library(tidyverse)
library('RColorBrewer')

####################################

base_folder <- '../../data'
fig_folder <- '../../figures'

plasmid_aln <- read.csv(paste0(base_folder, '/all_plasmids.csv'), stringsAsFactors=F)
plasmid_summary <- read_csv(paste0(base_folder, '/long_plasmid_summary.csv'))

################
methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler',
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid',
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods

###############################################
# New plots for paper

selected_plasmids = c("ColRNAI", "IncFIA", "IncFIB_AP001918", "IncI1")

## Seq lengths
seq_data <- plasmid_aln %>%
  filter(plasmid_rep %in% selected_plasmids)

## Alignment lengths

align_data <- plasmid_summary %>%
  rename(plasmid_rep = name) %>%
  filter(plasmid_rep %in% selected_plasmids) %>%
  mutate(delta_len = roary_avg_len - avg_len)

## Both in one plot

ggplot() +
  geom_violin(data = seq_data, aes(x = method, y = seq_length, group = interaction(plasmid_rep, method))) +
  geom_boxplot(data = align_data, aes(x = method, y = roary_avg_len, group = interaction(plasmid_rep, method)),
               colour = 'red') +
  geom_boxplot(data = align_data, aes(x = method, y = avg_len, group = interaction(plasmid_rep, method)),
               colour = 'grey') +
  facet_wrap(vars(plasmid_rep), scale = 'free_y') +
  coord_trans(y = 'log') +
  scale_x_discrete(labels = method_labels) +
  labs(x = 'Assembly Methods', y = 'Contig Length') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(paste0(fig_folder, '/Fig2.pdf'), width = 14, height = 12, device="pdf")



