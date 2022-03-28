library(tidyverse)

fig_folder <- '../../figures'

################
methods = c('Illumina_SPAdes', 'NP_Canu', 'NP_Canu_Hybrid_Polished', 'NP_SPAdes', 'NP_Unicycler',
            'PB_Canu', 'PB_Canu_Hybrid_Polished', 'PB_SPAdes', 'PB_Unicycler')
method_labels = c('Illumina-SPAdes', 'NP-Canu', 'NP-Canu-Hybrid', 'NP-SPAdes-Hybrid', 'NP-Unicycler-Hybrid',
                  'PB-Canu', 'PB-Canu-Hybrid', 'PB-SPAdes-Hybrid', 'PB-Unicycler-Hybrid')
names(method_labels) <- methods

################
plasmids_to_keep <- read_csv(paste0('../../data','/plasmids_to_keep.csv'))
plasmids_removed <- c('Col_MGD2', 'Col8282', 'IncB', 'IncX1', 'IncX3', 'p0111')
plasmids <- plasmids_to_keep %>%
  filter(! (name %in% plasmids_removed)) %>% pull(name)

################

all_roary_len <- tibble()

for (method in methods){
  for (plasmid in plasmids){
    plasmid_df <- read_delim(paste0("/Volumes/Extreme_SSD/sequences/Summaries/", method,
                      "/plasmid_by_rep/", plasmid, "_roary_length.txt"), delim = ' ', col_names = c('name', 'roary_length'))
    plasmid_df['method'] = method
    plasmid_df['plasmid_rep'] = plasmid
    
    all_roary_len <- bind_rows(all_roary_len, plasmid_df)
    }
}


all_roary_len
################

plasmid_aln <- read.csv(paste0('../../data/all_plasmids.csv'), stringsAsFactors=F)

rel_length <- plasmid_aln %>%
  left_join(all_roary_len, by = c('name', 'plasmid_rep', 'method')) %>%
  filter(plasmid_rep %in% plasmids) %>%
  separate(name, into = c('sample', 'contig')) %>%
  mutate(method_name = method_labels[method],
         roary_length = ifelse(is.na(roary_length), 0, roary_length),
         ratio = roary_length/seq_length)

ggplot(rel_length) +
  geom_boxplot(aes(x = method_name, y = ratio)) +
  geom_point(aes(x = method_name, y = ratio), alpha = 0.5) +
  facet_wrap(vars(plasmid_rep), ncol = 2) +
  labs(y = 'Ratio of annotated to total seq. length', x = 'Method') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing.x = unit(3, 'lines'),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(paste0(fig_folder, '/FigS_rel_length.pdf'), width = 10, height = 14, device="pdf")

# by plasmid (method = colour)
ggplot(rel_length) +
  geom_point(aes(x = seq_length, y = roary_length, colour = method_name)) +
  geom_line(aes(x = seq_length, y = roary_length, colour = method_name)) +
  facet_wrap(vars(plasmid_rep), ncol = 4, scale = 'free') +
  labs(y = 'Length found by Roary', x = 'Sequence Length', colour = 'Method') +
  #coord_trans(x = 'log10', y = 'log10') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.spacing.x = unit(3, 'lines'),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'bottom',
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))

ggsave(paste0(fig_folder, '/FigS_annot_scaling.pdf'), width = 14, height = 10, device="pdf")

# by method (plasmid = colour)
ggplot(rel_length) +
  geom_point(aes(x = seq_length, y = roary_length, colour = plasmid_rep)) +
  geom_line(aes(x = seq_length, y = roary_length, colour = plasmid_rep)) +
  facet_wrap(vars(method), ncol = 3, scale = 'free') +
  labs(y = 'Method', x = 'Sample') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        panel.spacing.x = unit(3, 'lines'),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

###
# average

rel_length %>%
  group_by(method, plasmid_rep) %>%
  summarise(mean = mean(ratio))

