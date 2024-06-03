# Age of SNPs - summarizing estimates from GEVA output
# 6.3.2024
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

# Libraries
library(data.table)
library(tidyverse)
library(cowplot)
library(viridis)
library(patchwork)

# Working directory
setwd("/project/berglandlab/connor/")

# SNP metadata
tot <- readRDS(file="/project/berglandlab/connor/data/classified_snps_filt.rds") %>% 
  mutate(chrom=as.integer(tstrsplit(chrom, "_")[[2]]))

# Betascan results
beta <- readRDS(file="/project/berglandlab/connor/data/betascan.daphnia.TSP.rds") %>% 
  mutate(chrom=as.integer(tstrsplit(chrom, "_")[[2]]))

# Save output
ds <- data.table(rbind(fread("for_daniel/nam.samps.TMRCA.txt"),
                 fread("for_daniel/euro.samps.TMRCA.txt")))


names(ds) = c("chrom","position","ID","ClockModel","Filter","N_Concordant",
              "N_Discordant", "PMean","PMode","PMedian","Set")

# Metadata
ds1 <- na.omit(data.table(ds %>% 
                left_join(tot, by=c("chrom", "position"))) %>% 
  mutate(classified1=case_when(classified=="shared_poly"~"TSP",
                               TRUE ~ "Not TSP")))

# Gene-wise average
ds2 <- data.table(ds1 %>% 
            group_by(classified1, gene, Set) %>% 
            summarize(mean=mean(PMean, na.rm = T),
                      sd=sd(PMean, na.rm = T),
                      n=n()) %>% 
            mutate(se=sd/sqrt(n)))

# Get summary stats
ds2 %>% 
  ggplot(.) +
  geom_pointrange(aes(x=classified1, 
                      y=mean, 
                      ymin=mean-2*se, 
                      ymax=mean+2*se, 
                      color=classified1)) +
  facet_wrap(~Set, nrow = 1)

# Plot age estimates
plot <- {
  
  ds1 %>% 
  ggplot(., 
         aes(x=simpleAnnot, 
             y=PMean,
             fill=classified1)) +
  geom_boxplot() +
  scale_y_log10() +
  coord_flip() +
  facet_wrap(~Set, nrow=1) +
  scale_fill_manual(values = c("TSP"="darkcyan",
                               "Not TSP"="gray")) +
  theme_bw() +
  labs(x = "SNP classification", 
       fill = "",
       y = "Mean TMRCA") +
  theme(legend.title = element_text(face="bold", size=18), 
        strip.text = element_text(face="bold", size=18), 
        axis.text.x = element_text(face="bold", size=18),
        axis.text.y = element_text(face="bold", size=18),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.title = element_text(face="bold", size=20))
  
}

# Save output
ggsave("candgene/geva.euro.nam.tmrca.pdf", plot)

# Testing differences between shared polymorphisms and control SNPs
t.test(ds1[classified1=="TSP"][simpleAnnot=="NS"][!Set %like% "nam"]$PMean,
       ds1[classified1=="Not TSP"][simpleAnnot=="NS"][!Set %like% "nam"]$PMean)

t.test(ds1[classified1=="TSP"][simpleAnnot=="NS"][!Set %like% "nam"]$PMean,
       ds1[classified1=="Not TSP"][simpleAnnot=="NS"][!Set %like% "nam"]$PMean)

# Add intervals - average of BP windows
upp=10
dt2 <- data.table(ds1[gene.x==gene.y] %>% 
                    #group_by(.id, SNP_A) %>% 
                    mutate(distc = as.factor(cut(dist, 
                                                 breaks = seq(from = 0,
                                                              to = max(dist) + upp,
                                                              by = upp), right = F))))

# Add Gene x Gene comps
dt3 <- data.table(dt2 %>% 
                    mutate(gene_comp = case_when(classified.x=="shared_poly" & 
                                                   classified.y=="shared_poly" ~ "TSP",
                                                 TRUE ~ "Not TSP")))

# Summarize
dt.fin <- data.table(dt3 %>% 
                       group_by(dist, cont, gene_comp) %>% 
                       summarise(mean = mean(R2, na.rm = T),
                                 uci.r2 = quantile(R2, probs = 0.95, na.rm = T),
                                 lci.r2 = quantile(R2, probs = 0.05, na.rm = T),
                                 median = median(R2, na.rm = T)))


{
  
  ds1 %>% 
    ggplot(., 
           aes(x=PMean, 
               y=Beta)) +
    geom_point() +
    scale_y_log10() +
    coord_flip() +
    scale_fill_manual(values = c("TSP"="darkcyan",
                                 "Not TSP"="gray")) +
    theme_bw() +
    labs(x = "SNP classification", 
         fill = "",
         y = "Mean TMRCA") +
    theme(legend.title = element_text(face="bold", size=18), 
          strip.text = element_text(face="bold", size=18), 
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=18),
          axis.title.y = element_text(face="bold", size=18),
          axis.title = element_text(face="bold", size=20))
  
}
