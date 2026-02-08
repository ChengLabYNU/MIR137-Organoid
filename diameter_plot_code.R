library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggsignif)

setwd('E:/007-MIR137/4.forbrain_organoid/000-experimental_data/001.Statistical_results/01.WHK_diameter')
diameter_data = read.csv('F_Day2550_DF_diameter.csv',header = T)

diameter_data$genotype = factor(diameter_data$genotype,levels = c('WT','Het','KO'))

color_map = c(WT = '#9f9f9f', Het = '#7ca4c5', KO = '#50b8e7')

Day25_D_plot = ggplot(filter(diameter_data,time == 'Day25',part == 'dfo'), 
                      aes(x= genotype,y = diameter))+
  stat_summary(fun = median,geom = "errorbar", 
    aes(ymax = ..y.., ymin = ..y..),color = "grey",
    width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = "Day25 diameter (μm)",title = 'Dorsal forebrain organoid')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")+
  ggplot(filter(diameter_data,time == 'Day25',part == 'vfo'), 
         aes(x= genotype,y = diameter))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = "",title = 'Ventral forebrain organoid')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_D_plot.png',Day25_D_plot,width = 5.5,height = 3.5)
ggsave('Day25_D_plot.pdf',Day25_D_plot,width = 5.5,height = 3.5)


Day50_D_plot = ggplot(filter(diameter_data,time == 'Day50',part == 'dfo'), 
                      aes(x= genotype,y = diameter))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = "Day50 diameter (μm)",title = 'Dorsal forebrain organoid')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")+
  ggplot(filter(diameter_data,time == 'Day50',part == 'vfo'), 
         aes(x= genotype,y = diameter))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = "",title = 'Ventral forebrain organoid')+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day50_D_plot.png',Day50_D_plot,width = 5.5,height = 3.5)
ggsave('Day50_D_plot.pdf',Day50_D_plot,width = 5.5,height = 3.5)
