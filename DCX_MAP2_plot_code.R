library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggsignif)

setwd('E:/007-MIR137/4.forbrain_organoid/012-IF_results/05.DCX_MAP2/')
DCX_MAP2_data = read.csv('day25/Area_Results.csv',header = T) %>%
  mutate('genotype' = strsplit(File,'-') %>% sapply('[',1) )

DCX_MAP2_data$genotype = factor(DCX_MAP2_data$genotype,levels = c('WT','Het'))

color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')

Day25_DCX_plot = ggplot(DCX_MAP2_data, aes(x= genotype,y = (AF555_Area/DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% DCX"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_DCX_plot.png',Day25_DCX_plot,width = 2.8,height = 3.5)
ggsave('Day25_DCX_plot.pdf',Day25_DCX_plot,width = 2.8,height = 3.5)

Day25_MAP2_plot = ggplot(DCX_MAP2_data, aes(x= genotype,y = (AF488_Area /DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% MAP2"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_MAP2_plot.png',Day25_MAP2_plot,width = 2.8,height = 3.5)
ggsave('Day25_MAP2_plot.pdf',Day25_MAP2_plot,width = 2.8,height = 3.5)


Day25_DCX_MAP2_plot = ggplot(DCX_MAP2_data, aes(x= genotype,y = (AF555_AF488_DAPI_Area/DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% DCX"^"+" *" MAP2"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_DCX_MAP2_plot.png',Day25_DCX_MAP2_plot,width = 2.8,height = 3.5)
ggsave('Day25_DCX_MAP2_plot.pdf',Day25_DCX_MAP2_plot,width = 2.8,height = 3.5)


################################################################################


D50_DCX_MAP2_data = read.csv('day50/Area_Results.csv',header = T) %>%
  mutate('genotype' = strsplit(File,'-') %>% sapply('[',1) )

D50_DCX_MAP2_data$genotype = factor(D50_DCX_MAP2_data$genotype,levels = c('WT','Het'))

color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')

Day50_DCX_plot = ggplot(D50_DCX_MAP2_data, aes(x= genotype,y = (AF555_Area/DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% DCX"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day50_DCX_plot.png',Day50_DCX_plot,width = 2.8,height = 3.5)
ggsave('Day50_DCX_plot.pdf',Day50_DCX_plot,width = 2.8,height = 3.5)

Day50_MAP2_plot = ggplot(D50_DCX_MAP2_data, aes(x= genotype,y = (AF488_Area /DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% MAP2"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day50_MAP2_plot.png',Day50_MAP2_plot,width = 2.8,height = 3.5)
ggsave('Day50_MAP2_plot.pdf',Day50_MAP2_plot,width = 2.8,height = 3.5)


Day50_DCX_MAP2_plot = ggplot(D50_DCX_MAP2_data, aes(x= genotype,y = (AF555_AF488_DAPI_Area/DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% DCX"^"+" *" MAP2"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 12,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day50_DCX_MAP2_plot.png',Day50_DCX_MAP2_plot,width = 2.8,height = 3.5)
ggsave('Day50_DCX_MAP2_plot.pdf',Day50_DCX_MAP2_plot,width = 2.8,height = 3.5)