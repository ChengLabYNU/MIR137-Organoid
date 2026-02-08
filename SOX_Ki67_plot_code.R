library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggsignif)

setwd('E:/007-MIR137/4.forbrain_organoid/012-IF_results/02.WHK_SOX2_Ki67_ipsc')
SOX2_ki67_data = read.csv('00.select_slices/Area_Results.csv',header = T) %>%
  mutate('genotype' = strsplit(File,'-') %>% sapply('[',1) )

SOX2_ki67_data$genotype = factor(SOX2_ki67_data$genotype,levels = c('WT','Het','KO'))

color_map = c(WT = '#9f9f9f', Het = '#7ca4c5', KO = '#50b8e7')

Day25_ki67_plot = ggplot(SOX2_ki67_data, aes(x= genotype,y = (AF488_and_DAPI_Area /DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% Ki67"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_ki67_plot.png',Day25_ki67_plot,width = 2.8,height = 3.5)
ggsave('Day25_ki67_plot.pdf',Day25_ki67_plot,width = 2.8,height = 3.5)

Day25_SOX2_plot = ggplot(SOX2_ki67_data, aes(x= genotype,y = (AF555_and_DAPI_Area /DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% SOX2"^"+" *"/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_SOX2_plot.png',Day25_SOX2_plot,width = 2.8,height = 3.5)
ggsave('Day25_SOX2_plot.pdf',Day25_SOX2_plot,width = 2.8,height = 3.5)


Day25_SOX2_Ki67_plot = ggplot(SOX2_ki67_data, aes(x= genotype,y = (AF555_AF488_DAPI_Area/DAPI_Area)))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% SOX2"^"+" *"Ki67"^"+" * "/DAPI area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_SOX2_Ki67_plot.png',Day25_SOX2_Ki67_plot,width = 2.8,height = 3.5)
ggsave('Day25_SOX2_Ki67_plot.pdf',Day25_SOX2_Ki67_plot,width = 2.8,height = 3.5)
