library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggsignif)

setwd('E:/007-MIR137/4.forbrain_organoid/000-experimental_data/001.Statistical_results/03.WHK_Ki67_Edu')
ki67_Edu_data = read.csv('Area_Results.csv',header = T) %>%
  mutate('genotype' = strsplit(File,'-') %>% sapply('[',1) )

ki67_Edu_data$genotype = factor(ki67_Edu_data$genotype,levels = c('WT','Het','KO'))

color_map = c(WT = '#9f9f9f', Het = '#7ca4c5', KO = '#50b8e7')

Day25_ki67_Edu_plot = ggplot(ki67_Edu_data, 
                      aes(x= genotype,y = (GreenandDAPI_Area -RedandGreenandDAPI_Area )/GreenandDAPI_Area ))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 0.8)+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.75)+
  geom_signif(comparisons =  list(c("WT", "Het"),c("Het", "KO"),c("WT", "KO")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression("% Ki67"^"-" *"Edu"^"+" *"/Edu"^"+" *"area"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('Day25_ki67_Edu_plot.png',Day25_ki67_Edu_plot,width = 2.8,height = 3.5)
ggsave('Day25_ki67_Edu_plot.pdf',Day25_ki67_Edu_plot,width = 2.8,height = 3.5)


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
