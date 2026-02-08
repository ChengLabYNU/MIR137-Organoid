# 免疫荧光统计
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)

setwd('E:/007-MIR137/4.forbrain_organoid/000-experimental_data/001.Statistical_results/04.WT_Het_ventral_Neuronal_differentiation')

Day25_NKX21 = read.csv('select_slices/day25/Area_Results.csv') %>% 
  mutate(genotype = strsplit(File,'-') %>% sapply('[',1),
         time = 'Day50')

Ventral_NKX21_data = Day25_NKX21 %>%
  select(c('File', 'DAPI_Area','RedandDAPI_Area','genotype','time'))
df_colnames = c('File','DAPI','NKX2.1_DAPI','genotype','time')
colnames(Ventral_NKX21_data) = df_colnames
Ventral_NKX21_data$genotype = factor(Ventral_NKX21_data$genotype,levels = c('WT','Het'))

aggregate(Ventral_NKX21_data$genotype  ,list(Ventral_NKX21_data$genotype  ,Ventral_NKX21_data$time),length)


color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')


ventral_nkx21_plot = ggplot(Ventral_NKX21_data, 
                         aes(x= genotype,y = NKX2.1_DAPI/DAPI))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 1  )+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.8)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression(" NKX2.1"^"+" * "/DAPI area"),
       title = "Day25 Ventral forebrain organoid")+
  theme_classic()+
  scale_y_continuous(labels = scales::label_percent())+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('ventral_nkx21_plot_25.pdf',ventral_nkx21_plot,width = 2.75,height = 3)



Day50_NKX21 = read.csv('select_slices/day50/Area_Results.csv') %>% 
  mutate(genotype = strsplit(File,'-') %>% sapply('[',1),
         time = 'Day50')

Day50_NKX21$genotype = factor(Day50_NKX21$genotype,levels = c('WT','Het'))
df_colnames = c('File','DAPI','NKX2.1_DAPI','genotype','time')
colnames(Day50_NKX21) = df_colnames

aggregate(Day50_NKX21$genotype  ,list(Day50_NKX21$genotype  ,Day50_NKX21$time),length)


color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')


ventral_nkx21_plot_50 = ggplot(Day50_NKX21, 
                            aes(x= genotype,y = NKX2.1_DAPI/DAPI))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 1  )+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.8)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression(" NKX2.1"^"+" * "/DAPI area"),
       title = "Day50 Ventral forebrain organoid")+
  theme_classic()+
  scale_y_continuous(labels = scales::label_percent())+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('ventral_nkx21_plot_50.pdf',ventral_nkx21_plot_50,width = 2.75,height = 3)


