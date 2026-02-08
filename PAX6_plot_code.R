# 免疫荧光统计
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)

setwd('E:/007-MIR137/4.forbrain_organoid/000-experimental_data/001.Statistical_results/03.WT_Het_dorsal_Neuronal_differentiation')

Day25_WT_pax6 = read.csv('selsct_slices/Day25/WT/Area_Results.csv') %>% 
  mutate(genotype = 'WT',time = 'Day25')
Day25_Het_pax6 = read.csv('selsct_slices/Day25/Het/Area_Results.csv') %>% 
  mutate(genotype = 'Het',time = 'Day25')
Dorsal_pax6_data = rbind(Day25_WT_pax6,Day25_Het_pax6) %>%
  select(c('File', 'DAPI_Area','GreenandDAPI_Area','genotype','time'))
df_colnames = c('File','DAPI','PAX6_DAPI','genotype','time')
colnames(Dorsal_pax6_data) = df_colnames
Dorsal_pax6_data$genotype = factor(Dorsal_pax6_data$genotype,levels = c('WT','Het'))

sample_count = aggregate(Dorsal_pax6_data$genotype  ,list(Dorsal_pax6_data$genotype  ,Dorsal_pax6_data$time),length)


color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')


dorsal_pax_plot = ggplot(Dorsal_pax6_data, 
                           aes(x= genotype,y = PAX6_DAPI/DAPI))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 1  )+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.8)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression(" PAX6"^"+" * "/DAPI area"),
       title = "Day25 Dorsal forebrain organoid")+
  theme_classic()+
  scale_y_continuous(labels = scales::label_percent())+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('dorsal_pax_plot_25.pdf',dorsal_pax_plot,width = 2.75,height = 3)



Day50_pax6 = read.csv('selsct_slices/Day50/Area_Results.csv') %>% 
  mutate(genotype = strsplit(File,'-') %>% sapply('[',1),
         time = 'Day50')

Day50_pax6$genotype = factor(Day50_pax6$genotype,levels = c('WT','Het'))
df_colnames = c('File','DAPI','PAX6_DAPI','genotype','time')
colnames(Day50_pax6) = df_colnames

aggregate(Day50_pax6$genotype  ,list(Day50_pax6$genotype  ,Day50_pax6$time),length)


color_map = c(WT = '#9f9f9f', Het = '#7ca4c5')


dorsal_pax_plot_50 = ggplot(Day50_pax6, 
                         aes(x= genotype,y = PAX6_DAPI/DAPI))+
  stat_summary(fun = median,geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y..),color = "grey",
               width = 0.8,linetype = 1,linewidth = 3,alpha = 1  )+
  geom_quasirandom(aes(color = genotype),size = 2,alpha = 0.8)+
  geom_signif(comparisons =  list(c("WT", "Het")),
              step_increase = 0.1,map_signif_level = T)+
  scale_color_manual(values = color_map)+
  labs(y = expression(" PAX6"^"+" * "/DAPI area"),
       title = "Day50 Dorsal forebrain organoid")+
  theme_classic()+
  scale_y_continuous(labels = scales::label_percent())+
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 10,color = "black",face = "bold"),
        legend.position = "none")
ggsave('dorsal_pax_plot_50.pdf',dorsal_pax_plot_50,width = 2.75,height = 3)


