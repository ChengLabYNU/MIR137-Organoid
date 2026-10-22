save_path = ""
load('F_data_241123.rda')

library(Seurat)
library(tidyverse)
library(ggsignif)
library(ggbreak)
library(ggh4x)
library(ggpubr)
library(patchwork)
library(stringr)
library(UpSetR)
library(reshape2)
library(pheatmap)
library(ggsci)
library(qvalue)
library(clusterProfiler)
library(org.Mm.eg.db)

cmap = c("OLG" = "#c5e0a5","OEG" = "#7bb241",
         "ASC" = "#da6d6d",
         "NEUR" = "#ffa626",
         "EPC" = "#7fcbc4","CPC" = "#009687",
         "EC" = "#cd92d8","PC" = "#ab46bb","VSMC" = "#8d24aa","Hb-VC" = "#8c6ec2","VLMC" = "#6639b7",
         "MG" = "#80d3f9","MNC" = "#28b6f6","MAC" = "#017dc2","DC" = "#BADEFA","NEUT" = "#9fa7d9","T cells" = "#5b6bbf","B cells" = "#3948ab")

Cmap = c("Oligodendrocyte" = "#c5e0a5","Olfactory ensheathing glia" = "#7bb241",
         "Astrocyte" = "#da6d6d",
         "Neuron" = "#ffa626",
         "Ependymocytes" = "#7fcbc4","Choroid plexus epithelial cells" = "#009687",
         "Endothelial cells" = "#cd92d8","Pericytes" = "#ab46bb","Vascular smooth muscle cells" = "#8d24aa","Hemoglobin-expressing vascular cells" = "#8c6ec2","Vascular and leptomeningeal cells" = "#6639b7",
         "Microglia" = "#80d3f9","Monocyte" = "#28b6f6","Macrophage" = "#017dc2","Dendritic cells" = "#BADEFA", "Neutrophils" = "#9fa7d9","T cells" = "#5b6bbf","B cells" = "#3948ab")

geno_color = c('WT' = "#9f9f9f","Het" = "#7faacd","KO" = "#4dc4f7")
ANM_col = cmap[c("ASC","MG","NEUR")]
gene_set1_color = list(celltype = cmap[c("ASC","NEUR","MG")], genotype = c("Het:WT" = "#9b7ebd", "KO:WT" = "#7ab2d3"))
miR137_target = read.csv('mouse_TargetScan8.0__miR-137-3p.predicted_targets.txt', sep = '\t', header = TRUE)$Target.gene
cell_list = c('ASC','MG','NEUR')
GO_info = read.csv('mouse_mm10_GO.txt', sep = '\t')
GO_info = GO_info[GO_info$GO.domain != "",]

result_plot = function(pd,group){
  pd = filter(pd, Count > 1, pvalue <0.05, GO.domain == group)
  if(dim(pd)[1] > 1){
    pd = pd[order(pd$Count,decreasing = T),]
    f_pd = head(pd,10)
    plot = ggplot(f_pd,aes(x=-log10(pvalue),y=reorder(GO.term.name,-log10(pvalue))))+
             geom_col()+
             geom_text(aes(x=-log10(pvalue),y=reorder(GO.term.name,-log10(pvalue)),
                           label=GeneRatio))+
             ylab("GO term")+
             theme_bw()+
             theme(axis.text = element_text(color = "black"))
    return(plot)
  }else{
    print("GO result <= 1")
  }
}

full_name_tsne = DimPlot(F_data, group.by="CellType",reduction = "tsne",cols = Cmap,raster = F)+
  guides(color = guide_legend(override.aes = list(size=3),ncol=2))+
  theme(plot.title = element_blank())
ggsave(paste0(save_path,'full_name_tsne.png'),full_name_tsne,width = 10,height = 4,dpi = 600)
ggsave(paste0(save_path,'full_name_tsne.pdf'),full_name_tsne,width = 10,height = 4)

DimPlot(F_data, group.by="celltype",split.by = 'genotype',reduction = "tsne",cols = cmap,ncol=3)+
  theme(plot.title = element_blank())
DimPlot(F_data, group.by="celltype", split.by = 'orig.ident', reduction = "tsne",cols = cmap,ncol=3)+
  theme(plot.title = element_blank())

marker_gene = c('Gapdh','Olig2','Npy','Cldn10','Tubb3',
                'Ccdc153','Sostdc1','Cldn5','Kcnj8','Acta2','Alas2',
                'Slc6a13','Tmem119','Plac8','Pf4','H2-Eb1','S100a9','Cd3e','Cd79a')
marker_vln = VlnPlot(F_data, features=marker_gene, stack =T,pt.size=0,
                     group.by='celltype',fill='ident',cols = cmap,
                     flip = T)+NoLegend()+
  theme(axis.title.x = element_blank())
ggsave(paste0(save_path,'marker_gene_vln.png'),marker_vln,width = 6,height = 5,dpi=600)
ggsave(paste0(save_path,'marker_gene_vln.pdf'),marker_vln,width = 6,height = 5)

genotype_cT = aggregate(F_data$orig.ident,list(F_data$genotype,F_data$celltype),length)
colnames(genotype_cT) = c('genotype','celltype','count')
genotype_cT$genotype = factor(genotype_cT$genotype,
                              levels = c('WT','Het','KO'))

celltype_stack = ggplot(genotype_cT,aes(x=genotype,y=count,fill = celltype))+
  geom_bar(stat="identity",position="fill",width = 0.75)+
  scale_fill_manual(values = cmap)+
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format(accuracy= 1))+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_text(color = "black",size=12),
        legend.title = element_blank())
ggsave('celltype_stack.pdf',celltype_stack,width = 3.5,height = 5)
ggsave('celltype_stack.png',celltype_stack,width = 3.5,height = 5,dpi=600)

sample_celltype = F_data@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(genotype = sapply(strsplit(orig.ident,'-'),'[',1))%>%
  select(orig.ident, celltype, count, proportion,genotype)
sample_celltype$genotype = factor(sample_celltype$genotype,
                                  levels = c('WT','Het','KO'))

p = ggplot(sample_celltype,aes(x = genotype, y = proportion,
                               colour = genotype, fill = genotype))+ 
  stat_summary(color=NA,
               fun = mean, geom = "bar",
               position = position_dodge(), width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar",color = "black", 
               position = position_dodge(0.7), width = 0.2)+
  geom_point(aes(group = genotype),size=1,color = "black", 
             position = position_dodge(width = 0.9))+
  stat_compare_means(comparisons = list(c("WT","Het"),c("Het","KO"),c("WT","KO")),
                     method = "t.test",hide.ns = T,
                     paired=F,label = "p.signif",
                     tip.length=0.005,step.increase = 0.005) +
  facet_wrap2(~ celltype, scales="free_x", nrow = 1,strip.position = "bottom",
              strip = strip_themed(background_x = elem_list_rect(fill = cmap,alpha=0.3))) +
  scale_fill_manual(values = geno_color)+
  scale_color_manual(values = geno_color)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5),
                     labels = scales::percent_format(accuracy= 1))+
  scale_y_break(c(0.12, 0.3),expand = F,scales=c(0.3,0.7))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8))

ggsave('genotype_proportion_1.pdf',p,width = 10,height = 3)

feature_plot = VlnPlot(F_data, features = c("nCount_RNA","nFeature_RNA","percent.mt"),group.by = "celltype",pt.size = 0,stack = T,flip = T,fill.by = "ident",cols = cmap)+ theme(legend.position = 'none',axis.title.x = element_blank())
feature_plot
ggsave('feature_plot.pdf',feature_plot,width = 7,height = 3)
ggsave('feature_plot.png',feature_plot,width = 7,height = 3)

gene_summary = data.frame(matrix(nrow = 0,ncol = 2))
colnames(gene_summary) = c('DetectedGenes','DetectedCells')
for(i in unique(F_data$celltype)){
  test_cell = subset(F_data,celltype == i)@assays[["RNA"]]@counts
  test_cell_count = rowSums(test_cell)
  test_cell_gene = rownames(test_cell)[test_cell_count>0]
  gene_summary[nrow(gene_summary)+1,] = c(length(test_cell_gene),ncol(test_cell))
  rownames(gene_summary)[nrow(gene_summary)] = i
}

gene_summary$celltype =factor(rownames(gene_summary),
                               levels = rev(levels(F_data$celltype)))
p1 = ggplot(gene_summary,aes(x = celltype,y=DetectedCells,
                             fill = celltype))+
  geom_col()+
  geom_label(aes(y=5,label = DetectedCells),size= 2.7,
             color='black',fill='white',alpha = 0.5,angle=90)+
  scale_fill_manual(values = cmap)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = 'black',size=8),
        axis.text.x = element_text(size=10,angle = 90),
        legend.position='none')+
  scale_y_log10(expand = c(0,0)) 

p2 = ggplot(gene_summary,aes(x = celltype,fill = celltype))+
  geom_col(aes(y=DetectedGenes))+
  geom_label(aes(y=4000,label = DetectedGenes),size= 2.7,
             color='black',fill='white',alpha = 0.5,angle=90)+
  scale_fill_manual(values = cmap)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = 'black',size=8),
        axis.text.x = element_text(size=10,angle = 90),
        legend.position='none')
ggsave('DetectedCell_Genes.pdf',p1/p2,width = 8,height = 4)
ggsave('DetectedCell_Genes.png',p1/p2,width = 8,height = 4,dpi=600)

DE_path = "DEsingle_res"
DE_list = list.files(DE_path)
DE_list_HW = DE_list[grepl("HW",DE_list)]

DE_HW_data = as.data.frame(matrix(,nrow = 0,ncol = 4))
colnames(DE_HW_data) = c("X","norm_foldChange","pvalue.adj.FDR","cell_type")
for(i in DE_list_HW){
  tmp = read.csv(paste0(DE_path,"/",i))
  cell_type = gsub("_"," ",gsub("_HW_DE_results.csv","",i))
  tmp_data = tmp[,c("X","norm_foldChange","pvalue.adj.FDR")]
  tmp_data$cell_type = cell_type
  DE_HW_data = rbind(DE_HW_data,tmp_data)
}

DE_HW_data$cell_type = factor(DE_HW_data$cell_type,levels = levels(F_data$celltype))

DE_HW_data$LFC = log2(DE_HW_data$norm_foldChange)
DE_HW_data$DEG = ""
DE_HW_data[abs(DE_HW_data$LFC)>0.5 & DE_HW_data$pvalue.adj.FDR<0.05,'DEG'] = T

DE_HW_up = DE_HW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC > 0.5 & pvalue.adj.FDR < 0.05 ))

DE_HW_down = DE_HW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC < -0.5 & pvalue.adj.FDR < 0.05 ))

P_DE_HW = ggplot()+
  geom_jitter(data = filter(DE_HW_data,DEG == ""),
              aes(x=cell_type,y=LFC),color='grey90',width=0.2)+
  geom_jitter(data = filter(DE_HW_data,DEG == T),
              aes(x=cell_type,y=LFC,color=cell_type),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_HW_up,aes(x=cell_type,y=8,label=count),color="darkred")+
  geom_text(data = DE_HW_down,aes(x=cell_type,y=-8,label=count),color="darkblue")+
  scale_color_manual(values = cmap)+
  labs(title = "Het:WT DEGs")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))

ggsave('HW_DE_jitter.pdf',P_DE_HW,height = 4,width = 6)
ggsave('HW_DE_jitter.png',P_DE_HW,height = 4,width = 6,dpi = 600)

DE_HW_DEG_count = merge(DE_HW_up,DE_HW_down,by="cell_type")
colnames(DE_HW_DEG_count) = c("cell_type","up","down")
DE_HW_DEG_count$down = -DE_HW_DEG_count$down
DE_HW_DEG_count = melt(DE_HW_DEG_count,variable.name = "deg")
P_DE_HW_gene_count = ggplot(DE_HW_DEG_count,aes(x= cell_type,y =value,fill = deg))+
  geom_col(color='black')+
  scale_fill_manual(values = c('up'= '#ffc0cb','down' = '#add8e6'))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title = element_blank(),
        legend.position = 'none')
ggsave('P_DE_HW_gene_count.pdf',P_DE_HW_gene_count,height = 2.5,width = 6)
ggsave('P_DE_HW_gene_count.png',P_DE_HW_gene_count,height = 2.5,width = 6,dpi = 600)

DE_HW_data$target = ""
DE_HW_data$target[DE_HW_data$X %in% miR137_target & DE_HW_data$DEG ==T] = T

DE_HW_target = DE_HW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(target == T))

DE_HW_data$target = ""
DE_HW_data$target[DE_HW_data$X %in% miR137_target & DE_HW_data$DEG ==T] = T

P_DE_HW_target = ggplot()+
  geom_jitter(data = filter(DE_HW_data,target == ""),
              aes(x=cell_type,y=LFC),color='grey90',width=0.2)+
  geom_jitter(data = filter(DE_HW_data,target == T),
              aes(x=cell_type,y=LFC,color=cell_type),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_HW_target,aes(x=cell_type,y=8,label=count),color="black")+
  scale_color_manual(values = cmap)+
  labs(title = "Het:WT DEGs")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))
ggsave('P_DE_HW_target_jitter.pdf',P_DE_HW_target,height = 4,width = 6)
ggsave('P_DE_HW_target_jitter.png',P_DE_HW_target,height = 4,width = 6,dpi = 600)

DE_HW_up_t = DE_HW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC > 0.5 & pvalue.adj.FDR < 0.05 & target == T))

DE_HW_down_t = DE_HW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC < -0.5 & pvalue.adj.FDR < 0.05 & target == T))

DE_HW_DEG_count_t = merge(DE_HW_up_t,DE_HW_down_t,by="cell_type")
colnames(DE_HW_DEG_count_t) = c("cell_type","up","down")
DE_HW_DEG_count_t$down = -DE_HW_DEG_count_t$down
DE_HW_DEG_count_t = melt(DE_HW_DEG_count_t,variable.name = "deg")
P_DE_HW_DEG_count_t = ggplot(DE_HW_DEG_count_t,aes(x= cell_type,y =value,fill = deg))+
  geom_col(color='black')+
  scale_fill_manual(values = c('up'= '#ffc0cb','down' = '#add8e6'))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title = element_blank(),
        legend.position = 'none')
ggsave('P_DE_HW_DEG_count_t.pdf',P_DE_HW_DEG_count_t,height = 2.5,width = 6)
ggsave('P_DE_HW_DEG_count_t.png',P_DE_HW_DEG_count_t,height = 2.5,width = 6,dpi = 600)

DE_list_KW = DE_list[grepl("KW",DE_list)]

DE_KW_data = as.data.frame(matrix(,nrow = 0,ncol = 4))
colnames(DE_KW_data) = c("X","norm_foldChange","pvalue.adj.FDR","cell_type")
for(i in DE_list_KW){
  tmp = read.csv(paste0(DE_path,"/",i))
  cell_type = gsub("_"," ",gsub("_KW_DE_results.csv","",i))
  tmp_data = tmp[,c("X","norm_foldChange","pvalue.adj.FDR")]
  tmp_data$cell_type = cell_type
  DE_KW_data = rbind(DE_KW_data,tmp_data)
}

DE_KW_data$LFC = log2(DE_KW_data$norm_foldChange)
DE_KW_data$DEG = ""
DE_KW_data[abs(DE_KW_data$LFC)>0.5 & DE_KW_data$pvalue.adj.FDR<0.05,'DEG'] = T

DE_KW_data$cell_type = factor(DE_KW_data$cell_type,levels = levels(F_data$celltype))

DE_KW_up = DE_KW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC > 0.5 & pvalue.adj.FDR < 0.05 ))

DE_KW_down = DE_KW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC < -0.5 & pvalue.adj.FDR < 0.05 ))

P_DE_KW = ggplot()+
  geom_jitter(data = filter(DE_KW_data,DEG == ""),
              aes(x=cell_type,y=LFC),color='grey90',width=0.2)+
  geom_jitter(data = filter(DE_KW_data,DEG == T),
              aes(x=cell_type,y=LFC,color=cell_type),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_KW_up,aes(x=cell_type,y=8,label=count),color="darkred")+
  geom_text(data = DE_KW_down,aes(x=cell_type,y=-8,label=count),color="darkblue")+
  scale_color_manual(values = cmap)+
  labs(title = "KO:WT DEGs")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))

ggsave('KW_DE_jitter.pdf',P_DE_KW,height = 4,width = 6)
ggsave('KW_DE_jitter.png',P_DE_KW,height = 4,width = 6,dpi = 600)

DE_KW_DEG_count = merge(DE_KW_up,DE_KW_down,by="cell_type")
colnames(DE_KW_DEG_count) = c("cell_type","up","down")
DE_KW_DEG_count$down = -DE_KW_DEG_count$down
DE_KW_DEG_count = melt(DE_KW_DEG_count,variable.name = "deg")
P_DE_KW_gene_count = ggplot(DE_KW_DEG_count,aes(x= cell_type,y =value,fill = deg))+
  geom_col(color='black')+
  scale_fill_manual(values = c('up'= '#ee5250','down' = '#41a5f4'))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title = element_blank(),
        legend.position = 'none')
ggsave('P_DE_KW_gene_count.pdf',P_DE_KW_gene_count,height = 2.5,width = 6)
ggsave('P_DE_KW_gene_count.png',P_DE_KW_gene_count,height = 2.5,width = 6,dpi = 600)

DE_KW_up_t = DE_KW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC > 0.5 & pvalue.adj.FDR < 0.05 & target == T))

DE_KW_down_t = DE_KW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC < -0.5 & pvalue.adj.FDR < 0.05 & target == T))

DE_KW_DEG_count_t = merge(DE_KW_up_t,DE_KW_down_t,by="cell_type")
colnames(DE_KW_DEG_count_t) = c("cell_type","up","down")
DE_KW_DEG_count_t$down = -DE_KW_DEG_count_t$down
DE_KW_DEG_count_t = melt(DE_KW_DEG_count_t,variable.name = "deg")
P_DE_KW_DEG_count_t = ggplot(DE_KW_DEG_count_t,aes(x= cell_type,y =value,fill = deg))+
  geom_col(color='black')+
  scale_fill_manual(values = c('up'= '#ee5250','down' = '#41a5f4'))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title = element_blank(),
        legend.position = 'none')
ggsave('P_DE_KW_DEG_count_t.pdf',P_DE_KW_DEG_count_t,height = 2.5,width = 6)
ggsave('P_DE_KW_DEG_count_t.png',P_DE_KW_DEG_count_t,height = 2.5,width = 6,dpi = 600)

DE_KW_data$target = ""
DE_KW_data$target[DE_KW_data$X %in% miR137_target & DE_KW_data$DEG ==T] = T
DE_KW_target = DE_KW_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(target == T))

P_DE_KW_target = ggplot()+
  geom_jitter(data = filter(DE_KW_data,target == ""),
              aes(x=cell_type,y=LFC),color='grey90',width=0.2)+
  geom_jitter(data = filter(DE_KW_data,target == T),
              aes(x=cell_type,y=LFC,color=cell_type),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_KW_target,aes(x=cell_type,y=8,label=count),color="black")+
  scale_color_manual(values = cmap)+
  labs(title = "KO:WT DEGs")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))

ggsave('P_DE_KW_target_jitter.pdf',P_DE_KW_target,height = 4,width = 6)
ggsave('P_DE_KW_target_jitter.png',P_DE_KW_target,height = 4,width = 6,dpi = 600)

DE_KH_path = "HetKO"
DE_list_KH = list.files(DE_KH_path)[grepl("KH",list.files(DE_KH_path))]

DE_KH_data = as.data.frame(matrix(,nrow = 0,ncol = 4))
colnames(DE_KH_data) = c("X","norm_foldChange","pvalue.adj.FDR","cell_type")
for(i in DE_list_KH){
  tmp = read.csv(paste0(DE_KH_path,"/",i))
  cell_type = gsub("_"," ",gsub("_KH_DE_results.csv","",i))
  tmp_data = tmp[,c("X","norm_foldChange","pvalue.adj.FDR")]
  tmp_data$cell_type = cell_type
  DE_KH_data = rbind(DE_KH_data,tmp_data)
}

DE_KH_data$LFC = log2(DE_KH_data$norm_foldChange)
DE_KH_data$DEG = ""
DE_KH_data[abs(DE_KH_data$LFC)>0.5 & DE_KH_data$pvalue.adj.FDR<0.05,'DEG'] = T

DE_KH_data$cell_type = factor(DE_KH_data$cell_type,levels = levels(F_data$celltype))

DE_KH_up = DE_KH_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC > 0.5 & pvalue.adj.FDR < 0.05 ))

DE_KH_down = DE_KH_data %>%
  group_by(cell_type) %>%
  summarize(count = sum(LFC < -0.5 & pvalue.adj.FDR < 0.05 ))

P_DE_KH = ggplot()+
  geom_jitter(data = filter(DE_KH_data,DEG == ""),
              aes(x=cell_type,y=LFC),color='grey90',width=0.2)+
  geom_jitter(data = filter(DE_KH_data,DEG == T),
              aes(x=cell_type,y=LFC,color=cell_type),width=0.2,show.legend =FALSE)+
  geom_text(data = DE_KH_up,aes(x=cell_type,y=8,label=count),color="darkred")+
  geom_text(data = DE_KH_down,aes(x=cell_type,y=-8,label=count),color="darkblue")+
  scale_color_manual(values = cmap)+
  labs(title = "KO:Het DEGs")+
  xlab("")+
  ylab(expression(log[2]("Fold Change")))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0.5),
        axis.text.x = element_text(color = "black",size = 10,angle = 90),
        axis.text.y = element_text(color = "black",size=8),
        axis.title.y = element_text(color = "black",size = 12))

ggsave('KH_DE_jitter.pdf',P_DE_KH,height = 4,width = 6)
ggsave('KH_DE_jitter.png',P_DE_KH,height = 4,width = 6,dpi = 600)

for(i in levels(F_data$celltype)){
  DE_HW_tmp = filter(DE_HW_data,cell_type==i,abs(log2(norm_foldChange))>0.5,pvalue.adj.FDR<0.05)
  if(nrow(DE_HW_tmp)>0){DE_HW_tmp$group = "HW"}
  DE_KW_tmp = filter(DE_KW_data,cell_type==i,abs(log2(norm_foldChange))>0.5,pvalue.adj.FDR<0.05)
  if(nrow(DE_KW_tmp)>0){DE_KW_tmp$group = "KW"}
  DE_KH_tmp = filter(DE_KH_data,cell_type==i,abs(log2(norm_foldChange))>0.5,pvalue.adj.FDR<0.05)
  if(nrow(DE_KH_tmp)>0){DE_KH_tmp$group = "KH"}
  tmp = rbind(DE_HW_tmp,DE_KW_tmp,DE_KH_tmp)
  tmp_up = filter(tmp,log2(norm_foldChange)>0.5)
  tmp_down = filter(tmp,log2(norm_foldChange)< -0.5)
  tmp_up_genelist = group_by(tmp_up,group) %>%
    summarise(gene_list = list(X))%>%
    deframe()
  names(tmp_up_genelist) = paste0(names(tmp_up_genelist),'_up')
  tmp_down_genelist = group_by(tmp_down,group) %>%
    summarise(gene_list = list(X))%>%
    deframe()
  names(tmp_down_genelist) = paste0(names(tmp_down_genelist),'_down')
  bar_color = c(KH_up = "#E57272",KH_down = "#64B4F6",
                       KW_up = "#EE9999",KW_down = "#90CAF8",
                       HW_up = "#FFCCD2",HW_down = "#BADEFA")
  plot_list = c(tmp_up_genelist,tmp_down_genelist)
  set_order = names(plot_list)
  tmp_plot = upset(fromList(plot_list), sets = set_order,
        nsets = 6,mb.ratio = c(0.4, 0.6),nintersects=NA,shade.color=cmap[i],
      order.by="freq", decreasing = TRUE,keep.order=T,set_size.show	=T,
      sets.bar.color=bar_color[names(plot_list)])
  tmp_plot
  pdf(file = paste0('DEG_upset/',i,'_DEGs_upset.pdf'),width = 5,height = 3,onefile=FALSE)
  print(tmp_plot)
  dev.off()
  png(file = paste0('DEG_upset/',i,'_DEGs_upset.png'),width = 5,height = 3,res=300,units="in")
  print(tmp_plot)
  dev.off()
}

for(i in cell_list){
  DE_HW_tmp = filter(DE_HW_data,cell_type==i,abs(log2(norm_foldChange))>0.5,pvalue.adj.FDR<0.05)
  if(nrow(DE_HW_tmp)>0){DE_HW_tmp$group = "HW"}
  DE_KW_tmp = filter(DE_KW_data,cell_type==i,abs(log2(norm_foldChange))>0.5,pvalue.adj.FDR<0.05)
  if(nrow(DE_KW_tmp)>0){DE_KW_tmp$group = "KW"}
  tmp = rbind(DE_HW_tmp,DE_KW_tmp)
  tmp_up = filter(tmp,log2(norm_foldChange)>0.5)
  tmp_down = filter(tmp,log2(norm_foldChange)< -0.5)
  tmp_up_genelist = group_by(tmp_up,group) %>%
    summarise(gene_list = list(X))%>%
    deframe()
  names(tmp_up_genelist) = paste0(names(tmp_up_genelist),'_up')
  tmp_down_genelist = group_by(tmp_down,group) %>%
    summarise(gene_list = list(X))%>%
    deframe()
  names(tmp_down_genelist) = paste0(names(tmp_down_genelist),'_down')
  bar_color = c(KW_up = "#EE9999",KW_down = "#90CAF8",
                HW_up = "#FFCCD2",HW_down = "#BADEFA")
  plot_list = c(tmp_up_genelist,tmp_down_genelist)
  set_order = names(plot_list)
  tmp_plot = upset(fromList(plot_list), sets = set_order,
        nsets = 6,mb.ratio = c(0.4, 0.6),nintersects=NA,shade.color=cmap[i],
      order.by="freq", decreasing = TRUE,keep.order=T,set_size.show	=T,
      sets.bar.color=bar_color[names(plot_list)])
  tmp_plot
  pdf(file = paste0('DEG_upset/',i,'_HK_DEGs_upset.pdf'),width = 3,height = 3,onefile=FALSE)
  print(tmp_plot)
  dev.off()
  png(file = paste0('DEG_upset/',i,'_HK_DEGs_upset.png'),width = 3,height = 3,res=300,units="in")
  print(tmp_plot)
  dev.off()
}

for(i in cell_list){
  tmp_HW_up = filter(DE_HW_data, cell_type == i,
                       log2(norm_foldChange) > 0.5, pvalue.adj.FDR < 0.05)
  tmp_HW_down = filter(DE_HW_data, cell_type == i,
                       log2(norm_foldChange) < -0.5, pvalue.adj.FDR < 0.05)
  tmp_KW_up = filter(DE_KW_data, cell_type == i,
                       log2(norm_foldChange) > 0.5, pvalue.adj.FDR < 0.05)
  tmp_KW_down = filter(DE_KW_data ,cell_type == i,
                         log2(norm_foldChange)< -0.5, pvalue.adj.FDR < 0.05)
  assign(paste0(i,'_HW_updeg'),tmp_HW_up)
  assign(paste0(i,'_HW_downdeg'),tmp_HW_down)
  assign(paste0(i,'_KW_updeg'),tmp_KW_up)
  assign(paste0(i,'_KW_downdeg'),tmp_KW_down)
}

for(i in cell_list){
  tmp_HW_up = get(paste0(i,'_HW_updeg'))
  tmp_HW_down = get(paste0(i,'_HW_downdeg'))
  tmp_KW_up = get(paste0(i,'_KW_updeg'))
  tmp_KW_down = get(paste0(i,'_KW_downdeg'))
  if(length(rownames(tmp_HW_up)) > 1){
    tmp_HW_up_go = enrich_go_local(tmp_HW_up$X,GO_info)
    assign(paste0(i,'_HW_updeg_go'),tmp_HW_up_go)
    assign(paste0(i,'_HW_updeg_go_bp'),result_plot(tmp_HW_up_go,"biological_process"))
    assign(paste0(i,'_HW_updeg_go_cc'),result_plot(tmp_HW_up_go,"cellular_component"))
    assign(paste0(i,'_HW_updeg_go_mf'),result_plot(tmp_HW_up_go,"molecular_function"))
    write.csv(tmp_HW_up_go,paste0(i,'_HW_updeg_go.csv'))
    }
  if(length(rownames(tmp_HW_down)) > 1){
    tmp_HW_down_go = enrich_go_local(tmp_HW_down$X,GO_info)
    assign(paste0(i,'_HW_downdeg_go'),tmp_HW_down_go)
    assign(paste0(i,'_HW_downdeg_go_bp'),result_plot(tmp_HW_down_go,"biological_process"))
    assign(paste0(i,'_HW_downdeg_go_cc'),result_plot(tmp_HW_down_go,"cellular_component"))
    assign(paste0(i,'_HW_downdeg_go_mf'),result_plot(tmp_HW_down_go,"molecular_function"))
    write.csv(tmp_HW_down_go,paste0(i,'_HW_downdeg_go.csv'))
    }
  if(length(rownames(tmp_KW_up)) > 1){
    tmp_KW_up_go = enrich_go_local(tmp_KW_up$X,GO_info)
    assign(paste0(i,'_KW_updeg_go'),tmp_KW_up_go)
    assign(paste0(i,'_KW_updeg_go_bp'),result_plot(tmp_KW_up_go,"biological_process"))
    assign(paste0(i,'_KW_updeg_go_cc'),result_plot(tmp_KW_up_go,"cellular_component"))
    assign(paste0(i,'_KW_updeg_go_mf'),result_plot(tmp_KW_up_go,"molecular_function"))
    write.csv(tmp_KW_up_go,paste0('DEG_upset/',i,'_KW_updeg_go.csv'))
    }
  if(length(rownames(tmp_KW_down)) > 1){
    tmp_KW_down_go = enrich_go_local(tmp_KW_down$X,GO_info)
    assign(paste0(i,'_KW_downdeg_go'),tmp_KW_down_go)
    assign(paste0(i,'_KW_downdeg_go_bp'),result_plot(tmp_KW_down_go,"biological_process"))
    assign(paste0(i,'_KW_downdeg_go_cc'),result_plot(tmp_KW_down_go,"cellular_component"))
    assign(paste0(i,'_KW_downdeg_go_mf'),result_plot(tmp_KW_down_go,"molecular_function"))
    write.csv(tmp_KW_down_go,paste0('DEG_upset/',i,'_KW_downdeg_go.csv'))
    }
}

plot_list = ls()[grepl('deg_go_',ls())]

for(i in plot_list){
  if(length(class(get(i))) == 2){
    ggsave(paste0('DEG_upset/',i,'.pdf'),get(i),width = 6,height = 5)
  }
}

enrich_go_local = function(gene, GO_info){
  gene = unique(gene)
  Q_GO_term = unique(filter(GO_info,Gene.name %in% gene))
  GO_term_gene = filter(GO_info,GO.term.accession %in% Q_GO_term$GO.term.accession)
  GO_term_gene = unique(GO_term_gene)
  Q_Gene2Term = with(Q_GO_term, split(as.character(GO.term.accession), 
                                      as.character(Gene.name)))
  Q_Term2Gene = with(Q_GO_term, split(as.character(Gene.name), 
                                      as.character(GO.term.accession)))
  Q_Term2GS = with(Q_GO_term, split(as.character(Gene.name), 
                                    as.character(GO.term.accession)))
  all_Term2Gene = with(GO_term_gene, split(as.character(Gene.name), 
                                           as.character(GO.term.accession)))
  range_index = sapply(all_Term2Gene, length) > 10 & sapply(all_Term2Gene, length) < 500
  Q_Term2Gene = Q_Term2Gene[range_index]
  all_Term2Gene = all_Term2Gene[range_index]
  Q_Term2GS = Q_Term2GS[range_index]
  k = sapply(Q_Term2Gene, length)
  k = as.numeric(k)
  M = sapply(all_Term2Gene, length)
  M = as.numeric(M)
  N = as.numeric(rep(length(unique(GO_info$Gene.name)), length(M)))
  n = as.numeric(rep(length(Q_Gene2Term), length(M)))
  args.df = data.frame(numWdrawn = k - 1, numW = M, numB = N - 
                         M, numDrawn = n)
  pvalues = apply(args.df, 1, function(n) phyper(n[1], n[2], 
                                                 n[3], n[4], lower.tail = FALSE))
  GeneRatio = apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], 
                                                                   "/", x[2], sep = "", collapse = ""))
  BgRatio = apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], 
                                                                 "/", x[2], sep = "", collapse = ""))
  FoldChange = apply(data.frame(a = k, b = n ,c = M, d = N), 1, 
                     function(x) (x[1]/x[2])/(x[3]/x[4]))
  Over = data.frame(ID = as.character(names(Q_Term2Gene)), GeneRatio = GeneRatio, 
                    BgRatio = BgRatio, FoldChange = FoldChange, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj = p.adjust(Over$pvalue, method = "BH")
  qobj = tryCatch(qvalue(p = Over$pvalue, lambda = 0.05, 
                         pi0.method = "bootstrap"), error = function(e) NULL)
  qvalues = qobj$qvalues
  geneID = sapply(Q_Term2Gene, function(i) paste(i, collapse = "/"))
  genename = sapply(Q_Term2GS, function(i) paste(i, collapse = "/"))
  Over = data.frame(Over, p.adjust = p.adj, qvalue = qvalues, 
                    geneID = geneID, genename=genename, Count = k, stringsAsFactors = FALSE)
  all_Description = GO_info[,c("GO.term.accession","GO.term.name","GO.term.definition","GO.domain")]
  all_Description = unique(all_Description)
  Over = merge(Over,all_Description,by.x='ID',by.y="GO.term.accession",all.x=T)
  nc = ncol(Over)
  Over = Over[, c(1, nc, 2:(nc - 1))]
  Over = Over[order(pvalues), ]
  Over$ID = as.character(Over$ID)
  rownames(Over) = as.character(Over$ID)
  Over = filter(Over,pvalue<0.05,Count>1)
  return(Over)
}

library(CellChat)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

load('WT_cellchat.rda')
load('Het_cellchat.rda')
load('KO_cellchat.rda')

cmap = c("ASC" = "#da6d6d","NEUR" = "#ffa626","MG" = "#80d3f9")
geno_color = c('WT' = "#9f9f9f","Het" = "#7faacd","KO" = "#4dc4f7")

cellchat = mergeCellChat(list(WT = WT_cellchat, Het = Het_cellchat, KO = KO_cellchat),
                         add.names = c("WT","Het","KO"))
cellchat_WH  = mergeCellChat(list(WT = WT_cellchat, Het = Het_cellchat), add.names = c("WT","Het"))
cellchat_WK  = mergeCellChat(list(WT = WT_cellchat, KO  = KO_cellchat),  add.names = c("WT","KO"))
cellchat_WHK = mergeCellChat(list(WT = WT_cellchat, Het = Het_cellchat, KO = KO_cellchat),
                             add.names = c("WT","Het","KO"))

par(mfrow = c(2,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c('Het','WT'), measure = "count",
                          label.edge = T, color.use = cmap, sources.use = 'ASC',
                          title.name = "Het vs WT - ASC as sender")
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c('KO','WT'),  measure = "count",
                          label.edge = T, color.use = cmap, sources.use = 'ASC',
                          title.name = "KO vs WT - ASC as sender")
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c('Het','WT'), measure = "count",
                          label.edge = T, color.use = cmap, targets.use = 'ASC',
                          title.name = "Het vs WT - ASC as receiver")
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c('KO','WT'),  measure = "count",
                          label.edge = T, color.use = cmap, targets.use = 'ASC',
                          title.name = "KO vs WT - ASC as receiver")

only_diff_net = function(X){
  tmp_data = filter(X$signaling.contribution, pvalues < 0.05, contribution.relative.1 != 1)
  ggplot(tmp_data, aes(x = contribution.scaled, y = name, fill = group))+
    geom_bar(position = "fill", stat = "identity")+
    geom_vline(xintercept = 0.5, linetype = 2)+
    scale_fill_manual(values = geno_color)+
    labs(x = "Relative number of interactions")+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(color = 'black', size = 12),
          axis.text = element_text(color = 'black'))
}

rankNet_HW_sA_p = only_diff_net(rankNet(cellchat_WHK, comparison = c(2,1), mode = "comparison",
                                        stacked = T, do.stat = TRUE, sources.use = "ASC",
                                        measure = 'count', return.data = T))
rankNet_KW_sA_p = only_diff_net(rankNet(cellchat_WHK, comparison = c(3,1), mode = "comparison",
                                        stacked = T, do.stat = TRUE, sources.use = "ASC",
                                        measure = 'count', return.data = T))
rankNet_HW_tA_p = only_diff_net(rankNet(cellchat_WHK, comparison = c(2,1), mode = "comparison",
                                        stacked = T, do.stat = TRUE, targets.use = "ASC",
                                        measure = 'count', return.data = T))
rankNet_KW_tA_p = only_diff_net(rankNet(cellchat_WHK, comparison = c(3,1), mode = "comparison",
                                        stacked = T, do.stat = TRUE, targets.use = "ASC",
                                        measure = 'count', return.data = T))

rankNet_HW_sA_p + rankNet_KW_sA_p + rankNet_HW_tA_p + rankNet_KW_tA_p

pathways.show = "SEMA4"
weight.max = getMaxWeight(cellchat_WHK, slot.name = "netP", attribute = pathways.show)

par(mfrow = c(1,3), xpd = TRUE)
for (i in 1:length(cellchat_WHK)) {
  netVisual_aggregate(cellchat_WHK[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      color.use = cmap,
                      signaling.name = paste(pathways.show, c("WT","Het","KO")[i]))
}

## (I) SEMA4 通路气泡图
netbubble_sema4_h = netVisual_bubble(cellchat_WHK, signaling = "SEMA4", remove.isolate = T,
                                     comparison = c(2,1), sources.use = c(3,4,12), targets.use = c(3,4,12))
netbubble_sema4_k = netVisual_bubble(cellchat_WHK, signaling = "SEMA4", remove.isolate = T,
                                     comparison = c(3,1), sources.use = c(3,4,12), targets.use = c(3,4,12))
netbubble_sema4 = rbind(netbubble_sema4_h$data, netbubble_sema4_k$data)
netbubble_sema4$dataset = factor(netbubble_sema4$dataset, levels = c('WT','Het','KO'))

p_netbubble_sema4 = ggplot(netbubble_sema4, aes(y = group.names, x = dataset))+
  geom_point(aes(colour = prob), size = 5)+
  facet_grid(~interaction_name)+
  scale_color_material('red')+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_text(color = 'black'),
        strip.background = element_rect(fill = NA))
ggsave('p_netbubble_sema4.pdf', p_netbubble_sema4, width = 12, height = 2)
ggsave('p_netbubble_sema4.png', p_netbubble_sema4, width = 3,  height = 2.5)