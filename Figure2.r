library(tidyverse)
library(reshape2)
library(Mfuzz)
library(ggsci)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(patchwork)
library(pheatmap)

trans = read.csv('E:/015-Reference/Mus_musculus.GRCm38.102.biomart.txt',sep='\t')
mir137_target = read.csv('E:/015-Reference/TargetScan8.0__miR-137-3p.Mouse.predicted_targets.txt',sep='\t')
mir137_target$transcript = sapply(strsplit(mir137_target$Representative.transcript,'\\.'),'[',1)


Count =  read.csv('E:/125-MouseEmbryo/_NewReference_/FB_all_time/count/BP_gene_count.tsv',sep='\t')
Count = Count[,!grepl('MB',colnames(Count))]
Count = merge(Count,trans[,c('Gene.stable.ID','Gene.name','Transcript.stable.ID')],
              by.x='gene_id',by.y ='Gene.stable.ID')
Count = Count[!duplicated(Count),]

time_order = c("E10.5.W1","E10.5.W2","E10.5.W3",
               "E10.5.H1","E10.5.H2","E10.5.H3",
               "E10.5.K1","E10.5.K2","E10.5.K3",
               "E11.5.W1","E11.5.W2","E11.5.W3",
               "E11.5.H1","E11.5.H2","E11.5.H3",
               "E11.5.K1","E11.5.K2","E11.5.K3",
               "E12.5.W1","E12.5.W2","E12.5.W3",
               "E12.5.H1","E12.5.H2","E12.5.H3",
               "E12.5.K1","E12.5.K2","E12.5.K3",
               "E13.5.W1","E13.5.W2","E13.5.W3",
               "E13.5.H1","E13.5.H2","E13.5.H3",
               "E13.5.K1","E13.5.K2","E13.5.K3",
               "E14.5.W1","E14.5.W2","E14.5.W3",
               "E14.5.H1","E14.5.H2","E14.5.H3",
               "E14.5.K1","E14.5.K2","E14.5.K3",
               "E15.5.W1","E15.5.W2","E15.5.W3",
               "E15.5.H1","E15.5.H2","E15.5.H3", 
               "E15.5.K1","E15.5.K2","E15.5.K3")

mir137_target_count = filter(Count,Transcript.stable.ID %in% mir137_target$transcript)
mir137_target_count[mir137_target_count$Transcript.stable.ID == 'ENSMUST00000102903','Gene.name'] = 'Akap2'
mir137_target_count[mir137_target_count$Transcript.stable.ID == 'ENSMUST00000150412','Gene.name'] = 'Palm2Akap2'
rownames(mir137_target_count) = mir137_target_count$Gene.name
mir137_target_count = mir137_target_count[,colnames(mir137_target_count) != "gene_id"]
mir137_target_count = mir137_target_count[,colnames(mir137_target_count) != "Gene.name"]
mir137_target_count = mir137_target_count[,colnames(mir137_target_count) != "Transcript.stable.ID"]
mir137_target_count = mir137_target_count[,time_order]

TPM =  read.csv('E:/125-MouseEmbryo/_NewReference_/FB_all_time/count/BP_gene_TPM.tsv',sep='\t')
TPM = merge(TPM,trans[,c('Gene.stable.ID','Gene.name','Transcript.stable.ID')],
              by.x='gene_id',by.y ='Gene.stable.ID')
TPM = TPM[!duplicated(TPM),]
mir137_target_TPM = filter(Count,Transcript.stable.ID %in% mir137_target$transcript)
mir137_target_TPM[mir137_target_TPM$Transcript.stable.ID == 'ENSMUST00000102903','Gene.name'] = 'Akap2'
mir137_target_TPM[mir137_target_TPM$Transcript.stable.ID == 'ENSMUST00000150412','Gene.name'] = 'Palm2Akap2'
rownames(mir137_target_TPM) = mir137_target_TPM$Gene.name
mir137_target_TPM = mir137_target_TPM[,colnames(mir137_target_TPM) != "gene_id"]
mir137_target_TPM = mir137_target_TPM[,colnames(mir137_target_TPM) != "Gene.name"]
mir137_target_TPM = mir137_target_TPM[,colnames(mir137_target_TPM) != "Transcript.stable.ID"]
mir137_target_TPM = mir137_target_TPM[,time_order]

time_color = c(E10.5 = "green", E11.5 = "teal", E12.5 = "cyan",
               E13.5 = "blue", E14.5 ="indigo", E15.5 = "deep-purple")
time_point = c('E10.5','E11.5','E12.5','E13.5','E14.5','E15.5')
for(i in time_point){
  set.seed(999)
  count_tmp = mir137_target_count[,grep(i,colnames(mir137_target_count))]
  count_tmp = data.frame(WT = apply(count_tmp[,grep('W',colnames(count_tmp))], 1, mean),
                         Het = apply(count_tmp[,grep('H',colnames(count_tmp))], 1, mean),
                         KO = apply(count_tmp[,grep('K',colnames(count_tmp))], 1, mean))
  
  mfuzz_data = as.matrix(count_tmp)
  mfuzz_mx = ExpressionSet(assayData = mfuzz_data)
  mfuzz_mx = filter.NA(mfuzz_mx,thres = 0.25)
  mfuzz_mx = filter.std(mfuzz_mx,min.std=0)

  gene.s = standardise(mfuzz_mx)
  m = mestimate(gene.s)
  cl = mfuzz(gene.s, c = 4, m = m)
  
  assign(paste0(i,'_gene.s'),gene.s)
  assign(paste0(i,'_m'),m)
  assign(paste0(i,'_cl'),cl)
  
  save(list = ls(pattern = i), file = paste0(i,"_mfuzz.RData"))
  
  plot_data = exprs(gene.s) %>%
    cbind(cl$cluster, sapply(names(cl$cluster), function(g) cl$membership[g, cl$cluster[g]])) %>%
    `colnames<-`(c(colnames(.)[1:(ncol(.) - 2)], "cluster", "membership")) %>%
    data.frame() %>%
    mutate(gene = rownames(.)) %>%
    melt(measure.vars = 1:3) %>%
    arrange(membership) %>%
    mutate(gene = factor(gene, unique(gene)))
  plot_data$cluster = paste0('Cluster',plot_data$cluster)
  main_line = cl$centers %>%
    as.data.frame() %>%
    mutate(cluster = rownames(.)) %>%
    melt(id.vars = "cluster")
  main_line$cluster = paste0('Cluster',main_line$cluster)
  
  # 绘图并保存
  plot =  ggplot(plot_data, aes(x = variable, y = value, group = gene, color = membership)) +
    geom_line(size = 1) +
    geom_line(data = main_line, aes(x = variable, y = value, group = cluster), color = 'black', size = 1) +
    facet_grid(~cluster) +
    theme_classic() +
    scale_x_discrete(expand = c(0.1,0.1))+
    scale_color_material("blue-grey", values = c(0, 1)) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_blank())
  assign(paste0(i,'_cluster_plot_grey'),plot)
  ggsave(paste0(i,'_cluster_plot_grey.pdf'), plot, width = 12, height = 2)
}

plot_list <- mget(ls(pattern = "_cluster_plot_grey$"))
combined_plot <- Reduce(`/`, plot_list)
ggsave('ALL_cluster_plot_grey.pdf', combined_plot, width = 8, height = 10)


attern1_gene = E12.5_C3_gene[E12.5_C3_gene %in% E14.5_C1_gene]
pattern1_GO = enrichGO(pattern1_gene,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.1)

PD_pattern1_GO = arrange(pattern1_GO@result,p.adjust)[1:10,]
PD_pattern1_GO$Description = factor(PD_pattern1_GO$Description,
                                    levels = rev(PD_pattern1_GO$Description))
P_pattern1_GO = ggplot(PD_pattern1_GO,aes(x=-log10(p.adjust),y=Description))+
  geom_col(fill = "#A5D6A6",width = 0.6)+
  geom_text(x = 0.05,aes(label = Description),hjust = 0,size = 6,size.unit = 'pt')+
  geom_text(x = 0.05,aes(label = geneID),color =  "#2D7D32",
            hjust = 0,vjust = 2.3,size = 6,fontface ='italic',size.unit = 'pt')+
  scale_x_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.text.y = element_blank())
ggsave('pattern1_2_GO.pdf',P_pattern1_GO,width = 6,height = 6, units = 'cm')


### pattern2
# patten2 het和KO都上调的gene
# 包含E10.5 Cluster1,E11.5 Cluster1和E12.5 Cluster1


pattern2_gene = names(which(table(c(E10.5_C1_gene, E11.5_C1_gene, E12.5_C1_gene)) > 1))
pattern2_GO = enrichGO(pattern2_gene,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.1)

PD_pattern2_GO = arrange(pattern2_GO@result,p.adjust)[1:10,]
PD_pattern2_GO$Description = factor(PD_pattern2_GO$Description,
                                    levels = rev(PD_pattern2_GO$Description))
P_pattern2_GO = ggplot(PD_pattern2_GO,aes(x=-log10(p.adjust),y=Description))+
  geom_col(fill = "#7FDEEA",width = 0.6)+
  geom_text(x = 0.05,aes(label = Description),hjust = 0,size = 6,size.unit = 'pt')+
  geom_text(x = 0.05,aes(label = geneID),color =  "#00838E",
            hjust = 0,vjust = 2.3,size = 6,fontface ='italic',size = 6,size.unit = 'pt')+
  scale_x_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"),
        axis.text.y = element_blank())
ggsave('pattern2_2_GO.pdf',P_pattern2_GO,width = 6,height = 6, units = 'cm')


### pattern3
pattern3_gene = names(which(table(c(E10.5_C3_gene, E11.5_C4_gene, E14.5_C1_gene)) > 1))
pattern3_GO = enrichGO(pattern3_gene,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",
                       ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.1)

PD_pattern3_GO = arrange(pattern3_GO@result,p.adjust)[1:10,]
PD_pattern3_GO$Description = factor(PD_pattern3_GO$Description,
                                    levels = rev(PD_pattern3_GO$Description))
P_pattern3_GO = ggplot(PD_pattern3_GO,aes(x=-log10(p.adjust),y=Description))+
  geom_col(fill = "#9FA7D9",width = 0.6)+
  geom_text(x = 0.05,aes(label = Description),hjust = 0,size = 6,size.unit = 'pt')+
  geom_text(x = 0.05,aes(label = geneID),color =  "#273492",
            hjust = 0,vjust = 2.3,size = 6,fontface ='italic',size = 6,size.unit = 'pt')+
  scale_x_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(colour = "black"),
        axis.text.y = element_blank())
ggsave('pattern3_2_GO.pdf',P_pattern3_GO,width = 6,height = 6, units = 'cm')

pattern1_venn = venn.diagram(list("E12.5 Cluster3" = E12.5_C3_gene,"E14.5 Cluster1" = E14.5_C1_gene), 
             filename = NULL, fill = c('#ADD480', '#548B2E'), 
             alpha = c(0.5), cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'sans', 
             cat.cex = 1.5, cat.fontfamily = 'sans')
ggsave(pattern1_venn, file="pattern1_venn.svg", device = "svg")

pattern2_venn = venn.diagram(list("E10.5 Cluster1" = E10.5_C1_gene,
                                  "E11.5 Cluster1" = E11.5_C1_gene,
                                  "E12.5 Cluster1" = E12.5_C1_gene), 
             filename = NULL, fill = c('#7FDEEA', '#00BBD3',"#00838E"), 
             alpha = c(0.5), cat.col = rep('black', 3), 
             col = 'black', cex = 1.5, fontfamily = 'sans', 
             cat.cex = 1.5, cat.fontfamily = 'sans')
ggsave(pattern2_venn, file="pattern2_venn.svg", device = "svg")

pattern3_venn = venn.diagram(list("E10.5 Cluster3" = E10.5_C3_gene,
                                  "E11.5 Cluster4" = E11.5_C4_gene,
                                  "E14.5 Cluster1" = E14.5_C1_gene), 
             filename = NULL, fill = c('#9FA7D9', '#3F51B4',"#273492"), 
             alpha = c(0.5), cat.col = rep('black', 3), 
             col = 'black', cex = 1.5, fontfamily = 'sans', 
             cat.cex = 1.5, cat.fontfamily = 'sans')
ggsave(pattern3_venn, file="pattern3_venn.svg", device = "svg")
