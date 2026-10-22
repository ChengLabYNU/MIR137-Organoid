
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

RNAseq_DEG <- read.csv('dfo_D25_Het_WT_DEres.csv')
proteomics_DEP <- read.csv('proteome_data.csv')

proteomics_DEP <- proteomics_DEP %>%
  separate_rows(Gene.name, sep = ";")

proteomics_DEP$protein_pvalue <- 10^(-proteomics_DEP$p.value)

up_DEP <- filter(proteomics_DEP, Difference > 1, q.value < 0.05)
down_DEP <- filter(proteomics_DEP, Difference < -1, q.value < 0.05)

enrich_go_one <- function(genes, keyType, ont) {
  enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = keyType, ont = ont,
           pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")@result %>%
    filter(!is.na(p.adjust), p.adjust < 0.05) %>%
    slice_min(p.adjust, n = 5, with_ties = FALSE)
}

run_go <- function(genes, keyType, title, out, x_point = NULL) {
  df <- map_dfr(c("BP", "CC", "MF"), function(ont) {
    enrich_go_one(genes, keyType, ont) %>% mutate(ONTOLOGY = ont)
  })
  if (nrow(df) == 0) {
    message("No enriched terms: ", out)
    return(invisible(NULL))
  }
  df <- df %>%
    mutate(logP = -log10(p.adjust),
           Description = factor(Description, levels = arrange(df, desc(p.adjust))$Description))
  if (is.null(x_point)) x_point <- -0.03 * max(df$logP)
  p <- ggplot(df, aes(x = logP, y = Description, fill = ONTOLOGY)) +
    geom_col() +
    geom_point(aes(x = x_point, size = Count, color = ONTOLOGY)) +
    geom_text(aes(x = 0.1, label = Description), hjust = 0, vjust = 0.4) +
    geom_text(aes(x = x_point, label = Count)) +
    scale_fill_manual(values = c('#a5d6a6', '#7fdeea', '#9fa7d9')) +
    scale_color_manual(values = c('#a5d6a6', '#7fdeea', '#9fa7d9')) +
    facet_grid(ONTOLOGY ~ ., space = "free", scales = "free") +
    theme_classic() +
    labs(title = title,
         x = expression(-log[10]*paste('(', italic(P), '-adjust)'))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.x = element_text(color = 'black', size = 12),
          axis.title.y = element_blank(),
          axis.text.x = element_text(color = 'black', size = 8),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.length.y = unit(0, 'cm'))
  ggsave(out, p, width = 5, height = 6)
}

up_DEP_genes <- unique(unlist(strsplit(up_DEP$Protein.Group, ";")))
down_DEP_genes <- unique(unlist(strsplit(down_DEP$Protein.Group, ";")))

run_go(up_DEP_genes, 'UNIPROT',
       "Day 25 Het dorsal forebrain orgaoid \n up-regulated protein function enrich",
       'Padj_qvalue_Day25_Het_WT_up_DEP_GO_plot.pdf')

run_go(down_DEP_genes, 'UNIPROT',
       "Day 25 Het dorsal forebrain orgaoid \n down-regulated protein function enrich",
       'Padj_qvalue_Day25_Het_WT_down_DEP_GO_plot.pdf')

up_DEG <- filter(RNAseq_DEG, log2FoldChange > 1, padj < 0.05)
down_DEG <- filter(RNAseq_DEG, log2FoldChange < -1, padj < 0.05)

run_go(up_DEG$X, 'SYMBOL',
       "Day 25 Het dorsal forebrain orgaoid \n up-regulated gene function enrich",
       'Padj_Day25_Het_WT_up_DEG_GO_plot.pdf')

run_go(down_DEG$X, 'SYMBOL',
       "Day 25 Het dorsal forebrain orgaoid \n down-regulated gene function enrich",
       'Padj_Day25_Het_WT_down_DEG_GO_plot.pdf')




colnames(RNAseq_data) = c('gene','baseMean','RNA_lfc','lfcSE','stat','pvalue','padj')
proteomics_data = filter(proteomics_DEP, Gene.name %in% RNAseq_data$gene)
proteomics_data = proteomics_data[!duplicated(proteomics_data$Gene.name),]
colnames(proteomics_data) = c("Protein.Group","Protein.Ids","Protein.Names","ENSG","ENSP","ENST",
                              "Gene.name","p.value","protein_pvalue","q.value","protein_Difference",
                              "statistic","protein_pvalue")

RNAseq_data_f = filter(RNAseq_DEG, X %in% proteomics_DEP$Gene.name)
proteomics_data_f = filter(proteomics_DEP, Gene.name %in% RNAseq_DEG$X)


plot_data = merge(RNAseq_data_f,proteomics_data_f,by.x = 'X',by.y = 'Gene.name')
plot_data <- plot_data %>%
  mutate(
    RNA_group = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "RNA_up",
      padj < 0.05 & log2FoldChange < -1 ~ "RNA_down",
      TRUE ~ 'RNA_ns'),
    protein_group = case_when(
      q.value < 0.05 & Difference > 1 ~ "protein_up",
      q.value < 0.05 & Difference < -1 ~ "protein_down",
      TRUE ~ 'protein_ns'),
    group = case_when(
      RNA_group == "RNA_up" & protein_group == "protein_up" ~"both up",
      RNA_group == "RNA_down" & protein_group == "protein_down" ~"both down",
      TRUE ~ "NS"
    )
  )

target_human = read.csv('Q:/015-Reference/TargetScan8.0__miR-137.Human.predicted_targets.txt',
                        sep ='\t')$Target.gene

plot_data$co_dif = plot_data$log2FoldChange+plot_data$Difference
plot_data$name = ''
up_top_10 = top_n(x=filter(plot_data,group == 'both up'),10,wt = co_dif)$X
down_top_10 = top_n(x=filter(plot_data,group == 'both down'),-10,wt = co_dif)$X
plot_data$name[plot_data$X %in% c(up_top_10,down_top_10)] = plot_data$X[plot_data$X %in% c(up_top_10,down_top_10)]

library(ggrepel)
co_dif_plot = ggplot()+
  geom_point(data = filter(plot_data,group == "NS"),
             aes(x = Difference, y = log2FoldChange),
             color = "grey")+
  geom_point(data = filter(plot_data,group == "both up"),
             aes(x = Difference, y = log2FoldChange),size = 2,
             shape = 21,fill = "pink",color = "#ee5250")+
  geom_point(data = filter(plot_data,group == "both down"),
             aes(x = Difference, y = log2FoldChange),size = 2,
             shape = 21,fill = "lightblue",color = "#41a5f4")+
  geom_text_repel(data = plot_data,aes(x = Difference, y = log2FoldChange,label = name),
                  inherit.aes = F,show.legend = F)+
  geom_vline(xintercept = c(-1,1),lty = 2)+
  geom_hline(yintercept = c(-1,1),lty = 2)+
  theme_classic()

ggsave('co_dif_plot.pdf',co_dif_plot,width = 5,height = 5)
