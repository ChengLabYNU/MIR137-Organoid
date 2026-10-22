library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(ggsignif)
library(patchwork)
library(scales)

color_map <- c(WT = '#9f9f9f', Het = '#7ca4c5', KO = '#50b8e7')
fmt_p <- function(p) ifelse(is.na(p), "NS.", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "NS."))))
save_plot <- function(p, name, width = 2.8, height = 3.5) {
  ggsave(paste0(name, ".png"), p, width = width, height = height)
  ggsave(paste0(name, ".pdf"), p, width = width, height = height)
}

plot_mixed <- function(df, yvar, ylab) {
  df <- df %>% mutate(genotype = droplevels(factor(genotype)))
  gs <- levels(df$genotype)
  df_org <- df %>% group_by(genotype, organoid_id) %>%
    summarise(mean_ratio = mean(.data[[yvar]], na.rm = TRUE), .groups = "drop")
  pairs <- if (length(gs) == 2) list(gs) else list(c("WT", "Het"), c("Het", "KO"), c("WT", "KO"))
  p_raw <- sapply(pairs, function(p) t.test(mean_ratio ~ genotype, data = subset(df_org, genotype %in% p))$p.value)
  get_p <- function(g1, g2) { i <- which(sapply(pairs, setequal, c(g1, g2))); if (length(i)) p_raw[i] else NA }
  labs_sig <- sapply(pairs, function(p) fmt_p(get_p(p[1], p[2])))
  y_max <- max(df[[yvar]], na.rm = TRUE); y_step <- y_max * 0.10
  ys <- y_max + y_step * seq_along(pairs); y_top <- y_max + y_step * (length(pairs) + 1)
  p <- ggplot() +
    stat_summary(data = df_org, aes(x = genotype, y = mean_ratio), fun = median, fun.min = median, fun.max = median,
                 geom = "errorbar", color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
    geom_quasirandom(data = df, aes(x = genotype, y = .data[[yvar]], color = genotype), size = 1.2, alpha = 0.25, width = 0.35) +
    geom_quasirandom(data = df_org, aes(x = genotype, y = mean_ratio, color = genotype), size = 3, alpha = 0.9, width = 0.25)
  for (i in seq_along(pairs)) {
    x1 <- match(pairs[[i]][1], gs); x2 <- match(pairs[[i]][2], gs)
    p <- p +
      annotate("segment", x = x1, xend = x1, y = ys[i] - y_step * 0.2, yend = ys[i], linewidth = 0.4) +
      annotate("segment", x = x2, xend = x2, y = ys[i] - y_step * 0.2, yend = ys[i], linewidth = 0.4) +
      annotate("segment", x = x1, xend = x2, y = ys[i], yend = ys[i], linewidth = 0.4) +
      annotate("text", x = (x1 + x2) / 2, y = ys[i] + y_step * 0.05, label = labs_sig[i], size = 5)
  }
  p + scale_color_manual(values = color_map) +
    scale_y_continuous(limits = c(0, y_top), expand = expansion(mult = c(0.02, 0.02)), labels = scales::percent_format()) +
    labs(y = ylab) + theme_classic() +
    theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          legend.position = "none")
}

# 00.Diameter
Diameter_data <- read.csv('00.Diameter/Data.csv', header = TRUE)
Diameter_data$Diameter <- as.numeric(Diameter_data$Diameter)
Diameter_data$Genotype <- factor(Diameter_data$Genotype, levels = c('WT','Het','KO'))
C_D10 <- filter(Diameter_data, Time == 'Day10', Region == 'Cerebral')
D_D25 <- filter(Diameter_data, Time == 'Day25', Region == 'Dorsal')
D_D50 <- filter(Diameter_data, Time == 'Day50', Region == 'Dorsal')
V_D25 <- filter(Diameter_data, Time == 'Day25', Region == 'Ventral')
V_D50 <- filter(Diameter_data, Time == 'Day50', Region == 'Ventral')
C_D10_plot <- ggplot(C_D10, aes(x = Genotype, y = Diameter)) +
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
  geom_quasirandom(aes(color = Genotype), size = 2, alpha = 0.75) +
  geom_signif(comparisons = list(c("WT", "KO")), step_increase = 0.1, map_signif_level = T) +
  scale_color_manual(values = color_map) + labs(y = "Diameter (μm)") + scale_y_continuous(limits = c(0, NA)) +
  theme_classic() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(color = "black"),
                          axis.text.x = element_text(size = 12, color = "black", face = "bold"), legend.position = "none")
save_plot(C_D10_plot, 'Cerebral_D10_plot')
D_D25_plot <- ggplot(D_D25, aes(x = Genotype, y = Diameter)) +
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
  geom_quasirandom(aes(color = Genotype), size = 2, alpha = 0.75) +
  geom_signif(comparisons = list(c("WT", "Het"), c("Het", "KO"), c("WT", "KO")), step_increase = 0.1, map_signif_level = T) +
  scale_color_manual(values = color_map) + labs(y = "Diameter (μm)") + scale_y_continuous(limits = c(0, NA)) +
  theme_classic() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(color = "black"),
                          axis.text.x = element_text(size = 12, color = "black", face = "bold"), legend.position = "none")
save_plot(D_D25_plot, 'Dorsal_D25_plot')
D_D50_plot <- ggplot(D_D50, aes(x = Genotype, y = Diameter)) +
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
  geom_quasirandom(aes(color = Genotype), size = 2, alpha = 0.75) +
  geom_signif(comparisons = list(c("WT", "Het"), c("Het", "KO"), c("WT", "KO")), step_increase = 0.1, map_signif_level = T) +
  scale_color_manual(values = color_map) + labs(y = "Diameter (μm)") + scale_y_continuous(limits = c(0, NA)) +
  theme_classic() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(color = "black"),
                          axis.text.x = element_text(size = 12, color = "black", face = "bold"), legend.position = "none")
save_plot(D_D50_plot, 'Dorsal_D50_plot')
V_D25_plot <- ggplot(V_D25, aes(x = Genotype, y = Diameter)) +
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
  geom_quasirandom(aes(color = Genotype), size = 2, alpha = 0.75) +
  geom_signif(comparisons = list(c("WT", "Het"), c("Het", "KO"), c("WT", "KO")), step_increase = 0.1, map_signif_level = T) +
  scale_color_manual(values = color_map) + labs(y = "Diameter (μm)") + scale_y_continuous(limits = c(0, NA)) +
  theme_classic() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(color = "black"),
                          axis.text.x = element_text(size = 12, color = "black", face = "bold"), legend.position = "none")
save_plot(V_D25_plot, 'Ventral_D25_plot')
V_D50_plot <- ggplot(V_D50, aes(x = Genotype, y = Diameter)) +
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), color = "grey", width = 0.8, linetype = 1, linewidth = 3, alpha = 0.8) +
  geom_quasirandom(aes(color = Genotype), size = 2, alpha = 0.75) +
  geom_signif(comparisons = list(c("WT", "Het"), c("Het", "KO"), c("WT", "KO")), step_increase = 0.1, map_signif_level = T) +
  scale_color_manual(values = color_map) + labs(y = "Diameter (μm)") + scale_y_continuous(limits = c(0, NA)) +
  theme_classic() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.text = element_text(color = "black"),
                          axis.text.x = element_text(size = 12, color = "black", face = "bold"), legend.position = "none")
save_plot(V_D50_plot, 'Ventral_D50_plot')

# 01.1.Day25_Dorsal_Neural_Rosettes
NR_data <- read.csv('01.Day25_Dorsal_SOX2.Ki67/Neural Rosettes.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = strsplit(File, '-') %>% sapply('[', 9),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"))
p1 <- plot_mixed(NR_data, "Value", expression("Neural Rosettes per field"))
save_plot(p1, '01.1.Day25_Neural_Rosettes_plot')

# 01.Day25_Dorsal_SOX2.Ki67
SOX2.Ki67_data <- read.csv('01.Day25_Dorsal_SOX2.Ki67/Default_Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File))
df_SOX2.Ki67 <- SOX2.Ki67_data %>%
  mutate(organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_sox2 = AF555_and_DAPI_Area / DAPI_Area,
         ratio_Ki67 = AF488_and_DAPI_Area / DAPI_Area,
         ratio_sox2Ki67 = All_Three_Area / DAPI_Area)
p1 <- plot_mixed(df_SOX2.Ki67, "ratio_sox2", expression("SOX2"^"+" * "DAPI"^"+"* "/DAPI"^"+"))
p2 <- plot_mixed(df_SOX2.Ki67, "ratio_Ki67", expression("Ki67"^"+" * "DAPI"^"+" * "/DAPI"^"+"))
p3 <- plot_mixed(df_SOX2.Ki67, "ratio_sox2Ki67", expression("SOX2"^"+" *"Ki67"^"+" * "DAPI"^"+" * "/DAPI"^"+"))
save_plot(p1, '01.Dorsal_Day25_SOX2_plot')
save_plot(p2, '01.Dorsal_Day25_Ki67_plot')
save_plot(p3, '01.Dorsal_Day25_SOX2.Ki67_plot')

# 02.Day25_Dorsal_Ki67.Edu
Ki67.Edu_data <- read.csv('02.Day25_Dorsal_Edu.Ki67/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File))
df_Ki67.Edu <- Ki67.Edu_data %>%
  mutate(organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_Ki67_Edu = (AF488_and_DAPI_Area - All_Three_Area) / AF488_and_DAPI_Area)
p1 <- plot_mixed(df_Ki67.Edu, "ratio_Ki67_Edu", expression("Ki67"^"-" *"Edu"^"+" * "DAPI"^"+" * "/Total Edu"^"+"* "DAPI"^"+"))
save_plot(p1, '02.Dorsal_Day25_Ki67.Edu_plot')

# 03.Day25_Dorsal_PAX6
PAX6_data <- read.csv('03.Day25_Dorsal_PAX6/New_same_antibody/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_PAX6 = AF555_and_DAPI_Area / DAPI_Area)
WH_df_PAX6 <- filter(PAX6_data, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_PAX6, "ratio_PAX6", expression("PAX6"^"+" * " DAPI"^"+"* " / DAPI"^"+"))
save_plot(p1, '03.Dorsal_Day25_PAX6_plot')

# 04.Dorsal_Day50_PAX6
PAX6_data <- read.csv('04.Day50_Dorsal_PAX6/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_PAX6 = AF555_and_DAPI_Area / DAPI_Area)
WH_df_PAX6 <- filter(PAX6_data, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_PAX6, "ratio_PAX6", expression("PAX6"^"+" * "DAPI"^"+"* "/DAPI"^"+"))
save_plot(p1, '04.Dorsal_Day50_PAX6_plot')

# 05.Day25_Ventral_NKX2.1
NKX2.1_data1 <- read.csv('05.Day25_Ventral_NKX2.1/Batch1/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_NKX2.1 = AF555_and_DAPI_Area / DAPI_Area)
NKX2.1_data2 <- read.csv('05.Day25_Ventral_NKX2.1/Batch23/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_NKX2.1 = AF555_and_DAPI_Area / DAPI_Area)
NKX2.1_data <- rbind(NKX2.1_data1, NKX2.1_data2[, colnames(NKX2.1_data1)])
WH_df_NKX2.1 <- filter(NKX2.1_data, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_NKX2.1, "ratio_NKX2.1", expression("NKX2.1"^"+"* "DAPI"^"+" *"/DAPI"^"+"))
save_plot(p1, '05.Dorsal_Day25_NKX2.1_plot')

# 06.Day50_Ventral_NKX2.1
NKX2.1_data <- read.csv('06.Day50_Ventral_NKX2.1/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_NKX2.1 = AF555_and_DAPI_Area / DAPI_Area)
WH_df_NKX2.1 <- filter(NKX2.1_data, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_NKX2.1, "ratio_NKX2.1", expression("NKX2.1"^"+"* "DAPI"^"+" *"/DAPI"^"+"))
save_plot(p1, '06.Dorsal_Day50_NKX2.1_plot')

# 07.Day25_Dorsal_DCX.MAP2
DCX.MAP2_data <- read.csv('07.Day25_Dorsal_DCX.MAP2/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File))
df_DCX.MAP2 <- DCX.MAP2_data %>%
  mutate(organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_dcx = AF555_Area / DAPI_Area,
         ratio_map2 = AF488_Area / DAPI_Area,
         ratio_dcxmap2 = AF555_and_AF488_Area / DAPI_Area)
WH_df_DCX.MAP2 <- filter(df_DCX.MAP2, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_DCX.MAP2, "ratio_dcx", expression("DCX"^"+"* "/DAPI"^"+"))
p2 <- plot_mixed(WH_df_DCX.MAP2, "ratio_map2", expression("MAP2"^"+"*"/DAPI"^"+"))
p3 <- plot_mixed(WH_df_DCX.MAP2, "ratio_dcxmap2", expression("DCX"^"+" *"MAP2"^"+" * "/DAPI"^"+"))
save_plot(p1, '07.Dorsal_Day25_DCX_plot')
save_plot(p2, '07.Dorsal_Day25_MAP2_plot')
save_plot(p3, '07.Dorsal_Day25_DCX.MAP2_plot')

# 08.Day50_Dorsal_DCX.MAP2
DCX.MAP2_data <- read.csv('08.Day50_Dorsal_DCX.MAP2/10_Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         section_id = gsub('-..czi','',File))
df_DCX.MAP2 <- DCX.MAP2_data %>%
  mutate(organoid_id = paste(batch, genotype, organoid, sep = "_"),
         ratio_dcx = AF555_Area / DAPI_Area,
         ratio_map2 = AF488_Area / DAPI_Area,
         ratio_dcxmap2 = AF555_and_AF488_Area / DAPI_Area)
WH_df_DCX.MAP2 <- filter(df_DCX.MAP2, genotype %in% c('WT','Het'))
p1 <- plot_mixed(WH_df_DCX.MAP2, organoid_id, "ratio_dcx", expression("DCX"^"+"* "/DAPI"^"+"))
p2 <- plot_mixed(WH_df_DCX.MAP2, organoid_id, "ratio_map2", expression("MAP2"^"+"*"/DAPI"^"+"))
p3 <- plot_mixed(WH_df_DCX.MAP2, organoid_id, "ratio_dcxmap2", expression("DCX"^"+" *"MAP2"^"+" * "/DAPI"^"+"))
save_plot(p1, '08.Dorsal_Day50_DCX_plot')
save_plot(p2, '08.Dorsal_Day50_MAP2_plot')
save_plot(p3, '08.Dorsal_Day50_DCX.MAP2_plot')

# Day75 / Day100 GFAP
V_d75_gfap_czi_data <- read.csv('11.Day75_Ventral_GFAP/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('_c1-2.jpg','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         batch = factor(batch),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         AA = AF555_Area / DAPI_Area)
p3 <- plot_mixed(V_d75_gfap_czi_data, "AA", expression("GFAP"^"+" * "Area / DAPI"^"+" * " Area"))
save_plot(p3, 'Ventral_Day75_GFAP_czi_AA_plot')

V_d100_gfap_czi_data <- read.csv('12.Day100_Ventral_GFAP/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('_c1-2.jpg','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         batch = factor(batch),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         AA = AF555_Area / DAPI_Area)
p5 <- plot_mixed(V_d100_gfap_czi_data, "AA", expression("GFAP"^"+" * "Area / DAPI"^"+" * " Area"))
save_plot(p5, 'Ventral_Day100_GFAP_czi_AA_plot')

D_d75_gfap_czi_data <- read.csv('09.Day75_Dorsal_GFAP/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('.czi','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         batch = factor(batch),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         AA = AF555_Area / DAPI_Area)
p7 <- plot_mixed(D_d75_gfap_czi_data, "AA", expression("GFAP"^"+" * "Area / DAPI"^"+" * " Area"))
save_plot(p7, 'Dorsal_Day75_GFAP_czi_AA_plot')

D_d100_gfap_czi_data <- read.csv('10.Day100_Dorsal_GFAP/Area_Results.csv', header = TRUE) %>%
  mutate(cellline = strsplit(File, '-') %>% sapply('[', 1),
         batch = strsplit(File, '-') %>% sapply('[', 2),
         time = strsplit(File, '-') %>% sapply('[', 3),
         region = strsplit(File, '-') %>% sapply('[', 4),
         genotype = strsplit(File, '-') %>% sapply('[', 6),
         organoid = strsplit(File, '-') %>% sapply('[', 7),
         section = strsplit(File, '-') %>% sapply('[', 8),
         field = gsub('_c1-2.jpg','',strsplit(File, '-') %>% sapply('[', 9)),
         genotype = factor(genotype, levels = c('WT', 'Het', 'KO')),
         batch = factor(batch),
         organoid_id = paste(batch, genotype, organoid, sep = "_"),
         AA = AF555_Area / DAPI_Area)
p8 <- plot_mixed(D_d100_gfap_czi_data, "AA", expression("GFAP"^"+" * "Area / DAPI"^"+" * " Area"))
save_plot(pp85, 'Dorsal_Day100_GFAP_czi_AA_plot')