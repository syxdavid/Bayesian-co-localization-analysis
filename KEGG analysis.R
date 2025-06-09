library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

file_path <- file.choose()
data <- read.table(file_path, header = TRUE, sep = ",")

clean_data <- subset(data, grepl("^[A-Za-z0-9]+$", gene))

gene_list <- clean_data$gene

genes_mapped <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ids <- genes_mapped$ENTREZID

file_path <- file.choose()
data <- read.table(file_path, header = TRUE, sep = ",")
gene_ids <- data$gene_ids

kegg_result <- enrichKEGG(gene = gene_ids, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)

kegg_df <- as.data.frame(kegg_result)
kegg_df$logP <- -log10(kegg_df$p.adjust)

kegg_df_top20 <- head(kegg_df[order(kegg_df$logP, decreasing = TRUE),], 20)

colors <- c("#3389bc", "#5496b0", "#6daba6", "#7cc2a0", "#91cd97", "#a2d296", 
            "#b5d99b", "#c7e197", "#dae797", "#e4e997", "#eeea93", "#f4e28c", 
            "#fedb89", "#fdc77d", "#fcb074", "#f89e6a", "#f78c62", "#f47f5c", 
            "#ed6e5b", "#e25b58", "#d44855")

my_colors <- colorRampPalette(colors)(100)

p <- ggplot(kegg_df_top20, aes(x = GeneRatio, y = reorder(Description, -logP), size = Count, color = logP)) +
  geom_point() +
  scale_size(range = c(3, 7)) + 
  scale_color_gradientn(colors = my_colors) +
  facet_grid(scales = "free_y", space = "free_y") +
  labs(title = "KEGG Analysis",
       x = "Gene Ratio",
       color = "-log10(q-value)",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(2, "lines"))

print(p)
