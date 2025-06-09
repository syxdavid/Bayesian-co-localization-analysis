library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

file_path <- file.choose()
data <- read.table(file_path, header = TRUE, sep = ",")

clean_data <- subset(data, grepl("^[A-Za-z0-9]+$", gene))
gene_list <- clean_data$gene

genes_mapped <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ids <- genes_mapped$ENTREZID

ego_bp <- enrichGO(gene = gene_ids, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_cc <- enrichGO(gene = gene_ids, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_mf <- enrichGO(gene = gene_ids, OrgDb = "org.Hs.eg.db", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

ego_bp_df <- head(as.data.frame(ego_bp), 10)
ego_cc_df <- head(as.data.frame(ego_cc), 10)
ego_mf_df <- head(as.data.frame(ego_mf), 10)
ego_bp_df$category <- "BP"
ego_cc_df$category <- "CC"
ego_mf_df$category <- "MF"
ego_df <- rbind(ego_bp_df, ego_cc_df, ego_mf_df)

ego_df <- na.omit(ego_df)
ego_df$logP <- -log10(ego_df$p.adjust)

ego_df <- na.omit(ego_df)

ego_df$logP <- -log10(ego_df$p.adjust)

colors <- c("#3389bc", "#5496b0", "#6daba6", "#7cc2a0", "#91cd97", "#a2d296", 
            "#b5d99b", "#c7e197", "#dae797", "#e4e997", "#eeea93", "#f4e28c", 
            "#fedb89", "#fdc77d", "#fcb074", "#f89e6a", "#f78c62", "#f47f5c", 
            "#ed6e5b", "#e25b58", "#d44855")

my_colors <- colorRampPalette(colors)(100)

p <- ggplot(ego_df, aes(x = GeneRatio, y = reorder(Description, -logP), size = Count, color = logP)) +
  geom_point() +
  scale_color_gradientn(colors = my_colors) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "GO Analysis",
       x = "Gene Ratio",
       y = "GO Term",
       color = "-log10(q-value)",
       size = "Gene Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(2, "lines"))

print(p)
