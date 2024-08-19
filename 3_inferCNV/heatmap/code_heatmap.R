library(pheatmap)

df = read.table("./heatmap_amplification.txt", header=T, row.names=1, sep="\t")
p = pheatmap(df, angle_col="45", color=colorRampPalette(c("white", "red"))(20), cluster_rows=F, cluster_cols=F, border_color="lightgray")
ggsave("heatmap_amplification.pdf", p, width=3, height=8)


df = read.table("./heatmap_deletion.txt", header=T, row.names=1, sep="\t")
p = pheatmap(df, angle_col="45", color=colorRampPalette(c("white", "blue"))(20), cluster_rows=F, cluster_cols=F, border_color="lightgray")
ggsave("heatmap_deletion.pdf", p, width=3, height=8)


