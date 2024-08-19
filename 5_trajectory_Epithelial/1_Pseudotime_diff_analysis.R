# 
# 
# 
# (1)Basic Differential Analysis
# (2)Finding Genes that Distinguish Cell Type or State
# (3)Finding Genes that Changes as a Function of Pseudotime
# 
Time_diff = differentialGeneTest(cds[ordergene, ], cores=1, fullModelFormulaStr="~sm.ns(Pseudotime)")
## 
Time_diff = subset(Time_diff, qval<0.01)
Time_diff = arrange(Time_diff, qval)
Time_diff = Time_diff[, c(5,2,3,4,1,6,7)]  #
write.table(Time_diff, "./Time_diff.txt", col.names=T, row.names=T, quote=F, sep="\t")

Time_genes = Time_diff %>% pull(gene_short_name) %>% as.character()  #
plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters=5, show_rownames=F, return_heatmap=T)

### 
p = plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters=5, show_rownames=F, return_heatmap=T)
ggsave("./heatmap_pseudotime.pdf", p, width=3, height=6)
clusters = cutree(p$tree_row, k=5)
clustering = data.frame(clusters)
clustering[, 1] = as.character(clustering[, 1])
colnames(clustering) = "Gene_Clusters"   #
write.table(clustering, "./pseudotime_clusters.txt", col.names=T, row.names=T, sep="\t", quote=F)


