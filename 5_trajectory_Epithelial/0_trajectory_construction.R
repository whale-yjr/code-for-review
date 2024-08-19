# 
# 
# 
# 
# monocle
library(monocle)
mydata = readRDS("./Epithelial_trajectory.rds") #
# 
expr_matrix = as(as.matrix(mydata$SCT@counts), "sparseMatrix")
# 
p_data = mydata@meta.data

# 
f_data = data.frame(gene_short_name=row.names(mydata), row.names=row.names(mydata))
# e
# 
# 
pd = new("AnnotatedDataFrame", data=p_data)
fd = new("AnnotatedDataFrame", data=f_data)

cds = newCellDataSet(expr_matrix, phenoData=pd, featureData=fd, lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())
# 
# negbinomial.size() negbinomial(): 
# negbinomial

######################################################################################
# 
# 
# 
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

# 
# 
# 
cds = detectGenes(cds, min_expr=0.1)
expressed_genes = row.names(subset(fData(cds), num_cells_expressed>=10)) #

########### 

# 
# step 1: choosing genes that define progress
# step 2: reducing the dimensionality of the data
# step 3: ordering the cells in pseudotime

# Step 1: 
# 
# 
# 
#
# disp_table = dispersionTable(cds)
# disp.genes = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# 
# cds = setOrderingFilter(cds, disp.genes)

##
# DEG = FindMarkers(mydata, group.by="Type", ident.1="tumor tissues", ident.2="peri-tumor tissues", logfc.threshold=0, min.pct=0)
# express_genes = subset(DEG, abs(avg_logFC)>0.5&p_val_adj<0.05)$gene
# cds = setOrderingFilter(cds, express_genes)
########################################################################
## 
## 
##
## expressed_genes = VariableFeatures(mydata)

diff = differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr="~cell_type", cores=1)
## ~
## qval<1e-20,decreasing=F
deg = subset(diff, qval<1e-20)
deg = deg[order(deg$qval, decreasing=F),]
write.table(deg, "deg.txt", col.names=T, row.names=T, quote=F, sep="\t")
# 
ordergene = rownames(deg)
cds = setOrderingFilter(cds, ordergene)
## 
## cds@featureData@data[["use_for_ordering"]]
## table(cds@featureData@data[["use_for_ordering"]])
## 
## !!!!!!!!! 
## Step 2:
## 
cds = reduceDimension(cds, max_components=2, method="DDRTree")

## Step 3: 

## 
cds = orderCells(cds, root_state=2)
# root_state
########################################################################
## 
# pseudotime
p = plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=T)+scale_color_distiller(palette="RdYlBu")+theme(text=element_text(family="Times"))
ggsave("./trejactory_pesudotime.pdf", p, width=6, height=6)
# 
# 

# 
colors = c("aquamarine", "cornflowerblue")
p = plot_cell_trajectory(cds, color_by="Type", size=1, show_backbone=T)+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("./trajectory_Type.pdf", p, width=6, height=6)

colors = c("aquamarine", "#F0E442", "#009E73", "#E69F00", "#56B4E9")
p = plot_cell_trajectory(cds, color_by="cell_type")+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("./trajectory_cellType.pdf", p, width=6, height=6)

#######################################################################################
# 
library(ggpubr)
df <- pData(cds) 
## pData(cds)cds@phenoData@data
colors = c("aquamarine", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
p = ggplot(df, aes(x=Pseudotime, y=cell_type, color=cell_type, fill=cell_type))+
scale_color_manual(values=colors)+
scale_fill_manual(values=colors)+
geom_density_ridges(alpha=0.8)+
theme_classic()+
theme(text=element_text(family="Times"), axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, face="bold"), axis.title.y=element_blank(), legend.text=element_text(size=12, face="bold"), legend.title=element_blank(), legend.position="top")
ggsave("./density_Type.pdf", p, height=3, width=6)


colors = c("aquamarine", "cornflowerblue")
set = c("CDCA8", "CDK6", "DCTN1", "MAP4", "MAEA")
plot_genes_in_pseudotime(cds[set, ], color_by="Type", cell_size=1)+scale_color_manual(values=colors)+theme(text=element_text(family="Times", face="bold", size=13))


set = c("COL18A1", "MMP2", "CAV1", "RBPJ", "CALD1")
plot_genes_in_pseudotime(cds[set, ], color_by="Type", cell_size=1)+scale_color_manual(values=colors)+theme(text=element_text(family="Times", face="bold", size=13))



