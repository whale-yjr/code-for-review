library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)

Epithelial = subset(mydata, cell_type=="Epithelial cells")
Epithelial = subset(Epithelial, Type=="tumor tissues")
Epithelial = SCTransform(Epithelial, vars.to.regress="percent.mito", verbose=FALSE)
Epithelial = RunPCA(Epithelial, verbose=FALSE, assay="SCT")
Epithelial = RunHarmony(Epithelial, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(Epithelial)
Epithelial <- FindNeighbors(Epithelial, dims=1:20, reduction="harmony")
Epithelial <- RunUMAP(Epithelial, dims=1:20, reduction="harmony")

colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
######################################  
DimPlot(Epithelial, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)
#################################################################################
#############################  
mydata <- FindClusters(Epithelial, resolution=0.2)
UMAPPlot(mydata, pt.size=2, label=T, cols=colors)+NoLegend()
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "markers.txt", col.names=T, row.names=T, quote=F, sep="\t")
#################################################################################
# 
VlnPlot(mydata, features=c("KRT6C"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())

cell_label = c(
"KRT6C+ Epithelial cells", "UBD+ Epithelial cells", "GSTA1+ Epithelial cells",
"GAS2L3+ Epithelial cells", "T cells"
)
#################################################################################
## 
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
p = UMAPPlot(mydata, pt.size=2, label=T, label.size=5)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("UMAP_Epithelial_sub.pdf", p, width=6, height=6)

set = c("KRT6C", "UBD", "GSTA1", "GAS2L3")
p = VlnPlot(mydata, features=set, pt.size=0, cols=colors, ncol=2)+NoLegend()+theme(text=element_text(family="Times"), axis.title.x=element_blank())
ggsave("Violin_Epithelial_sub.pdf", p, width=8, height=7)
############################################################################


