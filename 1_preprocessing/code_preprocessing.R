library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)

## 
### 
dir <- dir("./")
# 
samples_name = dir
# 
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}
### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
##############################################################################################################
### 
scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<10000&percent.mito<10)
scRNA = SCTransform(scRNA, vars.to.regress="percent.mito", verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
colnames(scRNA@meta.data)[1] = "Sample"
scRNA = RunHarmony(scRNA, group.by.vars="Sample", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
##############################################################################################################
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
colors = c4a("bold", 12)
colors = c("skyblue", "springgreen", "steelblue", "#F2B701", "pink", "#11A579", "#E68310", "#008695", "yellow", "#F97B72", "blueviolet", "#A5AA99")
######################################  
p = DimPlot(scRNA, reduction="umap", group.by="Sample", pt.size=1)+
theme(legend.position="right", plot.title=element_blank())+
scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
ggsave("UMAP_Sample.pdf", p, width=7, height=6)

p = VlnPlot(scRNA, features=c("nFeature_RNA"), cols=colors, pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))
ggsave("./violin_nFeature_RNA.pdf", p, width=6, height=6)
p = VlnPlot(scRNA, features=c("percent.mito"), cols=colors, pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))
ggsave("./violin_percent_mito.pdf", p, width=6, height=6)

#############################  
mydata <- FindClusters(scRNA, resolution=0.1)
mydata = subset(mydata, seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "./All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
#################################################################################
# 
FeaturePlot(mydata, features=c("Ly6d"), cols=c("lightgray", "red"))+NoLegend()
VlnPlot(mydata, features=c("AGER"), pt.size=0, cols=colors)+NoLegend()+theme(axis.title.x=element_blank())

# 0: NK/T cells: CCL5, GNLY, NKG7, CD3D, GZMA, TRAC
# 1: Mast cells: TPSB2, CPA3, CTSG
# 2: B cells: MS4A1, CD79A
# 3: Fibroblasts: DCN, COL1A1, LUM
# 4: Macrophages: LYZ, C1QA, C1QB, AIF1, CTSB
# 5: Neutrophils: NAMPT, CXCL8, G0S2
# 6: Epithelial cells: KRT14, KRT13, KRT5, KRT6A, KRT17, 
# 7: Plasma B cells: IGKC, MZB1, JCHAIN, IGHG1
# 8: Proliferative B cells: MKI67, TOP2A
# 9: Endothelial cells: VWF, PECAM1, CDH5
#10: Smooth muscle cells: TAGLN, ACTA2, MYL9, MYLK

cell_label = c(
"NK/T cells", "Mast cells", "B cells", "Fibroblasts", "Macrophages", "Neutrophils", "Epithelial cells",
"Plasma B cells", "Proliferative B cells", "Endothelial cells", "Smooth muscle cells"
)
#################################################################################
##
names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
UMAPPlot(mydata, pt.size=1, label=T, label.size=4)+NoLegend()+scale_color_manual(values=colors)+theme(text=element_text(family="Times"))
Type_arr = c()
for (i in 1:nrow(mydata@meta.data)) {
	if (mydata@meta.data$Sample[i]=="GSM5900215") {
		Type_arr[i] = "tumor tissues"
	} else if (mydata@meta.data$Sample[i]=="GSM5900216") {
		Type_arr[i] = "peri-tumor tissues"
	} else if (mydata@meta.data$Sample[i]=="GSM5900217") {
		Type_arr[i] = "tumor tissues"
	} else if (mydata@meta.data$Sample[i]=="GSM5900218") {
		Type_arr[i] = "peri-tumor tissues"
	} else if (mydata@meta.data$Sample[i]=="GSM5900219") {
		Type_arr[i] = "tumor tissues"
	} else if (mydata@meta.data$Sample[i]=="GSM5900220") {
		Type_arr[i] = "peri-tumor tissues"
	}
}
mydata@meta.data[["Type"]] = Type_arr
saveRDS(mydata, "./mydata_cluster.rds")

genes = c("CD3D", "NKG7", "GNLY", "TPSB2", "CPA3", "MS4A1", "CD79A", "DCN", "LUM", "LYZ", "AIF1", "C1QA", "CXCL8", "G0S2", "KRT14", "KRT5", "MZB1", "JCHAIN", "MKI67", "TOP2A", "VWF", "PECAM1", "CDH5", "TAGLN", "ACTA2", "MYL9")
p = DotPlot(mydata, features=genes)+coord_flip()+scale_color_distiller(palette="RdYlBu")+
theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave("marker_dotplot.pdf", p, width=7, height=6)

p = DoHeatmap(mydata, features=genes, group.colors=colors, label=FALSE)+scale_fill_gradientn(colors=c("white", "snow", "firebrick3"))+theme(text=element_text(family="Times"))
ggsave("marker_heatmap.pdf", p, width=10, height=5)

####################################################################### 
bar = mydata@meta.data %>% group_by(Type, cell_type) %>% count()
bar$cell_type = factor(bar$cell_type, levels=cell_label)
p = ggplot(data=bar, aes(x=Type, y=n, fill=cell_type))+ 
geom_bar(stat="identity", position=position_fill())+
scale_fill_manual(values=colors)+theme_bw()+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=15, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("barplot_cellType.pdf", p, width=6, height=6)

####################################################################### 
bar = mydata@meta.data %>% group_by(Sample, cell_type) %>% count()
bar$cell_type = factor(bar$cell_type, levels=cell_label)
p = ggplot(data=bar, aes(x=Sample, y=n, fill=cell_type))+ 
geom_bar(stat="identity", position=position_fill())+
scale_fill_manual(values=colors)+theme_classic()+
theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=15, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("barplot_Sample.pdf", p, width=6, height=6)



