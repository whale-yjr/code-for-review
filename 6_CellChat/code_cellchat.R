library(CellChat)
#################### 
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human
#####  
> unique(CellChatDB$interaction$annotation)
[1] "Secreted Signaling"	"ECM-Receptor"	"Cell-Cell Contact"  #
CellChatDB.use <- subsetDB(CellChatDB, search="Cell-Cell Contact")
cellchat@DB <- CellChatDB.use
#####  
cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)
## SeuratFindMarkers，
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)  #
## cellchat@LR$LRsig
cellchat = projectData(cellchat, PPI.human)
## 
## 
#####################################################################################
### 
cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)
# 
cellchat = filterCommunication(cellchat, min.cells=10)
#####  
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
#####################################################################################
#######  
## 

netVisual_bubble(cellchat, sources.use=c(5, 1, 2, 3), targets.use=c(4))+
theme(axis.text.x=element_text(angle=15, hjust=0.5, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))

# 

