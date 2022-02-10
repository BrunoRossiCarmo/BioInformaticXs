#Bruno Rossi Carmo
#03/11/2021
#gse154659: 

#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\scRNA\\gse154659" #Meu
setwd(dir.main)
id <- "gse154659"

#Libraries:
#Biology:
library("BiocManager")
library("Seurat")
library(sctransform)

#Others:
library(devtools)
library(dplyr)
library(tidyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(readxl)

#For memory:
memory.limit(9999999999)

#Reading data:
ATK <- readRDS("rds/Atf3_WT_KO_Raw_counts.RDS")
C57 <- readRDS("rds/C57_Raw_counts.RDS")
c57_true_Paclitaxel <-C57[,grep("Paclitaxel",C57@Dimnames[[2]])]
c57_true_Control <- C57[,grep("Naive",C57@Dimnames[[2]])]

#Creating Seurat Object 1:
drg_C57_Paclitaxel <- CreateSeuratObject(counts = c57_true_Paclitaxel , project = "C57_pa",min.features = 500)

#Creating Seurat Object 2:
drg_C57_Control <- CreateSeuratObject(counts = c57_true_Control , project = "C57_pa",min.features = 500)

#Merging:
drg_C57_Control@meta.data$dataset <- "Ctrl"
drg_C57_Paclitaxel@meta.data$dataset <- "Paclitaxel"
total_drg <- merge(drg_C57_Paclitaxel,drg_C57_Control)
meta <- total_drg@meta.data

#Cell_Types:
names <- strsplit(row.names(total_drg@meta.data), "_")
types <- vector()
for(i in 1:length(names)){
        types[i] <- names[[i]][6]
}
total_drg@meta.data$celltype <- types
table(types)

#Filtering:
total_drg <- PercentageFeatureSet(total_drg, pattern = "^mt", col.name = "percent.mt") 

#Data Exploration:
out <- paste("images/",id,"_ViolinPlot_PreFilt_Tot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
VlnPlot(total_drg , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#ScTransform: Normalization
total_drg <- SCTransform(total_drg, vars.to.regress = "percent.mt")
saveRDS(total_drg, file = "rds/total_drg_filter_norm.rds")
total_drg <- readRDS("rds/total_drg_filter_norm.rds")

#PCA: Principal Component Analysis.
total_drg   <- RunPCA(total_drg  , verbose = FALSE)

#Create UMAP
#PCS: 35 PCS.
npc <- 1:30
total_drg  <- FindNeighbors(total_drg  , dims = npc, verbose = FALSE)
total_drg  <- FindClusters(total_drg  , verbose = FALSE, resolution = 0.8)
levels(total_drg $seurat_clusters)

#UMAP geralmente fica melhor.
total_drg  <- RunUMAP(total_drg , dims = npc, verbose = FALSE)

#Image_UMAP:
"Plot Total"
out <- paste("images\\",id,"_TotalCluster_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(total_drg, reduction = "umap",pt.size=0.7,label = T)  
dev.off()

#Identification:
"0 ------------------- Hip = PEP1"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==0])

"1 ------------------- Hip = Satglia "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==1])

"2 ------------------- Hip =  Ab-LTMR"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==2])

"3 ------------------- Hip =  NP "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==3])

"4 ------------------- Hip =  Ab-LTMR  "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==4])

"5 ------------------- Hip = Schwann "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==5])

"6 ------------------- Hip = CLTMR1 "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==6])

"7 ------------------- Hip = NP "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==7])

"8 ------------------- Hip =  Ab-LTMR  "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==8])

"9 ------------------- Hip = PEP2 "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==9])

"10 ------------------- Hip = SST"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==10])

"11 ------------------- Hip = NP/PEP1"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==11])

"12 ------------------- Hip = NP"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==12])

"13 ------------------- Hip = Schwann"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==13])

"14 ------------------- Hip = Fibroblast "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==14])

"15 ------------------- Hip = Fibroblast "
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==15])

"16 ------------------- Hip = pCLTMR"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==16])

"17 ------------------- Hip = Ad-LTMR"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==17])

"18 ------------------- Hip = Endothelial"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==18])

"19 ------------------- Hip = pCLTMR"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==19])

"20 ------------------- Hip = NP"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==20])

"21 ------------------- Hip = Fibroblast"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==21])

"22 ------------------- Hip = Macrophage"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==22])

"23 ------------------- Hip = Pericyte"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==23])

"24 ------------------- Hip = Neutrophil"
table(total_drg@meta.data$celltype[total_drg@meta.data$seurat_clusters==24])

#Idents:
names <- c("PEP1","Satglia","Ab-LTMR","NP","Ab-LTMR",
           "Schwann","CLTMR1","NP","Ab-LTMR","PEP2",
           "SST","NP/PEP1","NP","Schwann","Fibroblast",
           "Fibroblast","pCLTMR","Ad-LTMR","Endothelial","pCLTMR",
           "NP","Fibroblast","Macrophage","Pericyte","Neutrophil")
names(names) <- levels(total_drg)
total_drg <- RenameIdents(total_drg, names)

#charac:
neuron_cells <- c("NF2","PEP1","NP","NF1","CLTMR1","PEP2","SST","p","NF3")
total_drg@meta.data$Cell_Type <- "non_neuronal"
result <- grep(paste(neuron_cells, collapse="|"),total_drg@meta.data$celltype)
total_drg@meta.data$Cell_Type[result]<- "neuronal"
View(total_drg@meta.data)

#Image_UMAP:
"Plot Total"
out <- paste("images\\",id,"_TotalCluster_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(total_drg, reduction = "umap",pt.size=0.7,label = T)  
dev.off()

#Image_UMAP:
"Plot Total"
out <- paste("images\\",id,"_CellsCluster_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(total_drg, reduction = "umap",pt.size=0.7,label = T,group.by = "Cell_Type")  
dev.off()

#Image_UMAP:
"Plot Total"
out <- paste("images\\",id,"_CondCluster_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(total_drg, reduction = "umap",pt.size=0.7,label = T,group.by = "dataset")  
dev.off()

#Conts:
length(total_drg@meta.data$orig.ident[total_drg@meta.data$dataset=="Ctrl"])
length(total_drg@meta.data$orig.ident[total_drg@meta.data$dataset=="Paclitaxel"])

#DEG Analyses Prep:
obj.list <- SplitObject(total_drg, split.by = "dataset")
table(obj.list[[1]]@meta.data$dataset)
table(obj.list[[2]]@meta.data$dataset)

#Rename Idents:
levels(obj.list[[2]])
table(Idents(obj.list[[2]]))
names <- c("Ctrl_PEP1","Ctrl_Satglia","Ctrl_Ab-LTMR","Ctrl_NP",
           "Ctrl_Schwann","Ctrl_CLTMR1","Ctrl_PEP2",
           "Ctrl_SST","Ctrl_NP/PEP1","Ctrl_Fibroblast",
           "Ctrl_pCLTMR","Ctrl_Ad-LTMR","Ctrl_Endothelial","Ctrl_Macrophage",
           "Ctrl_Pericyte","Ctrl_Neutrophil")
names(names) <- levels(obj.list[[2]])
obj.list[[2]] <- RenameIdents(obj.list[[2]], names)
table(Idents(obj.list[[2]]))

#DEG:
total2_drg <- merge(obj.list[[1]],obj.list[[2]])
table(Idents(total2_drg))

#New Rds:
saveRDS(total2_drg, file = "rds/total_drg_end.rds")
total2_drg <- readRDS("rds/total_drg_end.rds")

#A-Beta:
deg_abeta <- FindMarkers(total2_drg, ident.1 = "Ab-LTMR", ident.2 = "Ctrl_Ab-LTMR", verbose = FALSE)
colnames(deg_abeta) <- c("P-Value","Avg_LogFC","Pct-Paclitaxel","Pct-Controle","P-Value_Ajd")
write.csv(deg_abeta,file = "table\\deg_abeta.csv",sep = "\t",row.names=T,col.names=T)

#A-Delta:
deg_adelta <- FindMarkers(total2_drg, ident.1 = "Ad-LTMR", ident.2 = "Ctrl_Ad-LTMR", verbose = FALSE)
colnames(deg_adelta) <- c("P-Value","Avg_LogFC","Pct-Paclitaxel","Pct-Controle","P-Value_Ajd")
write.csv(deg_adelta,file = "table\\deg_adelta.csv",sep = "\t",row.names=T,col.names=T)
