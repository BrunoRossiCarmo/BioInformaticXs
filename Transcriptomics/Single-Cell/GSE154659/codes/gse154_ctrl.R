#Bruno Rossi Carmo
#03/11/2021
#gse154659: 

#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\scRNA\\gse154659" #Meu
setwd(dir.main)
id <- "gse154659_Ctrl"

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

#Reading Data:
ATK <- readRDS("rds/Atf3_WT_KO_Raw_counts.RDS")
C57 <- readRDS("rds/C57_Raw_counts.RDS")
C57@Dimnames[[2]]
c57_true_Paclitaxel <-C57[,grep("Paclitaxel",C57@Dimnames[[2]])]
c57_true_Control <- C57[,grep("Naive",C57@Dimnames[[2]])]

#Creating Seurat Object 2:
drg_C57_Control <- CreateSeuratObject(counts = c57_true_Control , project = "C57_pa",min.features = 500)
drg_C57_Control@meta.data

#Analysing Just The Treat Data
"Mitochondrial genes:"
drg_C57_Control <- PercentageFeatureSet(drg_C57_Control, pattern = "^mt", col.name = "percent.mt") 
x <- drg_C57_Control@meta.data

#Data Exploration:
out <- paste("images/",id,"_ViolinPlot_PreFilt_CTRL.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
VlnPlot(drg_C57_Control , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#More data exploration:
out <- paste("images/",id,"_QC_ViolinPlot_crtl.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
plot1 <- FeatureScatter(drg_C57_Control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(drg_C57_Control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#ScTransform: Normalization
drg_C57_Control <- SCTransform(drg_C57_Control, vars.to.regress = "percent.mt")
saveRDS(drg_C57_Control, file = "rds/filter_norm_ctrl.rds")
drg_C57_Control <- readRDS("rds/filter_norm_ctrl.rds")

#PCA: Principal Component Analysis.
drg_C57_Control  <- RunPCA(drg_C57_Control , verbose = FALSE)

#Create UMAP
#PCS: 35 PCS.
npc <- 1:30
drg_C57_Control <- FindNeighbors(drg_C57_Control , dims = npc, verbose = FALSE)
drg_C57_Control <- FindClusters(drg_C57_Control , verbose = FALSE, resolution = 0.8)
levels(drg_C57_Control$seurat_clusters)

#UMAP geralmente fica melhor.
drg_C57_Control <- RunUMAP(drg_C57_Control, dims = npc, verbose = FALSE)

#Cell_Types:
names <- strsplit(row.names(drg_C57_Control@meta.data), "_")
types <- vector()
for(i in 1:length(names)){
        types[i] <- names[[i]][6]
}
drg_C57_Control@meta.data$celltype <- types
table(types)

#Image_UMAP:
"Plot Ctrl"
out <- paste("images\\",id,"_CTRLCluster_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(drg_C57_Control, reduction = "umap",pt.size=0.7,label = T)  
dev.off()

#Identification:
"0 ------------------- Hip = NF2"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==0])

"1 ------------------- Hip = Satglia"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==1])

"2 ------------------- Hip = PEP1"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==2])

"3 ------------------- Hip = NP"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==3])

"4 ------------------- Hip = NF1"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==4])

"5 ------------------- Hip = NP"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==5])

"6 ------------------- Hip = Schwann"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==6])

"7 ------------------- Hip = CLTMR1"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==7])

"8 ------------------- Hip = NP"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==8])

"9 ------------------- Hip = Schwann"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==9])

"10 ------------------- Hip = PEP2"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==10])

"11 ------------------- Hip = SST"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==11])

"12 ------------------- Hip = Fibroblast"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==12])

"13 ------------------- Hip = NP"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==13])

"14 ------------------- Hip = Fibroblast"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==14])

"15 ------------------- Hip = Endothelial"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==15])

"16 ------------------- Hip = p"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==16])

"17 ------------------- Hip = NF3"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==17])

"18 ------------------- Hip = Fibroblast"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==18])

"19 ------------------- Hip = Macrophage"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters==19])

"20 ------------------- Hip = Pericyte"
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters== 20])

"21 ------------------- Hip = Neutrophil "
table(drg_C57_Control@meta.data$celltype[drg_C57_Control@meta.data$seurat_clusters== 21])

#Idents:
names <- c("NF2","Satglia","PEP1","NP","NF1",
           "NP","Schwann","CLTMR1","NP","Schwann",
           "PEP2","SST","Fibroblast","NP","Fibroblast","Endothelial",
           "p","NF3","Fibroblast","Macrophage","Pericyte","Neutrophil")
names(names) <- levels(drg_C57_Control)
drg_C57_Control <- RenameIdents(drg_C57_Control, names)

#Plot2:
"Plot CTRL"
out <- paste("images\\",id,"_CTRLCluster_UMAP_annot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(drg_C57_Control, reduction = "umap",pt.size=0.7,label = T,label.size = 4)
dev.off()

#Plots:
DotPlot(drg_C57_Control,feature= c("Sucnr1","Ntrk2"))

#Resave:
saveRDS(drg_C57_Control, file = "rds/filter_norm_ctrl_end.rds")

