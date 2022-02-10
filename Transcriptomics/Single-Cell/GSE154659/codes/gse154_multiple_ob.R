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
library(ggpubr)
library(devtools)
library(dplyr)
library(tidyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(readxl)

#For memory:
memory.limit(9999999999)

#Read the data:
NORM <- readRDS("rds/filter_norm_ctrl.rds")
pacli <- readRDS("rds/filter_norm_end.rds")
View(NORM@meta.data)

#Create new id:
NORM@meta.data$type = "ctrl"
pacli@meta.data$type = "pacli"
dim(NORM@meta.data)

#Create total:
total <- merge(NORM,pacli)
View(total@meta.data)
saveRDS(total, file = "rds/total_end.rds")

#Create the full plot:
out <- paste("images\\",id,"_TotalCluster_UMAP_annot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(NORM, reduction = "umap",pt.size=0.7,label = T)
dev.off() #Deu erro.

#Create Diff:
ctrl <- DimPlot(NORM, reduction = "umap",pt.size=0.7,label = T)
paclitax <-  DimPlot(pacli, reduction = "umap",pt.size=0.7,label = T)
out <- paste("images\\",id,"_TotalCluster_UMAP_annot.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(ctrl, paclitax, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#charac:
neuron_cells <- c("NF2","PEP1","NP","NF1","CLTMR1","PEP2","SST","p","NF3")
NORM@meta.data$Cell_Type <- "non_neuronal"
pacli@meta.data$Cell_Type <- "non_neuronal"
result <- grep(paste(neuron_cells, collapse="|"),NORM@meta.data$celltype)
result2 <- grep(paste(neuron_cells, collapse="|"),pacli@meta.data$celltype)
NORM@meta.data$Cell_Type[result]<- "neuronal"
pacli@meta.data$Cell_Type[result2]<- "neuronal"

#Create Diff2:
ctrl <- DimPlot(NORM, reduction = "umap",pt.size=0.7,label = T,group.by = "Cell_Type")
paclitax <-  DimPlot(pacli, reduction = "umap",pt.size=0.7,label = T,group.by = "Cell_Type")
out <- paste("images\\",id,"_TotalNeuronCluster_UMAP_annot.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(ctrl, paclitax, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Create Diff3:
ctrl <- DimPlot(NORM, reduction = "umap",pt.size=0.7,label = T,group.by = "orig.ident")
paclitax <-  DimPlot(pacli, reduction = "umap",pt.size=0.7,label = T,group.by = "orig.ident")
out <- paste("images\\",id,"_TotalSexCluster_UMAP_annot.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(ctrl, paclitax, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#DotPlots:
genes <- c("Scn5a","Cacna1h","Kcnj4","Gabra1","Gabra2")
p1 <- DotPlot(NORM,feature= genes,idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes,idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Individual Plots:
p1 <- DotPlot(NORM,feature= genes[1],idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes[1],idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT1.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Individual Plots:
p1 <- DotPlot(NORM,feature= genes[2],idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes[2],idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT2.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Individual Plots:
p1 <- DotPlot(NORM,feature= genes[3],idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes[3],idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT3.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Individual Plots:
p1 <- DotPlot(NORM,feature= genes[4],idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes[4],idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT4.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Individual Plots:
p1 <- DotPlot(NORM,feature= genes[5],idents = neuron_cells,cols = c("blue4","red","gold"))
p2 <- DotPlot(pacli,features = genes[5],idents = neuron_cells,cols = c("blue4","red","gold"))
out <- paste("images\\",id,"_Total_DOT5.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(p1,p2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Violin Plot:
x1 <- VlnPlot(NORM,feature= genes[1],idents = neuron_cells)
x2 <- VlnPlot(pacli,feature= genes[1],idents = neuron_cells)
out <- paste("images\\",id,"_Total_Vln1.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(x1,x2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Violin Plot:
x1 <- VlnPlot(NORM,feature= genes[2],idents = neuron_cells)
x2 <- VlnPlot(pacli,feature= genes[2],idents = neuron_cells)
out <- paste("images\\",id,"_Total_Vln2.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(x1,x2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Violin Plot:
x1 <- VlnPlot(NORM,feature= genes[3],idents = neuron_cells)
x2 <- VlnPlot(pacli,feature= genes[3],idents = neuron_cells)
out <- paste("images\\",id,"_Total_Vln3.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(x1,x2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Violin Plot:
x1 <- VlnPlot(NORM,feature= genes[4],idents = neuron_cells)
x2 <- VlnPlot(pacli,feature= genes[4],idents = neuron_cells)
out <- paste("images\\",id,"_Total_Vln4.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(x1,x2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()

#Violin Plot:
x1 <- VlnPlot(NORM,feature= genes[5],idents = neuron_cells)
x2 <- VlnPlot(pacli,feature= genes[5],idents = neuron_cells)
out <- paste("images\\",id,"_Total_Vln5.png",sep="")
png(filename = out,units="in",width = 15,height = 7.5,res=300)
ggarrange(x1,x2, labels = c("CTRL","PACLITAXEL"),widths = c(20,20))
dev.off()