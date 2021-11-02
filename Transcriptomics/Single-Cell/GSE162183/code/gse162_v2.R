#Bruno Rossi Carmo
#08/09/2021
#gse162183: 

#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\scRNA\\GSE162183" #Meu
setwd(dir.main)
id <- ""

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

#Organizando os Dados:
rawcounts <- read.table("/home/crid/Documents/bruno/analises/zeca-skin/scrna-seq/gse162183//data/GSE162183_Raw_gene_counts_matrix.tab", header = T,sep = "\t",row.names = 1)

#Sparse Matrix:
ps_sparse <-  Seurat::as.sparse(rawcounts)
saveRDS(ps_sparse,"rds/drg_sparse.rds")
ps_sparse <- readRDS("rds/drg_sparse.rds")

#Iniciando análise:
ps_x <- CreateSeuratObject(counts = ps_sparse , 
                           project = "ps_ZECA",
                            min.cells = 1000, 
                           min.features = 500)

#Mitocondrial:
ps_x  <- PercentageFeatureSet(ps_x, pattern = "^MT-", col.name = "percent.mt") #Something wrong here.
ps_x@meta.data
"Other way:"
mito.genes <- grep(pattern = "^MT-", rownames(ps_x), value = TRUE)
percent.mito <- (Matrix::colSums(ps_x[mito.genes, ])/Matrix::colSums(ps_x))*100
ps_x@meta.data$percent.mt <- percent.mito
"Comparar os valores de percent.mt e percent.mito"

#Data Exploration:
out <- paste("images/",id,"_ViolinPlot_PreFilt.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
VlnPlot(ps_x , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#More data exploration:
out <- paste("images/",id,"_QC_ViolinPlot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
plot1 <- FeatureScatter(ps_x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ps_x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#Filtering by cells and mt.mitochondrial:
ps_filtr <- subset(ps_x, subset = percent.mt < 5 & nCount_RNA < 75000 & nFeature_RNA < 7500)
VlnPlot(ps_filtr , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#ScTransform: Normalization
ps_filtr<- SCTransform(ps_filtr, vars.to.regress = "percent.mt")
saveRDS(ps_filtr, file = "rds/filter_norm.rds")
ps_filtr <- readRDS("code\\filter_norm_2.rds")

#PCA: Principal Component Analysis.
ps_filtr <- RunPCA(ps_filtr, verbose = FALSE)

#Create UMAP
#PCS: 35 PCS.
npc <- 1:35
ps_filtr <- FindNeighbors(ps_filtr, dims = npc, verbose = FALSE)
ps_filtr<- FindClusters(ps_filtr, verbose = FALSE, resolution = 1.7)
levels(ps_filtr$seurat_clusters)

#UMAP geralmente fica melhor.
ps_filtr <- RunUMAP(ps_filtr, dims = npc, verbose = FALSE)
ps_filtr <- RunTSNE(ps_filtr, dims = npc, verbose = FALSE)
ps_filtr@meta.data

#Image_UMAP:
"Plot Geral (Psor e Crtl juntas)"
out <- paste("images\\",id,"_Clusters_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.7,label = T)  
dev.off()

#Image_UMAP-2:
"Plot Ctrl"
ps_filtr@meta.data$cell <- "na"
ps_filtr@meta.data$cell[grep("Psor",ps_filtr@meta.data$orig.ident)] <- "Psor"
ps_filtr@meta.data$cell[grep("Ctrl",ps_filtr@meta.data$orig.ident)] <- "Ctrl"
out <- paste("images/",id,"_Clusters_UMAP_Ctrl.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.7,label = T,cells = ps_filtr@meta.data$cell == c("Ctrl"))  
dev.off()

#Image_UMAP-3:
"Plot Psor"
out <- paste("images/",id,"_Clusters_UMAP_Psor.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.7,label = T,cells = ps_filtr@meta.data$cell == c("Psor"))  
dev.off()

#Seguir com Análise:
"Identificação dos Clusters:"
Idents <- c("Mes","End","SchM","Mast Cells","IM_Mgl","T Cells",
            "T Cells","IM_DC1","IM_DC2","IM_Lg","T Cells",
            "Endo_Lym","Epd_Basal1","Epd_Basal2","Epd_Basal3",
           "Epd_Corneum","Epd_Foli","Epd_Granular","Epd_Granular/spinous",
           "Epd_SG","Epd_Spinous")
dim(ps_filtr)

#Possíveis Identificações:
"Mesenchymal -> GENES: (PDGFRB,ACTA2,DCN,ACTC1,LUM)"
"Hipótese: Cluster 0,28,27,11,4,13"
genes <- c("PDGFRB","ACTA2","DCN","ACTC1","LUM")
DotPlot(ps_filtr,features = genes)
FeaturePlot(ps_filtr,features=genes)

"End_Endo:"
genes <- c("PECAM1","CD36","LYVE1")
DotPlot(ps_filtr,features = genes)
FeaturePlot(ps_filtr,features=genes)
#Hip = 3,20,19,32,29,18,1

"Endo_Lym:"
genes <- c("PECAM1","CD36","LYVE1","FABP4","CD200")
DotPlot(ps_filtr,features = genes)
FeaturePlot(ps_filtr,features=genes)
#Hip = 25

"SchM:"
genes <- c("SOX10","NGFR","MLANA","PLP1")
DotPlot(ps_filtr,features = genes)
FeaturePlot(ps_filtr,features=genes)
#Hip = 31.

"Epd_Spinous"
genes <- c("KRT14","KRT5","KRT10","KRT1","GRHL1","KRT17","MGST1","KRT73")
DotPlot(ps_filtr,features = genes)
#hIP = 15

"Epd_Corneum"
"KRT14*,KRT5*,KRT10,KRT1,GRHL1,APOE,IEF,DRTDAP,KPRP-,FLG-,CHP2,GSTA3,MT5"
genes <- c("KRT10","KRT1","GRHL1","APOE","CHP2","GSTA3","MT5","KPRP","FLG")
DotPlot(ps_filtr,features = genes)
#hIP = 10/7

"Epd_Foli"
genes <- c("GRHL1","LGR5","KRT19","CD200","SOX9","KRT17","MGST1","NOTUM","CFTR","FOXC1","RUNX3")
DotPlot(ps_filtr,features = genes)
#hIP = 34

"Epd_Granular"
genes <- c("KRT14","KRT5","KRT10","KRT1","ITGA6","ITGB1GRHL1","APOE","KRT6","POSTN","KRT15")
DotPlot(ps_filtr,features = genes)
#hIP = 21

"Epd_Granular/Spinous"
genes <- c("KRT14","KRT5","KRT10","KRT1","GRHL1","KPRP","FLG")			
DotPlot(ps_filtr,features = genes)
#hIP = 26

"Epd_SG"
genes <- c("ABCC3","ACSL5","APOC1","SOX9","MYC")			
DotPlot(ps_filtr,features = genes)
#hIP = ?

"Epd_Basal"
genes <- c("KRT15","LGR6","KRT19","CD200","SOX9","MGST1",'CFTR',"FOXC1","RUNX3","KRT10","KRT1","LGR5")	
DotPlot(ps_filtr,features = genes)
#hIP = 2/21/34

"IM_CD4+T"
genes <- c("PTPRC","CD3G","CD3E","CD4")			
DotPlot(ps_filtr,features = genes)
#hIP = 5

"IM_CD8+T"
genes <- c("PTPRC","CD3G","CD3E","CD8A","CTLA4")			
DotPlot(ps_filtr,features = genes)
#hIP = 16

"IM_DC1"
genes <- c("ID1","ID3","CCN1","CSKMT","CD1C","THBD","CDH1","GPR157")			
DotPlot(ps_filtr,features = genes)
#hIP = 

"IM_DC2"
genes <- c("IL1B","CCR7","DUSP4","DUSP5","CXCL8","CD80","CD86","GPR158")		
DotPlot(ps_filtr,features = genes)
#hIP = 9

"IM_Lg"
genes <- c("CD207","PLEK2","CD1B","CD1A","CD1C","CXCL8")		
DotPlot(ps_filtr,features = genes)
#hIP = 9

"IM_Mast"
genes <- c("PTPRC","KIT","RAB38","MYO16")		
DotPlot(ps_filtr,features = genes)
#hIP = 30

"IM_Mgl1"
genes <- c("CD68S","CD163","LGMN","MS4A4A")		
DotPlot(ps_filtr,features = genes)
#hIP = 24

"IM_TREG"
genes <- c("PTPRC","CD3G","CD3E","CD4","CTLA4","FOXP3","IL2RA")		
DotPlot(ps_filtr,features = genes)
#hIP = 14

"Mes_Per"
genes <- c("PDGFRB","ACTA2","CCL26","EGFL6","FMO3","SEMA5B")		
DotPlot(ps_filtr,features = genes)
#hIP = 22/23/6

"Mes_Per2"
genes <- c("PDGFRB","ACTA2","CCL8","ADAMTS4","ACTN2","FOXC1")		
DotPlot(ps_filtr,features = genes)
#hIP = 12

"Mes_Fibro"
genes <- c("PDGFRB","LUM","DCN","C3","C7","APOD")		
DotPlot(ps_filtr,features = genes)
#hIP = 28/4/0

"Mes_SM"
genes <- c("ACTA2","ACTC1","PNCK","PCP4")		
DotPlot(ps_filtr,features = genes)
#hIP = 33

#Alterando Cluster Names:
Idents <- c("Mesenchymal Cells","Venular ECs","Basal Keratinocytes","Venular ECs","Mesenchymal Cells","T Cells","Mesenchymal Cells","Corneum Keratinocytes","Mesenchymal Cells","Myeloid cells","Corneum Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Mesenchymal Cells","T Cells","Spinous Keratinocytes","T Cells","Mesenchymal Cells","Venular ECs","Venular ECs","Venular ECs","Basal Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Myeloid cells","Lymphatic ECs","Granular Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Venular ECs","Mast Cells","Schwann Cells","Venular ECs","Mesenchymal Cells","Basal Keratinocytes")
names(Idents) <- levels(ps_filtr)
ps_filtr <- RenameIdents(ps_filtr, Idents)

#Image_UMAP_with_Annots:
"Plot Geral (Psor e Crtl juntas)"
library(tidyverse)
out <- paste("images\\",id,"_Clusters_UMAP_Annots.tiff",sep="")
tiff(filename = out,units="in",width = 12,height = 8,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.8) + 
ggtitle("Total")  +  theme(plot.title = element_text(hjust = 0.55)) +
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off()





#Com Anotações PT.1 (PKM2):
out <- paste("images\\",id,"_PKM_Dot.tiff",sep="")
tiff(filename = out,units="in",width = 10,height = 7,res=300)
DotPlot(ps_filtr,features = "PKM")
dev.off() 
out <- paste("images\\",id,"_PKM_Vln.tiff",sep="")
tiff(filename = out,units="in",width = 10,height = 7,res=300)
VlnPlot(ps_filtr,features="PKM")
dev.off() 

library(ggplot2)
#Plots Individuais:
#Ctrl:
out <- paste("images/",id,"_Clusters_UMAP_Ctrl.tiff",sep="")
tiff(filename = out,units="in",width = 12,height = 8,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.8,cells = ps_filtr@meta.data$cell == c("Ctrl")) + ggtitle("Ctrl")  +  theme(plot.title = element_text(hjust = 0.55)) +
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off()

dev.off()
#Psor:
out <- paste("images/",id,"_Clusters_UMAP_Psor.tiff",sep="")
tiff(filename = out,units="in",width = 12,height = 8,res=300)
DimPlot(ps_filtr, reduction = "umap",pt.size=0.8,cells = ps_filtr@meta.data$cell == c("Psor"))+
ggtitle("Psor")  +  theme(plot.title = element_text(hjust = 0.55)) +
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off()

dev.off()

#Plots:
#Violin's:
out <- paste("images\\",id,"_PKM_Vln.png",sep="")
png(filename = out,units="in",width = 12,height = 7,res=300)
VlnPlot(ps_filtr,features="PKM",split.by = "cell")
dev.off()
out <- paste("images\\",id,"_S100A9_Vln.png",sep="")
png(filename = out,units="in",width = 12,height = 7,res=300)
VlnPlot(ps_filtr,features="S100A9",split.by = "cell")
dev.off()
out <- paste("images\\",id,"_LCN2_Vln.png",sep="")
png(filename = out,units="in",width = 12,height = 7,res=300)
VlnPlot(ps_filtr,features="LCN2",split.by = "cell")
dev.off()
out <- paste("images\\",id,"_KRT14_Vln.png",sep="")
png(filename = out,units="in",width = 12,height = 7,res=300)
VlnPlot(ps_filtr,features="KRT14",split.by = "cell")
dev.off()
out <- paste("images\\",id,"_KRT17_Vln.png",sep="")
png(filename = out,units="in",width = 12,height = 7,res=300)
VlnPlot(ps_filtr,features="KRT17",split.by = "cell")
dev.off()
?VlnPlot
#Violin's Without Spinous (EDITAR):

#Feature Plot:
library("autoimage")
out <- paste("images\\",id,"_PKM_Feature.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="PKM",split.by = "cell", cols = c("grey","brown4"), pt.size = 1) +
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

dev.off() 
out <- paste("images\\",id,"_S100A9_Feature.png",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="S100A9",split.by = "cell", cols = c("grey","brown4") , pt.size = 1)+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off() 
out <- paste("images\\",id,"_LCN2_Feature.png",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="LCN2",split.by = "cell", cols = c("grey","brown4") , pt.size = 1)+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off() 
out <- paste("images\\",id,"_KRT14_Feature.png",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="KRT14",split.by = "cell", cols = c("grey","brown4") , pt.size = 1)+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off() 
out <- paste("images\\",id,"_KRT17_Feature.png",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="KRT17",split.by = "cell", cols = c("grey","brown4") ,  pt.size = 1)+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)
dev.off() 



#NEW PLOTS:
where <- FeaturePlot(obj.list[2]$Psor,features="PKM", cols = c("grey","brown4"), pt.size = 1)  + ggtitle("Psor PKM")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

wheter <- FeaturePlot(obj.list[1]$Ctrl,features="PKM", cols = c("grey","brown4"), pt.size = 1)  +
        ggtitle("Ctrl PKM")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

out <- paste("images\\",id,"_PKM_Feature.tiff",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
ggarrange(where,wheter)
dev.off()

where <- FeaturePlot(obj.list[2]$Psor,features="KRT14", cols = c("grey","brown4"), pt.size = 1)  + ggtitle("Psor KRT14")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

wheter <- FeaturePlot(obj.list[1]$Ctrl,features="KRT14", cols = c("grey","brown4"), pt.size = 1)  +
        ggtitle("Ctrl KRT14")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

out <- paste("images\\",id,"_KRT14_Feature.tiff",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
ggarrange(where,wheter)
dev.off()

where <- FeaturePlot(obj.list[2]$Psor,features="LCN2", cols = c("grey","brown4"), pt.size = 1)  + ggtitle("Psor LCN2")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

wheter <- FeaturePlot(obj.list[1]$Ctrl,features="LCN2", cols = c("grey","brown4"), pt.size = 1)  +
        ggtitle("Ctrl LCN2")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

out <- paste("images\\",id,"_LCN2_Feature.tiff",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
ggarrange(where,wheter)
dev.off()

where <- FeaturePlot(obj.list[2]$Psor,features="S100A9", cols = c("grey","brown4"), pt.size = 1)  + ggtitle("Psor S100A9")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

wheter <- FeaturePlot(obj.list[1]$Ctrl,features="S100A9", cols = c("grey","brown4"), pt.size = 1)  +
        ggtitle("Ctrl S100A9")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

out <- paste("images\\",id,"_S100A9_Feature.tiff",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
ggarrange(where,wheter)
dev.off()

where <- FeaturePlot(obj.list[2]$Psor,features="KRT17", cols = c("grey","brown4"), pt.size = 1)  + ggtitle("Psor KRT17")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -3, 11, label='Spinous Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

wheter <- FeaturePlot(obj.list[1]$Ctrl,features="KRT17", cols = c("grey","brown4"), pt.size = 1)  +
        ggtitle("Ctrl KRT17")  +  theme(plot.title = element_text(hjust = 0.5))+
        annotate("label", 6, 7 , label='Mesenchymal Cells',col='black',size = 4)+ 
        annotate("label", 9, -3 , label='Mesenchymal Cells',col='black',size = 4)+
        annotate("label", 1.4, 2.7 , label='Lymphatic ECs',col='black',size = 4)+
        annotate("label", -0.3, -8 , label='Venular ECs',col='black',size = 4)+
        annotate("label", -2.3, 0.6 , label='Schwann Cells',col='black',size = 4)+
        annotate("label", -6.3, -5, label='Mast Cells',col='black',size = 4)+
        annotate("label", -8.5, -1.5, label='T Cells',col='black',size = 4)+
        annotate("label", -9, -9.6, label='Myeloid cells',col='black',size = 4)+
        annotate("label", -3, 7, label='Basal Keratinocytes',col='black',size = 4)+
        annotate("label", -9.5, 10, label='Granular Keratinocytes',col='black',size = 4)+
        annotate("label", -9, 6, label='Corneum Keratinocytes',col='black',size = 4)

out <- paste("images\\",id,"_KRT17_Feature.tiff",sep="")
png(filename = out,units="in",width = 17,height = 7,res=300)
ggarrange(where,wheter)
dev.off()


dev.off()
?FeaturePlot

#FeaturePlot2:
out <- paste("images\\",id,"_PKM_FeatureT.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="PKM", cols= c("grey","brown4") , label = T, pt.size = 1)
dev.off() 
out <- paste("images\\",id,"_S100A9_FeatureT.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="S100A9", cols = c("grey","brown4") , label = T, pt.size = 1)
dev.off() 
out <- paste("images\\",id,"_LCN2_FeatureT.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="LCN2", cols = c("grey","brown4") , label = T, pt.size = 1)
dev.off() 
out <- paste("images\\",id,"_KRT14_FeatureT.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="KRT14", cols = c("grey","brown4") , label = T, pt.size = 1)
dev.off() 
out <- paste("images\\",id,"_KRT17_FeatureT.tiff",sep="")
tiff(filename = out,units="in",width = 17,height = 7,res=300)
FeaturePlot(ps_filtr,features="KRT17", cols = c("grey","brown4") , label = T, pt.size = 1)
dev.off()


Idents <- c("Mesenchymal Cells","Venular ECs","Basal Keratinocytes","Venular ECs","Mesenchymal Cells","T Cells","Mesenchymal Cells","Corneum Keratinocytes","Mesenchymal Cells","Myeloid cells","Corneum Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Mesenchymal Cells","T Cells","Spinous Keratinocytes","T Cells","Mesenchymal Cells","Venular ECs","Venular ECs","Venular ECs","Basal Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Myeloid cells","Lymphatic ECs","Granular Keratinocytes","Mesenchymal Cells","Mesenchymal Cells","Venular ECs","Mast Cells","Schwann Cells","Venular ECs","Mesenchymal Cells","Basal Keratinocytes")




#DotPlot:
order1<- c("Basal Keratinocytes","Spinous Keratinocytes","Corneum Keratinocytes","Granular Keratinocytes","Mast Cells","Myeloid cells","T Cells","Mesenchymal Cells","Lymphatic ECs","Venular ECs","Schwann Cells")
order2 <- c("Schwann Cells","Venular ECs","Lynphatic ECs","Mesenchymal Cells","Myeloid cells","T Cells","Mast Cells","Granular Keratinocytes","Corneum Keratinocytes","Spinous Keratinocytes","Basal Keratinocytes")
ps_filtr@active.ident<- factor(x = ps_filtr@active.ident, levels = order1)
levels(ps_filtr@active.ident)
obj.list <- SplitObject(ps_filtr, split.by = "cell")
obj.list[1]
obj.list[2]
#DOTPLOT CTRL:
                           
genes <- c("PKM")
out <- paste("images\\",id,"_Dot_Ctrl.tiff",sep="")
tiff(filename = out,units="in",width = 6,height = 7,res=300)
DotPlot(obj.list[1]$Ctrl,features = genes) + ggtitle("Ctrl")  +  theme(plot.title = element_text(hjust = 0.5))
dev.off() 
out <- paste("images\\",id,"_Dot_Psor.tiff",sep="")
tiff(filename = out,units="in",width = 6,height = 7,res=300)
DotPlot(obj.list[2]$Psor,features = genes) + ggtitle("Psor")  +  theme(plot.title = element_text(hjust = 0.5))
dev.off() 
out <- paste("images\\",id,"_Dot_total.tiff",sep="")
tiff(filename = out,units="in",width = 6,height = 7,res=300)
DotPlot(ps_filtr,features = genes) + ggtitle("Total")  +  theme(plot.title = element_text(hjust = 0.5))
dev.off() 
