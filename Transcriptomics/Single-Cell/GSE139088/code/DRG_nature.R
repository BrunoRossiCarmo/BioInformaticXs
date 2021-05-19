#Bruno Rossi Carmo
#09/04/21
#brunorossicarmo@usp.br
#Available in: 


#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\sc2_GSE139088" #Meu
setwd(dir.main)
id <- "GSE139088"

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
library("rsvd")
library(Matrix)
library(ggplot2)
library(readxl)


#For memory:
memory.limit(9999999999)


#Carregando os Dados:
rawcounts <- read.csv("code/GSM4130750_WT_1.csv", header = T)
genes <- data.frame(rawcounts[,1])


#Organizando os Dados
rawcounts[2,] <- colnames(rawcounts)


#Metadados (Cell Identity e Barcode)
meta <- t(rawcounts[1:2,])
colnames(meta) <- meta[1,]
meta <- as.data.frame(meta[-1,])


#Contar os "-n" dos Barcodes:
table(sub(".*-", "", meta$`Cell Barcode`))


#Ajustes:
colnames(rawcounts) <- rawcounts[1,]
rawcounts <- rawcounts[- c(1,2),]
rawcounts <- rawcounts[!duplicated(rawcounts[,1]),]
row.names(rawcounts) <- rawcounts[,1]
rawcounts$`Cell Barcode` <- NULL
drg_sparse <- Seurat::as.sparse(rawcounts)


#Save rds:
saveRDS(drg_sparse,"drg_sparse.rds")
saveRDS(rawcounts,"rawcounts.rds")
saveRDS(meta,"meta.rds")


#Rds:
drg_sparse <- readRDS("drg_sparse.rds")
rawcounts <- readRDS("rawcounts.rds")
meta_rds <- readRDS("meta.rds")
class(rawcounts[1,])
barcode <- as.vector(rawcounts[1,])
class(barcode)
cell_ident = as.vector(colnames(rawcounts))


#Quantificando tipos celulares de cada condição:
types = c("Aβ-field-LTMR cells","Aβ-RA-LTMR cells","Aδ-LTMR cells","C-LTMR cells","CGRP-α cells","CGRP-ε cells","CGRP-η cells","CGRP-γ cells","CGRP-θ cells", "CGRP-ζ cells","MRGPRD cells","proprioceptors","SST cells","cold thermoceptor cells")
ncells_E11.5 <- sum(1951+5402+2781)
ncells_E11.5
ncells_E12.5 <- c(30,20,30,122,57,87,48,60,9,37,555,37,24,105)
ncells_E12.5
ncells_E15.5 <- c(61,33,96,383,144,45,26,97,208,63,670,40,61,128)
ncells_E15.5
ncells_p0 <- c(214,163,165,739,284,188,122,216,359,122,1704,103,397,284)
ncells_p0
ncells_p5 <- c(209 ,297,237,1392,445,473,153,334,640,243,3019,104,787,405)
ncells_p5
ncells_adult <- c(257,273,182,1554,1440,850,270,705,758,333,2817,234,761,488)
sum(ncells_adult)

#Fazendo análise de condição:
cell_types <- data.frame(types, ncells_E12.5,ncells_E15.5,ncells_p0,ncells_p5,ncells_adult)
table(gsub("\\.[0-9]*$", "", meta_rds$Cell.identity))
"É adulto!"


#Alterando Metadados:
meta <- meta[gsub("\\.[0-9]*$", "", meta$Cell.identity),]
row.names(meta) <- colnames(rawcounts)
meta[,1] <- NULL

#-------------------------------------------INICIO DA ANALISE-----------------------------------------#
#Iniciando análise:
drg_x <- CreateSeuratObject(counts = rawcounts , project = "drg",
                            min.cells = 0, min.features = 0,
                            meta.data = meta)
?CreateSeuratObject


                               #MT-Percent:
#Mitochondrial percents / Non-neuronal cells markers percents:
drg_x  <- PercentageFeatureSet(drg_x , pattern = "Sucnr1", col.name = "Suc.pc")
max(drg_x@meta.data$Suc.pc)
drg_x  <- PercentageFeatureSet(drg_x , pattern = "^mt", col.name = "percent.mt")


#Data Exploration:
out <- paste("data\\images3\\",id,"_ViolinPlot_PreFilt.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
VlnPlot(drg_x , features = c("Suc.pc"), ncol = 1)
VlnPlot(drg_x , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


#ViolinPlot ir? auxiliar a saber qual ser? o melhor m?todo de cuttoff.
#Usar estat?stica para definir cuttoff [pode ser interessante].
out <- paste("data\\images3\\",id,"_QC_ViolinPlot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
plot1 <- FeatureScatter(drg_x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(drg_x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()



#Filtering by cells and mt.mitochondrial:
drg_filtr <- subset(drg_x, subset = percent.mt < 5 & nCount_RNA < 60000 & nFeature_RNA < 6200)
VlnPlot(drg_filtr , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
drg_filtr <- subset(drg_filtr, Cell.identity != "Unknown")
drg_filtr <- PercentageFeatureSet(drg_filtr , pattern = "Sucnr1", col.name = "Suc.pc")
max(drg_filtr@meta.data$Suc.pc)


#Pos-Filtro.
out <- paste("data\\images3\\",id,"_QC_Filter_ViolinPlot.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
plot1 <- FeatureScatter(drg_filtr, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(drg_filtr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


#ScTransform: Normalization
drg_filtr<- SCTransform(drg_filtr, vars.to.regress = "percent.mt")
saveRDS(drg_filtr, file = "drg_filtr.rds")
drg_filtr <- readRDS("drg_filtr.rds")


#PCA: Principal Component Analysis.
drg_filtr <- RunPCA(drg_filtr, verbose = FALSE)
ElbowPlot(drg_filtr,ndims = 50)
print(drg_filtr[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(drg_filtr, dims = 1:2, reduction = "pca")
DimPlot(drg_filtr, reduction = "pca")
DimHeatmap(drg_filtr, dims = 1:30, cells = 500, balanced = TRUE)


#Create UMAP
#PCS: 35 PCS.
npc <- 1:45
drg_filtr <- FindNeighbors(drg_filtr, dims = npc, verbose = FALSE)
drg_filtr<- FindClusters(drg_filtr, verbose = FALSE, resolution = 0.8)


#Run non-linear dimensional reduction 
#UMAP geralmente fica melhor.
drg_filtr <- RunUMAP(drg_filtr, dims = npc, verbose = FALSE)
drg_filtr <- RunTSNE(drg_filtr, dims = npc, verbose = FALSE)


#Genes:
genes <- data.frame(row.names(drg_filtr))


#Image_UMAP:
out <- paste("data\\images3\\",id,"_Clusters_UMAP.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(drg_filtr, reduction = "umap",pt.size=0.7,label = T)  
dev.off()


#Image_Tsne:
out <- paste("data\\images2\\",id,"_Clusters_TSNE.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DimPlot(drg_filtr, reduction = "tsne",pt.size=0.7,label = T)  
dev.off()


#Agora veremos as anotações:
types = c("Aβ-field-LTMR cells","Aβ-RA-LTMR cells","Aδ-LTMR cells","C-LTMR cells","CGRP-α cells","CGRP-ε cells","CGRP-η cells","CGRP-γ cells","CGRP-θ cells", "","MRGPRD cells","proprioceptors","SST cells","cold thermoceptor cells","Nonpeptidergic.nociceptors")


#Find Markers:
drg.markers <- FindAllMarkers(drg_filtr, only.pos = TRUE, logfc.threshold = 0.25)
drg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- drg.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(data.frame(top10),file = "data\\tsv\\markers.tsv")
saveRDS(drg_filtr, file = "last_drg.rds")
meta <- drg_filtr@meta.data
?FindAllMarkers


#Comparative Markers between annotations:
my_markers <- read.table("data\\tsv\\markers.tsv",header = T, row.names = 1)
author_markers <- read_excel("data\\Markes\\41586_2019_1900_MOESM1_ESM.xlsx",col_names = T)
#
table(meta$Cell.identity[meta$SCT_snn_res.0.8==28])
table(meta$Cell.identity)


#Anotações das Indentidades Celulares:
FeaturePlot(drg_filtr, features = c("Ikzf1","Rgs6","P2rx6"))
VlnPlot(drg_filtr, features = c("Ikzf1","Rgs6","P2rx6"))
"Hipótese: Cluster19,27 == Aβ-field-LTMR cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==19])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==27])


FeaturePlot(drg_filtr, features = c("Gabrg1","Pcdh18","Hs3st4"))
VlnPlot(drg_filtr, features = c("Gabrg1","Pcdh18","Hs3st4"))
"Hipótese: Cluster 17  == Aβ-RA-LTMR cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==17])


FeaturePlot(drg_filtr, features = c("Foxp2","Tox","Trpm8"))
VlnPlot(drg_filtr, features = c("Foxp2","Tox","Trpm8"))
"Hipótese: Cluster 13,26,25 == Cold Thermoceptors Cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==13])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==26])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==25])


FeaturePlot(drg_filtr, features = c("9230105E05Rik","Sst","Ptafr"))
VlnPlot(drg_filtr, features = c("9230105E05Rik","Sst","Ptafr"))
"Hipótese: Cluster 1 == Somatostatin-positive"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==1])


FeaturePlot(drg_filtr, features = c("Nxph1","Runx3","Wnt7a"))
VlnPlot(drg_filtr, features = c("Nxph1","Runx3","Wnt7a"))
"Hipótese: Cluster 18 == Proprioceptors Cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==18])


FeaturePlot(drg_filtr, features = c("Mrgprd","Otoa","Prkcq","Fam167a"))
VlnPlot(drg_filtr, features = c("Mrgprd","Otoa","Prkcq","Fam167a"))
"Hipótese: Cluster 24,3,12,8,9,6 == Nonpeptidergic Nociceptors Cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==24])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==3])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==12])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==8])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==9])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==6])


FeaturePlot(drg_filtr, features = c("Smr2","Igsf21","Chst1"))
VlnPlot(drg_filtr, features = c("Smr2","Igsf21","Chst1"))
"Hipótese: Cluster 16 == CGRP-ζ cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==16])


FeaturePlot(drg_filtr, features = c("Mrgpra3","Adora2b","Hrh1"))
VlnPlot(drg_filtr, features = c("Mrgpra3","Adora2b","Hrh1"))
"Hipótese: Cluster 7,20 == CGRP-θ cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==7])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==20])


FeaturePlot(drg_filtr, features = c("6330403A02Rik","Greb1l","Chrnb4"))
VlnPlot(drg_filtr, features = c("6330403A02Rik","Greb1l","Chrnb4"))
"Hipótese: Cluster 4 == CGRP-γ cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==4])


FeaturePlot(drg_filtr, features = c("Bmpr1b","Tctex1d1","Ngef"))
VlnPlot(drg_filtr, features = c("Bmpr1b","Tctex1d1","Ngef"))
"Hipótese: Cluster  == CGRP-η cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==14])


FeaturePlot(drg_filtr, features = c("Oprk1","Traf3ip3","Tlr4"))
VlnPlot(drg_filtr, features = c("Oprk1","Traf3ip3","Tlr4"))
"Hipótese: Cluster 8,23 == CGRP-ε cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==5])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==23])


FeaturePlot(drg_filtr, features = c("Avpr1a","Slc6a7","Sstr2"))
VlnPlot(drg_filtr, features = c("Avpr1a","Slc6a7","Sstr2"))
"Hipótese: Cluster 0,11  == CGRP-α cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==0])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==11])


FeaturePlot(drg_filtr, features = c("Cacna1i","Nkd2","Ttc29"))
VlnPlot(drg_filtr, features = c("Cacna1i","Nkd2","Ttc29"))
"Hipótese: Cluster 10,2,15 == C-LTMR cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==10])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==2])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==15])


FeaturePlot(drg_filtr, features = c("Adtrp","Colq","Kirrel3os","Cox6a2"))
VlnPlot(drg_filtr, features = c("Adtrp","Colq","Kirrel3os","Cox6a2"))
"Hipótese: Cluster 21 == Aδ-LTMR cells cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==21])


#Unknown:
table(meta$Cell.identity[meta$SCT_snn_res.0.8==20])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==22]) "Hipótese: Cluster 22 == CGRP-β cells"
table(meta$Cell.identity[meta$SCT_snn_res.0.8==26])
table(meta$Cell.identity[meta$SCT_snn_res.0.8==6])


Cluster_Names <- c("CGRP-α1","Somatostatin-positive","C-LTMR1","Nonpeptidergic Nociceptors1",
                   "CGRP-γ","CGRP-ε1","Nonpeptidergic Nociceptors2",
                   "CGRP-θ1","Nonpeptidergic Nociceptors3","Nonpeptidergic Nociceptors4",
                   "C-LTMR2","CGRP-α2","Nonpeptidergic Nociceptors5",
                   "Cold Thermoceptors1","CGRP-η"," C-LTMR3",
                   "CGRP-ζ","Aβ-RA-LTMR","Proprioceptors",
                   "Aβ-field-LTMR1","CGRP-θ2","Aδ-LTMR",
                   "CGRP-β","CGRP-ε2","Nonpeptidergic Nociceptors6",
                   "Cold Thermoceptors2","Cold Thermoceptors3","Aβ-field-LTMR2")
Cluster_Names2 <- c("CGRP","Somatostatin-positive","C-LTMR","Nonpeptidergic Nociceptors",
                   "CGRP","CGRP","Nonpeptidergic Nociceptors",
                   "CGRP","Nonpeptidergic Nociceptors","Nonpeptidergic Nociceptors",
                   "C-LTMR","CGRP","Nonpeptidergic Nociceptors",
                   "Cold Thermoceptors","CGRP"," C-LTMR",
                   "CGRP","Aβ-RA-LTMR","Proprioceptors",
                   "Aβ-field-LTMR","CGRP","Aδ-LTMR",
                   "CGRP","CGRP","Nonpeptidergic Nociceptors",
                   "Cold Thermoceptors","Cold Thermoceptors","Aβ-field-LTMR")

names(Cluster_Names2) <- levels(drg_filtr)
drg_filtr <- RenameIdents(drg_filtr, Cluster_Names2)


#New UMAP:
out <- paste("data\\images3\\",id,"_Together_UMAP.png",sep="")
png(filename = out,units="in",width = 15,height = 7,res=300)
DimPlot(drg_filtr, reduction = "umap",pt.size=0.7,label = T)  
dev.off()
        
        
#Without names:
out <- paste("data\\images3\\",id,"_Ending_nonames_UMAP.png",sep="")
png(filename = out,units="in",width = 15,height = 7,res=300)
DimPlot(drg_filtr, reduction = "umap",pt.size=0.7,label = F)
dev.off()


#Plots:
targets <-c("Hcar1", "Sucnr1")
FeaturePlot(drg_filtr, features = targets)
VlnPlot(drg_filtr, features = targets)
DotPlot(drg_filtr, features = targets)
DoHeatmap(drg_filtr, features = targets) #Ajustar default (GEN) do SCT.

#Sucnr1:
out <- paste("data\\images3\\",id,"_Dot_SUCNR1.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DotPlot(drg_filtr, features = "Sucnr1")
dev.off()
out <- paste("data\\images3\\",id,"_Violin_SUCNR1.png",sep="")
png(filename = out,units="in",width = 20,height = 7,res=300)
VlnPlot(drg_filtr, features = "Sucnr1")
dev.off()
out <- paste("data\\images3\\",id,"_Feature_SUCNR1.png",sep="")
png(filename = out,units="in",width = 20,height = 7,res=300)
FeaturePlot(drg_filtr, features = "Sucnr1")
dev.off()


#Hcar1:
out <- paste("data\\images3\\",id,"_Dot_HCAR1.png",sep="")
png(filename = out,units="in",width = 10,height = 7,res=300)
DotPlot(drg_filtr, features = "Hcar1")
dev.off()
out <- paste("data\\images3\\",id,"_Violin_HCAR1.png",sep="")
png(filename = out,units="in",width = 20,height = 7,res=300)
VlnPlot(drg_filtr, features = "Hcar1")
dev.off()
out <- paste("data\\images3\\",id,"_Feature_HCAR1.png",sep="")
png(filename = out,units="in",width = 20,height = 7,res=300)
FeaturePlot(drg_filtr, features = "Hcar1")
dev.off()


#DEG GENES:
CGRP_gama <- FindMarkers(drg_filtr, ident.1 = "CGRP-γ", only.pos = T, min.pct = 0)
f_out <- paste(id,"_DEGtable_CGRP.tsv",sep="")
write.table(CGRP_gama,file = f_out,sep = "\t",row.names=T,col.names=T)
Hcar1 <- FindMarkers(drg_filtr, ident.1 = c("CGRP-ζ","CGRP-η"), only.pos = T, min.pct = 0)
f_out <- paste(id,"_DEGtable_hcar.tsv",sep="")
write.table(CGRP_gama,file = f_out,sep = "\t",row.names=T,col.names=T)







