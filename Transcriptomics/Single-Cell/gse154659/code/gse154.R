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

#Creating Seurat Object 1:
drg_ATK <- CreateSeuratObject(counts = ATK , project = "ATK",min.features = 500)

#Organization aTF:
drg_ATK@meta.data$type <- "Atf"
drg_ATK@meta.data$model <- ""
drg_ATK@meta.data$model[grep("Naive",row.names(drg_ATK@meta.data))] <- "Naive"
drg_ATK@meta.data$model[grep("Crush",row.names(drg_ATK@meta.data))] <- "Crush"
drg_ATK@meta.data$trat <- ""
drg_ATK@meta.data$trat[grep("KO",row.names(drg_ATK@meta.data))] <- "KO"
drg_ATK@meta.data$trat[grep("WT",row.names(drg_ATK@meta.data))] <- "WT"
View(drg_ATK@meta.data)

#Creating Seurat Object 2:
drg_C57 <- CreateSeuratObject(counts = C57 , project = "C57",min.features = 500)

#Organization aTF:
drg_ATK@meta.data$type <- "Atf"
drg_ATK@meta.data$model <- ""
drg_ATK@meta.data$model[grep("Naive",row.names(drg_ATK@meta.data))] <- "Naive"
drg_ATK@meta.data$model[grep("Crush",row.names(drg_ATK@meta.data))] <- "Crush"
drg_ATK@meta.data$trat <- ""
drg_ATK@meta.data$trat[grep("KO",row.names(drg_ATK@meta.data))] <- "KO"
drg_ATK@meta.data$trat[grep("WT",row.names(drg_ATK@meta.data))] <- "WT"
View(drg_ATK@meta.data)
