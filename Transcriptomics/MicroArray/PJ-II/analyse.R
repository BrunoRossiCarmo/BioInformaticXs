##Bruno Rossi Carmo
##20/10/2021
#gseGSE161361
#Diretorio:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Grad\\6_sem\\PJII" #Coloque seu diretorio
setwd(dir.main)
getwd()

#Bibliotecas:
#Bio:
library("BiocManager")
library(GEOquery)
library(affy)
library(limma)
#Plots:
library("gplots")
library("pheatmap")
library(ggplot2)
#Dataset Manipulation:
library(plyr)
library(dplyr)

#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#Data:
gse <- getGEO("GSE161361") 
eset <- exprs(gse[[1]])
eset <- data.frame(eset)

#Table
table <- read.csv("table\\data.csv",header = T)

#Adjust Table:
table(remove) <- is.na(table$X)
eset <- subset(x = table, X != "NA")
row.names(eset) <- eset$X
eset$X <- NULL

#Renames:
names(eset) <- c(c(paste0("OP",seq(1,3))) ,  c(paste0("WT",seq(1,3))))

#Verify Normalization Log:
min(eset)
max(eset)
boxplot(eset, col = "red")

#Create PNG's:
png(paste0(dir.main, "\\images\\_boxplot_check.png"),
    width = 1500,
    height = 1500,
    res = 300,
    pointsize = 8)
boxplot(eset, col = "red",main="Norm Check")
dev.off()

#Normalization:
eset <- log2(eset)
png(paste0(dir.main, "\\images\\_boxplot_norm.png"),
    width = 1500,
    height = 1500,
    res = 300,
    pointsize = 8)
boxplot(eset, col = "red",main="Normalized")
dev.off()

#Verification:
min(eset)
max(eset)
head(summary(eset))

#Table Normalized:
eset <- read.csv("table/normalized.csv", row.names=1)

#MDS Plot:
png(paste0(dir.main, "\\images\\_mdsPLOT.png"),
    width = 1500,
    height = 1500,
    res = 300,
    pointsize = 8)
plotMDS(eset, col= c(rep(2,3),rep(1,3)))
dev.off()

#Making the Data:
annot <- fData(gse[[1]])[,c('ID','geneSymbol','EntrezID')]
annot <- subset(annot, ID != "NA" & geneSymbol != "")
annot <- read.csv("table/feature_data.csv", row.names=1)

#Reading the data (if necessary):
eset2 <- eset
eset2$ID <- row.names(eset2)

#New Eset:
library(dplyr)
feature_eset <- data.frame(inner_join(eset2,annot,by="ID"))
row.names(feature_eset) <- feature_eset$ID
feature_eset$ID <- NULL
f_out <- paste(dir.main,"\\table\\metadata.csv",sep="")
write.csv(feature_eset,file = f_out,row.names=T)
feature_eset <- read.csv("table/metadata.csv", row.names=1)

#Filtros Gerais:
dim(feature_eset) #7472 circRNA
feature_eset <- subset(feature_eset, !grepl('///', geneSymbol))
dim(feature_eset) #7472 circRNA
feature_eset <- subset(feature_eset, geneSymbol != "NA")
dim(feature_eset) #7472 circRNA
feature_eset <- subset(feature_eset, geneSymbol != "-")
dim(feature_eset) #7223 circRNA

#Filtragem das Probes de baixa express?o.
new_eset <- feature_eset
new_eset$EntrezID = NULL
#table(rowMeans(new_eset) >= 1)
#table(rowSums(new_eset >= 1) >= 2)
#filter <- rowSums(new_eset >= 1) >= 3
#new_eset <- new_eset[filter,]

#Criação dos fatores:
choices <- factor(c(  c(rep("OP",3))  ,  c(rep("WT",3)  )))
choices

#Desging Matrix:
design <- model.matrix(~0+choices)
colnames(design) <- levels(choices)

#Contrast and Fit:
fit <- lmFit(new_eset,design)
contrasts <- makeContrasts(OP-WT, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)

#Data:
summary(decideTests(fit2, p.value =0.05,lfc = 0.5))
results <- decideTests(fit2, p.value =0.05,lfc =0.5)
topTable(fit2, coef=1,confint = T)
all_data<-topTable(fit2, number=nrow(fit2), adjust="BH", coef= 1, confint = T)
table(table(row.names(all_data)))

#Ordenar tabelas:
all_data <- all_data[order(all_data[,'P.Value'],decreasing = FALSE),]   #Ordenar tabela P-Valor. 

#DEG Counts:
all_data$de <- "no_sig"                                                       
all_data$de[all_data$P.Value < 0.05 & all_data$logFC >= 0.5] <- "up"     
all_data$de[all_data$P.Value < 0.05 & all_data$logFC <= -0.5] <- "down"  
table(all_data$de)

#VolcanoPlot 1:
library(ggplot2)
all_data_new <- subset(all_data, P.Value > (10**(-10)))
png(filename= "images/Volcano_OP_vs_WT.PNG", width = 773, height =586)
ggplot(all_data_new, aes(x= logFC, y= -log10(P.Value), shape=de, color=de, size = de)) +
    geom_point() +
    scale_shape_manual(values=c(20,20,20))+ 
    scale_color_manual(values=c("cornflowerblue","black","Red"))+
    scale_size_manual(values=c(2,1,2))+
    theme_get()+
    theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+
    geom_hline(yintercept = -log10(.05), color= "black",size=1)+
    geom_vline(xintercept = c(-0.5,0.5), color= "black",size=1)+
    geom_vline(xintercept = c(0), linetype="dotted", color= "black", size= 0.8)+
    geom_vline(xintercept = c(-4,-2,2,4), linetype="dotted", color= "grey", size= 0.2)+
    labs(title = "Comparison: OP vs WT",x="Log2FoldChange", y = "-log10(P-Value)")
dev.off()


#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#Heatmap Preparation:
top50 <- subset(all_data,de=="up")
top50 <- top50[order(top50$logFC,decreasing = T),]
top50 <- head(top50,n=50)
eset_top50 <- new_eset[row.names(top50),]
geneS <- eset_top50$geneSymbol
eset_top50$geneSymbol = NULL

#Scale:
paletteLength <- 50
myColor <- colorRampPalette(c("mediumblue", "white", "red2"))(paletteLength)
myBreaks <- c(seq(min(scale(t(eset_top50)), na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(scale(t(eset_top50)), na.rm = T)/paletteLength, max(scale(t(eset_top50)), na.rm = T), length.out=floor(paletteLength/2)))

#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("OP","WT"),c(3,3)))
rownames(DF) = colnames(eset_top50)

#Colors for Annotation:
ann_colors = list(
    Samples = c("OP" = "black","WT" =" dimgray"))

#Plots:
png(paste0(dir.main, "\\images\\_HEAT.png"),
    width = 3000,
    height = 4000,
    res = 300,
    pointsize = 8)
pheatmap(eset_top50,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col = c(3,6),
         cellwidth = 37,
         cellheight = 17,
         treeheight_row = 50,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_legend = T,
         main = "",
         fontsize= 18,
         annotation_colors = ann_colors)
dev.off()

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
#Heatmap Preparation:
top50 <- subset(all_data,de=="up")
top50 <- top50[order(top50$logFC,decreasing = T),]
top50 <- head(top50,n=50)
eset_top50 <- new_eset[row.names(top50),]
row.names(eset_top50) <- geneS

#Heatmap_1:
paletteLength <- 50
myColor <- colorRampPalette(c("mediumblue", "white", "red2"))(paletteLength)
myBreaks <- c(seq(min(scale(t(eset_top50)), na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(scale(t(eset_top50)), na.rm = T)/paletteLength, max(scale(t(eset_top50)), na.rm = T), length.out=floor(paletteLength/2)))

#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("OP","WT"),c(3,3)))
rownames(DF) = colnames(eset_top50)

#Colors for Annotation:
ann_colors = list(
    Samples = c("OP" = "black","WT" =" dimgray"))

#Plots:
png(paste0(dir.main, "\\images\\_HEAT_SYMBOL.png"),
    width = 3000,
    height = 4000,
    res = 300,
    pointsize = 8)
pheatmap(eset_top50,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col = c(3,6),
         cellwidth = 37,
         cellheight = 17,
         treeheight_row = 50,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_legend = T,
         main = "",
         fontsize= 18,
         annotation_colors = ann_colors)
dev.off()

