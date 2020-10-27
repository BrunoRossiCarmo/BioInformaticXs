#Bruno Rossi Carmo
#27/10/2020
#brunorossicarmo@usp.br

#|--------------------------First Analyze--------------------------|
#Diretórios:
dir.main <- "C:your_pathway\\GSE54514" #Coloque seu diretorio do arquivo (Separe com "\\")
setwd(dir.main)

#Bibliotecas:
#Obs:Para instalar alguma, use "install.packages("name")
#Bio:
library("BiocManager")
library(GEOquery)
library(limma)
library(affy)
#Plots:
library("gplots")
library("pheatmap")
library(ggplot2)
#Dataset Manipulation:
library(plyr)
library(dplyr)
#Others:
library(annotate)
library(minfi)

#Carregando o arquivo de Series contendo as samples
GSE_1 <- getGEO(filename = "GSE54514_series_matrix.txt.gz") 

#Matriz normalizada
eset <- exprs(GSE_1)

#Raw Expression file:
write.table(eset,file = "raw_expression_set.tsv", sep = "\t",row.names=T,col.names=NA)

#Checar Normalização:
png(file = "raw_exp_Set.PNG", width = 5*300,height = 5*300,res = 300,pointsize = 8)
boxplot(eset, col = "red")
dev.off() #Fecha documento.

#DensityPlot para checar novamente:
png(file = "raw_dens_exp.PNG", width = 5*300,height = 5*300,res = 300,pointsize = 8)
densityPlot(eset, sampGroups = colnames(eset))
dev.off() #Fecha documento.

#Anotações:
probe_ids <- featureNames(GSE_1)
write.table(probe_ids, file = "probes_id.txt", row.names = F, col.names = F, quote = F)

#PhenoData
samp <- pData(GSE_1)
write.table(samp, file = "metadados.tsv",sep = "\t",row.names=F,col.names=T)

#Dados de Expressão e Design Matrix:
choices <- factor(samp$`group_day:ch1`)
design <- model.matrix(~0+choices)
colnames(design) <- levels(choices)



#--------------------------------------Comparação 1:--------------------------------------------#
#Fit:
fit <- lmFit(eset,design)
contrasts <- makeContrasts(NS_D1-S_D1, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)

#Data:
summary(decideTests(fit2, p.value =0.05,lfc = 1))
results <- decideTests(fit2, p.value =0.05,lfc =1)
topTable(fit2, coef=1,confint = T)
all_data<-topTable(fit2, number=nrow(fit2), adjust="BH", coef= 1, confint = T)

#Print about PlotMD:
png(file= "plotMD_NS_D1_Vs_S_D1.PNG",width = 850, height = 642)
plotMD(fit2, coef=1, status=results)
dev.off()

#Annotation File 1:
genes_annot <- read.delim("annotations_genes_probes.txt")
all_data_p2 <- all_data[rownames(all_data) %in% genes_annot$ILLUMINA.HumanHT.12.V3.probe ,]
table(duplicated(genes_annot$ILLUMINA.HumanHT.12.V3.probe))                                        #Infs duplicadas.
genes_annot <- genes_annot[!duplicated(genes_annot$ILLUMINA.HumanHT.12.V3.probe),]                 #Retira duplicadas.
all_data_p2 <- all_data_p2[match(genes_annot$ILLUMINA.HumanHT.12.V3.probe,rownames(all_data_p2)),] #Juntar Infs.
all_data_p2$symbol <- genes_annot$Gene.name

#Ordenar tabelas:
tab1 = all_data_p2[order(all_data_p2[,'P.Value'],decreasing = FALSE),]   #Ordenar tabela por ordem de Symbol e P.Value.
table(duplicated(tab1$symbol))                                           #Ver duplicados.
tab1_unique = tab1[!duplicated(tab1$symbol),]                            #Retirar Simbolos Duplicados.
tab1_unique <- tab1_unique[!is.na(tab1_unique$symbol),]                  #Retirar NA.
tab1_unique <- tab1_unique[!tab1_unique$symbol == "",]                   #Retirar valores vazios.
eset_v2 <- as.data.frame(eset[rownames(tab1_unique),])
eset_v2$symbol <- tab1_unique$symbol

#Inclusão de Tabela:
tab1_unique$de <- "no_sig"                                                        #Informação de Regulação.
tab1_unique$de[tab1_unique$P.Value < 0.05 & tab1_unique$logFC >= 0.5] <- "up"     #Genes Up Regulated (logFC > 0.5).
tab1_unique$de[tab1_unique$P.Value < 0.05 & tab1_unique$logFC <= -0.5] <- "down"  #Genes Down Regulated (logFC <= -0.5).

#Ajuste e Montagem: (Obs: "%>% - Pipe: Possibilita que o argumento da esquerda seja o primeiro parametro da direita.)
eset_v2_new <- eset_v2 %>% select(symbol, everything())
write.table(eset_v2_new,file = "Optimized_Express.tsv",sep = "\t",row.names=T,col.names=NA)
tab1_unique_new <- tab1_unique %>% select(symbol, everything())
write.table(tab1_unique_new,file = "All_Data_Reg.tsv",sep = "\t",row.names=T,col.names=NA)

#VolcanoPlot 1:
png(filename= "Volcano_NS_D1_vs_S_D1.PNG", width = 850, height = 642)
ggplot(tab1_unique_new, aes(x= logFC, y= -log10(P.Value), shape=de, color=de, size = de)) +
        geom_point() +
        scale_shape_manual(values=c(18,20,18))+ 
        scale_color_manual(values=c("blue","black","Red"))+
        #scale_color_manual(values=c("black","black","black"))+
        scale_size_manual(values=c(2,1,2))+
        theme_get()+ 
        geom_hline(yintercept = -log10(.05), color= "black")+
        geom_vline(xintercept = c(-1,1), color= "black")+
        geom_vline(xintercept = c(0), linetype="dotted", color= "black", size= 0.8)+
        labs(title = paste0("Comparison: NS_D1 vs S_D1"),x="Log2FoldChange", y = "-log10(P-Value)")
dev.off()
       

#Heatmap Preparation:
top10 <- head(tab1_unique_new, n=10)
eset_top10 <- eset_v2_new[row.names(top10),]
row.names(eset_top10) <- eset_top10$symbol
eset_top10 <- eset_top10[,-1]                #Torna os simbolos index e some ID-ILUMINA.

#Heatmap_1:
paletteLength <- 50
myColor <- colorRampPalette(c("navyblue", "white", "firebrick3"))(paletteLength)
myBreaks <- c(seq(min(scale(t(eset_top10)), na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(scale(t(eset_top10)), na.rm = T)/paletteLength, max(scale(t(eset_top10)), na.rm = T), length.out=floor(paletteLength/2)))

#Plots:
pheatmap(eset_top10,
         cluster_rows = T,cluster_cols = T,
         show_rownames = T,
         angle_col =90,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col =,
         cellwidth = 15,
         cellheight = 12,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 5,
         fontsize_col = 8,
         annotation =,
         annotation_colors =,
         main = "Top 10 NS_D1 vs S_D1 - GSE54514)")



