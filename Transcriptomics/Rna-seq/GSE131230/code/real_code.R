#Bruno Rossi Carmo
#024/05/21
#brunorossicarmo@usp.br
#Available in: 


#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\RNA\\GSE131230"
id <- "GSE131230"
setwd(dir.main)

#Biology:
library("BiocManager")
library(edgeR)


#Data_Manipulation:
library(dplyr)
library(tidyr)
library(purrr)

#Taking data:
#Aqui estou lendo a matriz de contagem.
expS <- read.csv("GSE131230_counts_official.csv",header = T,row.names = 1)
genes <- read.delim("mart_export.txt")
expS$Gene.stable.ID <- row.names(expS) 
exps <- data.frame(inner_join(expS,genes,by="Gene.stable.ID"))
exps$Gene.description = exps$Gene.stable.ID = NULL
row.names(exps) <- exps$Gene.name
exps <- exps[!duplicated(exps$Gene.name),] #Remove duplicated 
row.names(exps) <- exps$Gene.name
exps$Gene.name <- NULL
colnames(exps)
exps<- exps %>% dplyr::select(c(AÎ..LTMR_1,AÎ..LTMR_2,AÎ..LTMR_3), everything())#Passar o argumento que queremos como primeiro fator
#teste: 
genes_exp <- subset(exps, row.names(tab) == "Sucnr1" | row.names(tab) == "Ntrk2")
write.csv(data.frame(genes_exp),file = "files\\Raw_Expression.csv",sep = "\t",row.names=T,col.names=T)

#groups:
groups = factor(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,5),rep(6,3),rep(7,3),rep(8,3)))
"1-AI-LTMR/2-Nociceptor/3-Peptidergic/4-CLTMR/5-AI^2RA-LTMR/
6-AI^2SA-LTMR/7-AI-FIELD.LTMR/8-Proprioceptor"

#Filtering:
y <- DGEList(counts = as.matrix(exps),group = groups, remove.zeros = F)


#Normalization:
y <- calcNormFactors(y,method="TMM")  #Method = TMM by default (just few genes DE).
out <- paste("images\\",id,"_batcht_effects.png",sep="")
png(filename = out,units="in",width = 10, height = 7,res=300)
plotMDS(y, col=c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,5),rep(6,3),rep(7,3),rep(8,3)))
dev.off()

#Design:
design <- model.matrix(~groups)

#Dispersion:
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
sqrt(y$common.dispersion) #BCV - Biological Coefficient of Variation.
out <- paste("images\\",id,"_bcv_plot.png",sep="")
png(filename = out,units="in",width = 10, height = 7,res=300)
plotBCV(y)
dev.off()

#Fit/Contrast:
fit <- glmQLFit(y, design, robust=TRUE)
out <- paste("images\\",id,"_disp_plot.png",sep="")
png(filename = out,units="in",width = 10, height = 7,res=300)
plotQLDisp(fit)
dev.off()

#Realizar Comparação:
test <- glmQLFTest(fit, coef = 2:8)
topTags(test,n=30,adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
summary(decideTests(test))
FDR <- p.adjust(test$table$PValue, method="BH")

#Argumentação DEG:
tab <- data.frame(test$table)
tab$FDR <- FDR
tab$de <- "no_sign"
tab$de[tab$FDR < 0.053 & (tab$logFC.groups2 >= 0.5 & tab$logFC.groups3 >= 0.5 & tab$logFC.groups4 >= 0.5 & tab$logFC.groups5 >= 0.5 & tab$logFC.groups6 >= 0.5 & tab$logFC.groups7 >= 0.5 & tab$logFC.groups8 >= 0.5)] <- "up"
tab$de[tab$FDR < 0.053 & (tab$logFC.groups2 <= 0.5 & tab$logFC.groups3 <= 0.5 & tab$logFC.groups4 <= 0.5 & tab$logFC.groups5 <= 0.5 & tab$logFC.groups6 <= 0.5 & tab$logFC.groups7 <= 0.5 & tab$logFC.groups8 <= 0.5)] <- "down"
t
genes_op <- subset(tab, row.names(tab) == "Sucnr1" | row.names(tab) == "Ntrk2")

#Top_Up and Top_Down:

#UP:
top_up <- tab[tab$de=="up",]
mediaFC_up <- (top_up$logFC.groups2 + top_up$logFC.groups3 + top_up$logFC.groups4 + top_up$logFC.groups5 + top_up$logFC.groups6 + top_up$logFC.groups7 + top_up$logFC.groups8)/7
top_up$Media_FC <- mediaFC_up #Media dos FoldChange para classificação.
top_up <- top_up[order(top_up$Media_FC, decreasing = T),]
tup <- top_up
tup[,1:12] <- NULL
write.table(tup,file="files\\DOWN_GENES.tsv", sep = "\t",col.names= NA)

#DOWN:
top_down <- tab[tab$de=="down",]
mediaFC_down <- (top_down$logFC.groups2 + top_down$logFC.groups3 + top_down$logFC.groups4 + top_down$logFC.groups5 + top_down$logFC.groups6 + top_down$logFC.groups7 + top_down$logFC.groups8)/7
top_down$Media_FC <- mediaFC_down #Media dos FoldChange para classificação.
top_down <- top_down[order(top_down$Media_FC, decreasing = F),]
points <- c(rownames(top_down[1:10,]),rownames(top_up[1:10,]),"Sucnr1","Ntrk2")
write.csv(data.frame(tab[points,]),file = "files\\DEG_Points.csv",sep = "\t",row.names=T,col.names=T)

#topdown:
tdown <- top_down
tdown[1:12] <- NULL
tdown$Media_FC <- gsub("-","",tdown$Media_FC)
tdown <- subset(tdown, tdown$Media_FC >= 4.5)
write.table(tdown,file="files\\UP_GENES.tsv", sep = "\t",col.names= NA)

#Construir nova tabela de contagens normalizadas:
logcpm <- cpm(y, log=TRUE)
logcpm <- data.frame(logcpm)
tab_end <- logcpm
tab_end$gene_name <- row.names(tab_end)
sum(exps)/8

#Construir tabela final:
test <- subset(tab_end , row.names(tab_end) == "Sucnr1" | row.names(tab_end) == "Ntrk2" | row.names(tab_end) == "Hcar2")
tab_end_res <- tab_end[points,]
tab_end_res$gene_name <- NULL


#Library and other:
library(pheatmap)
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "grey99", "red"))(paletteLength)
"No MyBreaks, primeiro é realizado a operação de transpor a matriz de 
contagem com os genes específicos, realizando o scale para cada gene. Após isto,
é realizado a sequência do menor scalling até o valor determinado, no primeiro caso
é 0. Ou seja, serão as N/2 saídas pela paletteLenght desde o mínimo até o 0, depois do 0 até 
o máximo do scalling."
myBreaks <- c(seq(min(scale(t(tab_end_res)), na.rm = T), 0,length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale(t(tab_end_res)), na.rm = T)/paletteLength,max(scale(t(tab_end_res)),na.rm = T),
                  length.out=floor(paletteLength/2)))

#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("Aδ-LTMR","Peptidergic","C-LTMR","Nociceptor","Aβ RA-LTM","Aβ SA1-LTMR","Aβ Field-LTMR","Proprioceptor"),c(3,3,3,3,5,3,3,3)))
rownames(DF) = colnames(tab_end_res)

#Colors for Annotation:
ann_colors = list(
        Samples = c("Aδ-LTMR" = "chartreuse","Peptidergic" =" chocolate","C-LTMR" ="cornflowerblue","Nociceptor" ="darkmagenta", "Aβ RA-LTM" ="firebrick4","Aβ Field-LTMR"="gold3","Proprioceptor" ="deeppink1","Aβ SA1-LTMR" ="darkgray"))

#Tab-End Res:
tab_end_res
vls <- subset(tab_end_res, row.names(tab_end_res) == "Sucnr1")
mean_values <- c(as.numeric((vls[1] +vls[2]+vls[3])/3),as.numeric((vls[4] +vls[5]+vls[6])/3),
                 as.numeric((vls[7] +vls[8]+vls[9])/3), as.numeric((vls[10] +vls[11]+vls[12])/3),
                 as.numeric((vls[13] +vls[14]+vls[15]+vls[16]+vls[17])/5), as.numeric((vls[18] +vls[19]+vls[20])/3), as.numeric((vls[21] +vls[22]+vls[23])/3), as.numeric((vls[24] +vls[25]+vls[26])/3))
"8 1 7 5 4 3 6 2"

#Re-order:
tab_end_nv <- c(tab_end_res[4:6],tab_end_res[24:26],tab_end_res[18:20],tab_end_res[13:17],
                tab_end_res[10:12],tab_end_res[21:23],tab_end_res[7:9],tab_end_res[1:3])
tab_end_nv <- data.frame(tab_end_nv)

#Re-order 2:
DF = data.frame(Samples=rep(c("Nociceptor","Proprioceptor","Aβ SA1-LTMR","Aβ RA-LTMR","C-LTMR","Aβ Field-LTMR","Peptidergic","Aδ-LTMR"),c(3,3,3,5,3,3,3,3)))
rownames(DF) = colnames(tab_end_nv)
myBreaks <- c(seq(min(scale(t(tab_end_nv)), na.rm = T), 0,length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale(t(tab_end_nv)), na.rm = T)/paletteLength,max(scale(t(tab_end_nv)),na.rm = T),
                  length.out=floor(paletteLength/2)))

#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("Nociceptor","Proprioceptor","Aβ SA1-LTMR","Aβ RA-LTMR","C-LTMR","Aβ Field-LTMR","Peptidergic","Aδ-LTMR"),c(3,3,3,5,3,3,3,3)))
rownames(DF) = colnames(tab_end_nv )

#Colors for Annotation:
ann_colors = list(
        Samples = c("Nociceptor" = "darkorchid4","Proprioceptor" =" blue4","Aβ SA1-LTMR" ="deepskyblue","Aβ RA-LTMR" ="forestgreen", "C-LTMR" ="gold","Aβ Field-LTMR"="darkorange1","Peptidergic" ="darkgoldenrod4","Aδ-LTMR" ="firebrick4"))
row.names(tab_end_nv) <- row.names(tab_end_res)

#Heatmap Test:
out <- paste("images\\",id,"heatmap.png",sep="")
png(filename = out,units="in",width = 12, height = 7,res=300)
pheatmap(tab_end_nv,
         cluster_rows = F,
         cluster_cols = F,
         cutree_rows = c(10,20),
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col = c(3,6,9,14,17,20,23,26),
         cellwidth = 12,
         cellheight = 12,
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
         gaps_row = c(20),
         main = "DE Genes GSE131230",
         fontsize= 20,
         annotation_colors = ann_colors)
dev.off()

#Second_Plot:
vls <- subset(tab_end_nv, row.names(tab_end_nv) == "Sucnr1")
mean_values <- c(as.numeric((vls[1] +vls[2]+vls[3])/3),as.numeric((vls[4] +vls[5]+vls[6])/3),
                 as.numeric((vls[7] +vls[8]+vls[9])/3), as.numeric((vls[10] +vls[11]+vls[12])/3),
                 as.numeric((vls[13] +vls[14]+vls[15]+vls[16]+vls[17])/5), as.numeric((vls[18] +vls[19]+vls[20])/3), as.numeric((vls[21] +vls[22]+vls[23])/3), as.numeric((vls[24] +vls[25]+vls[26])/3))
mean_values <- 2^(mean_values)
mean_values <- data.frame(mean_values)
row.names(mean_values) <- c("Nociceptor","Proprioceptor","Aβ SA1-LTMR","Aβ RA-LTMR","C-LTMR","Aβ Field-LTMR","Peptidergic","Aδ-LTMR")
out <- paste("images\\",id,"barplot.png",sep="")
png(filename = out,units="in",width = 12, height = 7,res=300)
barplot(mean_values$mean_values,name = row.names(mean_values), col = "red",main = "Distribuição de Sucnr1 em Identidades Neuronais", ylab = "Contagens de Sucnr1")
dev.off()

#Heatmap Up:
points <- c(rownames(top_up[1:10,]),"Sucnr1","Ntrk2")
tab_end$gene_name <- row.names(tab_end)
tab_end_res <- tab_end[points,]
tab_end_res$gene_name <- NULL

#Heatmap Test:
out <- paste("images\\",id,"heatmap_down2.png",sep="")
png(filename = out,units="in",width = 12, height = 7,res=300)
pheatmap(tab_end_res,
         cluster_rows = F,
         cluster_cols = F,
         cutree_rows = c(10,20),
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col = c(3,6,9,12,17,20,23,26),
         cellwidth = 12,
         cellheight = 15,
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
         gaps_row = c(10),
         main = "DE Genes Down GSE131230",
         fontsize= 20)
dev.off()

#Heatmap Up:
points <- c(rownames(top_down[1:10,]),"Sucnr1","Ntrk2")
tab_end$gene_name <- row.names(tab_end)
tab_end_res <- tab_end[points,]
tab_end_res$gene_name <- NULL

#Heatmap Test:
out <- paste("images\\",id,"heatmap_up2.png",sep="")
png(filename = out,units="in",width = 12, height = 7,res=300)
pheatmap(tab_end_res,
         cluster_rows = F,
         cluster_cols = F,
         cutree_rows = c(10,20),
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col = c(3,6,9,12,17,20,23,26),
         cellwidth = 12,
         cellheight = 15,
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
         gaps_row = c(10),
         main = "DE Genes Up GSE131230",
         fontsize= 20)
dev.off()
