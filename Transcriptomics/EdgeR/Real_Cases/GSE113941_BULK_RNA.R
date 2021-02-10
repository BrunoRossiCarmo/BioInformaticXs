#Bruno Rossi Carmo
#02/02/2021
#brunorossicarmo@usp.br
#Available in: https: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113941

#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\GSE113941" #Meu Diretório.
setwd(dir.main)


#Additional:
id <- "GSE113941" #ID do DataSet no GEO.


#Libraries:
#Biology:
library("BiocManager")
library(edgeR)


#Data_Manipulation:
library(dplyr)
library(tidyr)
library(purrr)


#Taking data:
raw <- read.delim("GSE113941_uniform_process_tpm_vals.txt",row.names = 1)
expS <- data.frame(raw$DRG_normal_TRAP1._TPM,raw$DRG_normal_TRAP2_TPM,raw$DRG_normal_TRAP3_TPM,
                   raw$DRG_normal_TRAP4_TPM,raw$DRG_CIPN_TRAP5_TPM,raw$DRG_CIPN_TRAP6_TPM,
                   raw$MERGED_DRG_CIPN_TRAP7_TPM,raw$DRG_CIPN_TRAP8_TPM)
row.names(expS) <- row.names(raw)
names(expS)


#Filter (CPM):
#Tirar genes de baixa expressão.
keep <- rowSums(cpm(expS)>2) >=4  #CPM = Counts p/ million.
expS_filtr <- expS[keep,]
dim(expS)
dim(expS_filtr)


#Experimental Design:
choices <- factor(c(rep(1,4),rep(2,4))) #1 = Normal | 2 = CIPN.


#DGEList:
y <- DGEList(counts = expS_filtr, group = choices, remove.zeros = TRUE)


#Data Exploration 1:
pseudoCounts <- log2(y$counts + 1)
head(pseudoCounts)
par(mar = c(8,4,1,2))
boxplot(pseudoCounts, col = "red", las = 3, cex.names = 1)
f_out <- paste(id,"_normalization_boxplot.png",sep="")
#Make Image:
png(filename= f_out, width = 773, height =586)
boxplot(cpm(y,log = T), col = "red")
dev.off()


#Normalization:
y <- calcNormFactors(y,method="TMM")  #Method = TMM by default (just few genes DE).
colnames(y) <- c("DRG_normal_Trap1","DRG_normal_Trap2","DRG_normal_Trap3",
                 "DRG_normal_Trap4","DRG_CIPN_Trap1","DRG_CIPN_Trap2",
                 "DRG_CIPN_Trap3","DRG_CIPN_Trap4")


#Data Exploration:
f_out <- paste(id,"_batcht_effects.png",sep="")
#Make Image:
png(filename= f_out, width = 773, height =586)
plotMDS(y, col=rep(1:4,each=4))
dev.off()


#Design:
design <- model.matrix(~0+choices)
colnames(design) <- c("DRG_Norm_Trap","DRG_CIPN_Trap")
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
sqrt(y$common.dispersion) #BCV.
f_out <- paste(id,"_bcv_plot.png",sep="")
#Make Image:
png(filename= f_out, width = 773, height =586)
plotBCV(y)
dev.off()


#Fit/Contrast:
fit <- glmQLFit(y, design, robust=TRUE)
f_out <- paste(id,"_disp_plot.png",sep="")
#Make Image:
png(filename= f_out, width = 773, height =586)
plotQLDisp(fit)
dev.off()
contrast <- makeContrasts(G2vsG1 = DRG_CIPN_Trap - DRG_Norm_Trap, levels = design)  #CIPN trat vs Normal Trat.


#DE: Differential Expression:
qlf <- glmQLFTest(fit, contrast = contrast)
topTags(qlf,n=10,adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
FDR <- p.adjust(qlf$table$PValue, method="BH")
summary(decideTests(qlf))
top <- rownames(topTags(qlf))
cpm(y)[top,]


#MD Plot Average:
plotMD(qlf)
abline(h=c(-1,1), col="blue")


#Annotate:
tab1 <- data.frame(qlf$table)
tab1 = tab1[order(tab1[,'PValue'],decreasing = FALSE),]
genes_id <- row.names(tab1) 
write.table(genes_id,file = "Genes_Annot.txt",row.names = F, col.names = F, quote = F)
genes <- read.delim("mart_export.txt") #Pegar anotações direto do Ensembl.


#Annotate Pt.2:
tab1$Gene.stable.ID.version <- genes_id 
tab_end <- inner_join(tab1,genes,by="Gene.stable.ID.version") #Joining the tables by the ID gene.
dim(tab_end)


#Filtration:
tab_end <- tab_end[!duplicated(tab_end$Gene.stable.ID.version),] #Remove duplicated genes.
dim(tab_end)
row.names(tab_end) <- tab_end$Gene.stable.ID.version #Pass the ID gene to the row names.
tab_end$Gene.stable.ID.version <- NULL #Exclude de column ID gene.


#DEG Table:
#Modelo usado: pvalue < 0.5 e logFC entre -1.5 a 1.5.
tab_end$de <- "no_sig"
tab_end$de[tab_end$PValue < 0.05 & tab_end$logFC >= 1.5] <- "up"
tab_end$de[tab_end$PValue < 0.05 & tab_end$logFC <= -1.5] <- "down"


#Adjust Symbols Column.
tab1_unique_new <- tab_end %>% dplyr::select(Gene.name, everything())
f_out <- paste(id,"_DEGtable.tsv",sep="")
#Write table:
write.table(tab1_unique_new,file = f_out,sep = "\t",row.names=T,col.names=NA)
symbols <- row.names(tab1_unique_new) 
#Write table:
write.table(symbols,file = "symbols.txt",row.names = F, col.names = F, quote = F)


#Volcano plot com ggplot2
library(ggplot2)
f_out <- paste(id,"_volcanoplot_CIPNvsNormal_V2_pv05.png",sep="")
png(filename= f_out, width = 773, height =586)
ggplot(tab1_unique_new, aes(x= logFC, y= -log10(PValue), shape=de, color=de, size = de)) +
        geom_point() +
        scale_shape_manual(values=c(20,20,20))+ 
        scale_color_manual(values=c("cornflowerblue","black","Red"))+
        scale_size_manual(values=c(2,1,2))+
        theme_get()+
        theme(axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"))+
        geom_hline(yintercept = -log10(.05), color= "black",size=1)+
        geom_vline(xintercept = c(-1.5,1.5), color= "black",size=1)+
        geom_vline(xintercept = c(0), linetype="dotted", color= "black", size= 0.8)+
        geom_vline(xintercept = c(-4,-2,2,4), linetype="dotted", color= "grey", size= 0.2)+
        labs(title = "Comparison: DRG_CIPN vs DRG_normal",x="Log2FoldChange", y = "-log10(P-Value)")
dev.off()


#Filter for HeatMap:
#LogCPM = matrix of log2counts-per-million.
logcpm <- cpm(y, log=TRUE) #Creating the filtrated data to heatmap.
logcpm <- data.frame(logcpm)
logcpm$Gene.stable.ID.version <- row.names(logcpm)
logcpm <- inner_join(logcpm,genes,by="Gene.stable.ID.version") #Joining the tables by Gene ID.
dim(logcpm)


#Filter 2 for HeatMap:
logcpm <- logcpm[!duplicated(logcpm$Gene.stable.ID.version),] #Remove ID genes duplicated.
dim(logcpm)
row.names(logcpm) <- logcpm$Gene.stable.ID.version


#Exclusion of unnecessary column:
logcpm$Gene.stable.ID.version <- NULL
logcpm <- logcpm %>% dplyr::select(Gene.name, everything())
logcpm$Gene.type <- NULL
logcpm$Gene.description <- NULL
logcpm$NCBI.gene..formerly.Entrezgene..ID <- NULL
row.names(logcpm) <- logcpm$Gene.name
logcpm$Gene.name <- NULL


#Heatmap Preparation:
row.names(tab1_unique_new) <- tab1_unique_new$Gene.name
top25 <- head(tab1_unique_new, n=25)
exp_top25 <- logcpm[row.names(top25),]
colnames(exp_top25) <- c("DRG_Normal.1","DRG_Normal.2",
                         "DRG_Normal.3","DRG_Normal.4",
                         "DRG_CIPN.1","DRG_CIPN.2",
                         "DRG_CIPN.3","DRG_CIPN.4")


#Library and other:
library(pheatmap)
paletteLength <- 50
myColor <- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
myBreaks <- c(seq(min(scale(t(exp_top25)), na.rm = T), 0, 
                  length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale(t(exp_top25)), na.rm = T)/paletteLength, 
                  max(scale(t(exp_top25)), na.rm = T), 
                  length.out=floor(paletteLength/2)))


#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("DRG_Normal","DRG_CIPN"),c(4,4)))
rownames(DF) = colnames(exp_top25)


#Expression Comparison:
f_out <- paste(id,"_heatmap_toP25.png",sep="")
png(filename= f_out, width = 1200, height = 800)
pheatmap(exp_top25,
         cluster_rows = T,
         cluster_cols = T,
         cutree_cols = 2,
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col =,
         cellwidth = 30,
         cellheight = 20,
         treeheight_row = 50,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_colors = list(Samples = c(DRG_Normal = "gray72" , DRG_CIPN = "grey0")),
         annotation_legend = T,
         main = "Top 25 GSE113941",
         fontsize= 15)
dev.off()


#Seps (Tirar aqueles que não são DEG):
idents <- tab1_unique_new[!tab1_unique_new$de=="no_sig",]


#Quantifications:
dim(idents) #DE genes quantification (505 Genes).
exp_geral <- logcpm[row.names(idents),]


#Geral HeatMap (All DEG):
f_out <- paste(id,"_heatmap_TOTAL.png",sep="")
png(filename= f_out, width = 1200, height = 1200)
pheatmap(exp_geral,
         cluster_rows = T,
         cluster_cols = T,
         cutree_cols = 2,
         show_rownames = F,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col =,
         cellwidth = 48,
         cellheight = 2,
         treeheight_row = 100,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_colors = list(Samples = c(DRG_Normal = "gray72" , DRG_CIPN = "grey0")),
         annotation_legend = T,
         main = "DE Genes GSE113941",
         fontsize= 20)
dev.off()


#Separation between "up" and "down" for GO and Reactome:
idents_2 <- idents[!idents$de=="down",]
idents_3 <- idents[!idents$de=="up",]
genes_idents <- idents_2$Gene.name
genes_idents2 <- idents_3$Gene.name
write.table(genes_idents ,file = "Genes_Up.txt",row.names = F, col.names = F, quote = F)
write.table(genes_idents2 ,file = "Genes_Down.txt",row.names = F, col.names = F, quote = F)


#TCA GENES (REACTOME):
#Genes From BioPlanet 2019:
tca <- read.delim("Data\\Reactome\\BioPlanet\\BioPlanet_2019_table_up.txt",header = T)
Genes_TCA<-tca[2,]$Genes 
Genes_TCA <- strsplit(Genes_TCA,";")
Genes_TCA <- unlist(Genes_TCA)
Genes_TCA <- data.frame(Genes_TCA)
names(Genes_TCA) <- "symbol"


#Genes from Reactome 2016:
tca2 <- read.delim("Data\\Reactome\\Reactome_2016\\Reactome_2016_table_up.txt",header = T)
Genes_TCA2<-tca2[2,]$Genes 
Genes_TCA2 <- strsplit(Genes_TCA2,";")
Genes_TCA2 <- unlist(Genes_TCA2)
Genes_TCA2 <- data.frame(Genes_TCA2)
names(Genes_TCA2) <- "symbol"


#Joining all genes:
Genes_TCA_total <- full_join(Genes_TCA,Genes_TCA2, by= "symbol")
#Mesmo sendo poucos genes, bom deixar esse filtro apenas para casos futuros com
#maiores quantias de genes diferencialmente expressos.
Genes_TCA_total <- Genes_TCA_total[!duplicated(Genes_TCA_total$symbol),]
Genes_TCA_total <- data.frame(Genes_TCA_total)
row.names(Genes_TCA_total) <- Genes_TCA_total$Genes_TCA_total


#Take the genes Expression:
logcpm_tca_test <- logcpm
row.names(logcpm_tca_test) <- toupper(row.names(logcpm)) 
exp_tca <- logcpm_tca_test[row.names(Genes_TCA_total),]


#Aesthetics:
DF = data.frame(Samples=rep(c("DRG_Normal","DRG_CIPN"),c(4,4)))
rownames(DF) = colnames(exp_tca)


#Heatmap for TCA Genes:
f_out <- paste(id,"_TCA_DEG_Total.png",sep="")
png(filename= f_out, width = 1200, height = 1200)
pheatmap(exp_tca,
         cluster_rows = T,
         cluster_cols = T,
         cutree_cols = 2,
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col =,
         cellwidth = 60,
         cellheight = 70,
         treeheight_row = 100,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 20,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_colors = list(Samples = c(DRG_Normal = "gray72" , DRG_CIPN = "grey0")),
         annotation_legend = T,
         main = "TCA Genes Identity GSE113941",
         fontsize= 20)
dev.off()

