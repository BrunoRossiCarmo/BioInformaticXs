#Bruno Rossi Carmo
#30/01/2021
#brunorossicarmo@usp.br
#Available in: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7070/files/


#Directory:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\ESETS" #Meu Diretório.
path2 <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Projetos\\Isaac_Succinato\\MicroARRAY_E-MTAB" #Meu Diretório Secundário.
setwd(dir.main)


#Additional:
id <- "E-MTAB-7070" #ID do DataSet no ArrayExpress.


#Libraries:
#Biology:
library(mta10transcriptcluster.db)
library("BiocManager")
library(ArrayExpress)
library(annotate)
library(limma)
library(affy)


#Data_Manipulation:
library(readxl)
library(annotate)
library(plyr)
library(minfi)
library(dplyr)



#Get Raw Set:
rawset <- ArrayExpress(id)
ExpSet <- rma(rawset) #Normalização do arquivos brutos.
expF <- exprs(ExpSet) #Obtenção da matriz de expressão.


#Export ExpProfile:
output <- paste(id,"_raw_expression.tsv")
write.table(expF,file =output,sep = "\t",row.names=T,col.names=NA)


#Check if data is normalized (Boxplot)
png(paste(id,"_boxplot_check.png"),
    width = 5*300,
    height = 5*300,res = 300,pointsize = 8)
boxplot(expF, col = "red")
dev.off()


#Printing densityplot
png(paste(id,"_density_plot.png"),
    width = 5*300,
    height = 5*300,
    res = 300,pointsize = 8)
densityPlot(expF, sampGroups = colnames(expF))
dev.off()


#MDS Plot:
plotMDS(expF, col= c(rep(2,3),rep(1,3)))


#Annotation_Genes_Probes:
x <- mta10transcriptclusterENTREZID #Criar objeto que relaciona probe com EntrezId.
mapped_probes <- mappedkeys(x)
xx <- data.frame(x[mapped_probes])  #Objeto final como data frame.
#Write table about:
write.table(xx$gene_id, file = "gene_id.txt", row.names = F, col.names = F, quote = F)


#Probes_Cluster:
probe_ids <- featureNames(ExpSet)
probes <- data.frame(probe_ids)
colnames(probes) <- "probe_id"
#Write table about:
write.table(probe_ids, file = "probes_id.txt", row.names = F, col.names = F, quote = F)


###############
#DEG_From_Paper:
setwd(path2)
xls <- read_excel("Author_DEG.xlsx")
xls <- data.frame(xls)
#Creating just the Important Data:
xls <- data.frame(xls$ID,xls$Gene.Symbol,xls$Group,xls$Fold.Change)
#Adjust names:
names(xls)[1] <- "probe_id"


###############
#Filter DEG_Paper:
setwd(dir.main)
table(duplicated(xls$probe_id)) #Não tem Probes Duplicadas.
table(duplicated(xls$xls.Gene.Symbol)) #Tem Symbols Duplicados.
table(is.na(xls$probe_id)) #Não tem probes vazias.
table(is.na(xls$xls.Gene.Symbol)) #Tem Symbols Vazios.
xls <- xls[!duplicated(xls$xls.Gene.Symbol),]
xls <- xls[!is.na(xls$xls.Gene.Symbol),]
dim(xls)


#Filtr_By_My_Annots:
#Aqui eu estou usando as anotações obtidas pelo pacote mta10:
#Escolher este ou o de baixo (Filtr_By_Paper_Annots) para prosseguir.
genes <- read.delim("mart_export.txt",header = T)
colnames(genes) <- c("Gene_Stable_ID","Gene_Type","symbol","gene_id")
filtr <- semi_join(probes,xx,by="probe_id")
filtr1 <- inner_join(filtr,xx,by="probe_id")
filtr1$gene_id <- as.integer(filtr1$gene_id)
filtr2 <- inner_join(filtr1,genes,by="gene_id" )
annot_t2 <- filtr2 #Mesma variável.
dim(annot_t2) 
#20210 valores.
#annot_t1 == minhas anotações.


#Filtr_By_Paper_Annots:
#Aqui eu estou usando as anotações obtidas pelo artigo fornecido.
#Escolher este ou o de cima (Filtr_By_My_Annots) para prosseguir.
filtr_t2 <- semi_join(probes,xls,by="probe_id")
filtr_t3 <- inner_join(filtr_t2,xls,by="probe_id")
filtr_t3 <- filtr_t3[!is.na(filtr_t3$xls.Gene.Symbol),]
annot_t2 <- filtr_t3   #Mesma variável.
colnames(annot_t2) <- c("probe_id","symbol","group")
colnames(annot_t2) <- c("probe_id","symbol","type")
dim(annot_t2)
#39685 valores.
#annot_t2 == anotações do artigo.


#PhenoData:
pdat<- pData(ExpSet)
#Write table about:
output <- paste(id,"_metadata.tsv")
write.table(pdat,file =output,sep = "\t",row.names=T,col.names=NA)


#Factors:
pdata_end <- pdat
pdata_end$factors <- "na"
pdata_end$factors[pdata_end$Characteristics.phenotype. == "GCaMP3 positive, Tomato (NaV1.8) negative neurons"] <- "GCaMP3"
pdata_end$factors[pdata_end$Characteristics.phenotype. == "GCaMP3 positive, Tomato (NaV1.8) positive neurons"] <- "Tom"
table(pdata_end$factors)
choices <- factor(pdata_end$factors)


#Selection of probes - PT.1:
#Filtragem das Probeas
mean_all <- apply(expF,1,mean)   #Média de todos os genes,
mean_all  <- mean_all[order(mean_all,decreasing = TRUE)]
Expressed <- mean_all[1:round(0.75*length(mean_all))]
expF_filt_1 <- expF[names(Expressed),]


#Selection of probes - PT.2:
#Filtragem das Probes de baixa expressão.
table(rowMeans(expF_filt_1) >= 1)
table(rowSums(expF_filt_1 >= 1) >= 3)
filter <- rowSums(expF_filt_1 >= 1) >= 3
expF_filt_2 <- expF_filt_1[filter,]


#Test:
#Observar as diferenças de dimensão final.
table(rowMeans(expF_filt_2) >= 1)
table(rowSums(expF_filt_2 >= 1) >= 3)
dim(expF)
dim(expF_filt_1)
dim(expF_filt_2)


#Design Matrix:
design <- model.matrix(~0+choices)
colnames(design) <- substring(colnames(design),8,15)


#Contrast and Fit:
fit <- lmFit(expF_filt_2,design)
cont.matrix_1 <- makeContrasts(GCaMP3vcTOM = GCaMP3 - Tom, levels = design) #Treat vs Control.
#Interpret: Down <- Genes com expressão mais baixa do que expressão no controle.
#Interpret: Up <- Genes com expressão mais alta do que expressão no 


#Next Step: Initial Results:
fit2 <- contrasts.fit(fit, cont.matrix_1)
fit2 <- eBayes(fit2,trend=TRUE)
summary(deg <- decideTests(fit2, p.value=0.05, lfc=0.5))
all <- topTable(fit2, number=nrow(fit2), adjust="BH", coef= "GCaMP3vcTOM", confint = T)
head(all)


#Annotat_Filtr:
table(rownames(all) %in% annot_t2$probe_id)
all2 <- all[rownames(all) %in% annot_t2$probe_id,]


#Annotat_Filtr_Pt.2:
annot1 <- data.frame(Probes = annot_t2$probe_id[annot_t2$probe_id %in% rownames(all2)], Symbol = annot_t2$symbol[annot_t2$probe_id %in% rownames(all2)])


#how many probes/symbols are duplicated?
table(duplicated(annot1$Probes)) 
annot1 <- annot1[!duplicated(annot1$Probes),]
table(rownames(all2) %in% annot1$Probes)
all2 <- all2[match(annot1$Probes,rownames(all2)),]
table(rownames(all2) == annot1$Probes)
all2$symbol <- annot1$Symbol


#Order the tables by symbol and then by P.value
tab1 = all2[order(all2[,'P.Value'],decreasing = FALSE),]


#Remove the duplicated symbols by getting the first occurrence (lowest P-value)
table(duplicated(tab1$symbol))
tab1_unique = tab1[!duplicated(tab1$symbol),]
table(tab1_unique$symbol=="")


#Removing NAs and empty probes
dim(tab1_unique)
tab1_unique <- tab1_unique[!is.na(tab1_unique$symbol),]
dim(tab1_unique)
tab1_unique <- tab1_unique[!tab1_unique$symbol == "",]
dim(tab1_unique)


exp2 <- as.data.frame(expF_filt_2[rownames(tab1_unique),])
dim(exp2)
table(is.na(tab1_unique$symbol))
exp2$symbol <- tab1_unique$symbol
head(exp2)


#Including column DE in deg table
tab1_unique$de <- "no_sig"
tab1_unique$de[tab1_unique$P.Value < 0.05 & tab1_unique$logFC >= 0.5] <- "up"
tab1_unique$de[tab1_unique$P.Value < 0.05 & tab1_unique$logFC <= -0.5] <- "down"


#moving the last column (symbol) to the first position and writing tables
exp2_new <- exp2 %>% dplyr::select(symbol, everything())               #Passar para esquerda.
#Write table about:
f_out <- paste(id,"_expressionCollapsed_filter_br.tsv",sep="")
write.table(exp2_new,file = f_out,sep = "\t",row.names=T,col.names=NA)


tab1_unique_new <- tab1_unique %>% dplyr::select(symbol, everything()) #Passar para esquerda.
head(tab1_unique_new, n=2)
#Write table about:
f_out <- paste(id,"_DEGtable_filter_paper.tsv",sep="")
write.table(tab1_unique_new,file = f_out,sep = "\t",row.names=T,col.names=NA)


#########################################################
#Quantification (PAPER):
idents <- tab1_unique_new[!tab1_unique_new$de=="no_sig",]
#Write table about:
f_out <- paste(id,"_idents_paper.tsv",sep="")
write.table(idents ,file = f_out,sep = "\t",row.names=T,col.names=NA)
dim(idents) #13.663 DEG || Paper Annotation.
idents <- read.table("Data\\TSV\\E-MTAB-7070_idents_paper.tsv",header = T)
rownames(idents) <- idents$X
idents$X <- NULL
#########################################################


#########################################################
#Quantification (MY_ANNOT):
my_idents <- tab1_unique_new[!tab1_unique_new$de=="no_sig",]
#Write table about:
f_out <- paste(id,"_my_idents.tsv",sep="")
write.table(my_idents,file = f_out,sep = "\t",row.names=T,col.names=NA)
dim(my_idents) #6484 DEG || My Annotation.
my_idents <- read.table("Data\\TSV\\E-MTAB-7070_my_idents.tsv",header = T)
rownames(my_idents) <- my_idents$X
my_idents$X <- NULL
#########################################################


#SCATTERPLOT (LogFC):
#Transformar o FC da tabela em LogFC:
xls$signal <- "+"
xls$signal[xls$xls.Fold.Change<0] <- "-"
xls$xls.Fold.Change <- log2(ifelse(xls$xls.Fold.Change<0,xls$xls.Fold.Change*-1,xls$xls.Fold.Change))
xls$xls.Fold.Change[xls$signal=="-"] <- xls$xls.Fold.Change[xls$signal=="-"] * -1
#Reduzir a quantia pelos valores conhecidos:
#Pelas anotações do artigo (symbols genicos):
my_idents$probe_id <- row.names(my_idents)
xls_formated <- inner_join(xls,my_idents,by="probe_id")


#Q-Q PLOT (Scatter):
f_out <- paste(id,"_Q_QPLOT.png",sep="")
png(filename= f_out, width = 800, height =800)
plot(xls_formated$xls.Fold.Change, xls_formated$logFC , main="Q-Q Plot LogFC Comparison",
     xlab = "LogFC DEG Paper",ylab = "My DEG LogFC",pch=19,xlim=c(-10,5),ylim=c(-10,5))
lines(-10:5,-10:5, col = "darkred",lwd=3)
model.regression <- lm(logFC  ~  xls.Fold.Change, data= xls_formated)
abline(model.regression,col="darkblue")
legend(-10.12, 5.3, legend=c("Ideal Relation LogFC Trend", "Real Relation LogFC Trend"),
       col=c("darkred", "darkblue"),lty=1:1, cex=0.8)
dev.off()


#VENN Diagram:
library('VennDiagram')
cross <- semi_join(idents, my_idents, by = "symbol")
f_out <- paste(id,"_Venn_Annots.png",sep="")
png(filename= f_out, width = 800, height =800)
draw.pairwise.venn(length(idents$symbol), length(my_idents$symbol),length(cross$symbol),
                   category = c("Paper Annots","My Annots"),ty = rep("blank",2),
                   fill = c("light blue", "pink"),  cat.fontfamily = "sans", cex = 2, lwd = 4,
                   alpha = rep(0.5, 2), cat.pos = c(0,0), col=c("steelblue1", "maroon"),
                   cat.dist = rep(-0.025, 2),cat.cex=c(2,2),cat.fontface=c(2,2),ext.dist=c(-0.03))
dev.off()


#PLOTS:
#Volcano plot com ggplot2
library(ggplot2)
f_out <- paste(id,"_volcanoplot_TOMvsGC3_paper.png",sep="")
png(filename= f_out, width = 773, height =586)
ggplot(tab1_unique_new, aes(x= logFC, y= -log10(P.Value), shape=de, color=de, size = de)) +
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
    labs(title = "Comparison: TOM vs GC3",x="Log2FoldChange", y = "-log10(P-Value)")
dev.off()


#HeatMap Preparation:
top25 <- head(tab1_unique_new, n=20)
exp_top25 <- exp2_new[row.names(top25),]
row.names(exp_top25) <- exp_top25$symbol
exp_top25 <- exp_top25[,-1]


#HeatMap Adjust:
library(pheatmap)
paletteLength <- 50
myColor <- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
myBreaks <- c(seq(min(scale(t(exp_top25)), na.rm = T), 0, 
                  length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale(t(exp_top25)), na.rm = T)/paletteLength, 
                  max(scale(t(exp_top25)), na.rm = T), 
                  length.out=floor(paletteLength/2)))


#HeatMap Aesthetics:
DF = data.frame(Samples=rep(c("GC3","TOM"),c(3,3)))
rownames(DF) = colnames(exp_top25)


#HeatMap Expression Comparison:
f_out <- paste(id,"_heatmap_TOMvsGC3_paper.png",sep="")
png(filename= f_out, width = 800, height = 600)
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
         annotation_colors = list(Samples = c(GC3 = "gray72" , TOM = "grey0")),
         annotation_legend = T,
         main = "Top 25 E-MTAB-7070",
         fontsize= 15)
dev.off()


#HeatMap Total Preparation:
idents <- tab1_unique_new[!tab1_unique_new$de=="no_sig",]
exp_top <- exp2_new[row.names(idents),]
exp_top$symbol <- NULL
dim(exp_top)


#HeatMap Geral:
f_out <- paste(id,"_heatmap_TOTAL_paper.png",sep="")
png(filename= f_out, width = 1200, height = 1200)
pheatmap(exp_top,
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
         cellheight = 0.07,
         treeheight_row = 120,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_colors = list(Samples = c(GC3 = "gray72" , TOM = "grey0")),
         annotation_legend = T,
         main = "DE Genes E-MTAB-7070",
         fontsize= 20)
dev.off()

