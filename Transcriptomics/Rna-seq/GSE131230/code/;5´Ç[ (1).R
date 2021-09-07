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
raw <- read.csv("GSE131230_rlog_official.csv",header = T)


#Take Genes
write.table(raw$X,file = "Genes_Annot.txt",row.names = F, col.names = F, quote = F)
genes <- read.delim("mart_export.txt")
#Tem todos os genes!


#Raw:
raw$Gene.stable.ID <- raw$X
raw$X <- NULL


#Take the genes again:
"Igualar as tables com msm coluna (chave) para realizar o join"
join <- inner_join(raw,genes,by="Gene.stable.ID") #Joining the tables by the ID gene.
join$Gene.description = join$Gene.description = join$Gene.stable.ID = NULL


#Filter 2 for HeatMap:
join <- join[!duplicated(join$Gene.name),] #Remove ID genes duplicated.
row.names(join) <- join$Gene.name
join$Gene.name <-NULL


#Library and other:
library(pheatmap)
paletteLength <- 50
myColor <- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
myBreaks <- c(seq(min(scale(t(join)), na.rm = T), 0,length.out=ceiling(paletteLength/2) + 1), 
              seq(max(scale(t(join)), na.rm = T)/paletteLength,max(scale(t(join)),na.rm = T),
                  length.out=floor(paletteLength/2)))


#Aesthetics (Sample Legend):
DF = data.frame(Samples=rep(c("Non.Peptidergic Nocioceptor","Peptidergic Nocioceptor","CLTMR",
                              "A-LTMR","A.RA-LTMR","A.SA-LTMR","A.Field-LTMR","Proprioceptors"),c(3,3,3,3,5,3,3,3)))
rownames(DF) = colnames(join)


#Choose the genes:
test <- subset(join, row.names(join) == "Sucnr1" | row.names(join) == "Hcar1" | row.names(join) == "Hcar2")

#Genes Heatmap:
f_out <- paste(id,"_heatmap_GENES.png",sep="")
png(filename= f_out, width = 1200, height = 1200)
pheatmap(test,
         cluster_rows = T,
         cluster_cols = F,
         cutree_cols = 3,
         show_rownames = T,
         show_colnames = F,
         angle_col =90,
         legend = T,
         scale = "row",
         color = myColor,
         breaks =myBreaks,
         gaps_col =,
         cellwidth = 23,
         cellheight = 62,
         treeheight_row = 100,
         treeheight_col = 30,
         border_color = NA,
         display_numbers = F,
         fontsize_row = 10,
         fontsize_col = 10,
         annotation =,
         annotation_col = DF,
         annotation_row = NA,
         annotation_legend = T,
         main = "DE Genes GSE131230",
         fontsize= 20)
dev.off()
