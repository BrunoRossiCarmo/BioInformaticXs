#Bruno Rossi Carmo
#28/10/2020
#brunorossicarmo@usp.br

#Estudo oferecido em: (https://www.bioconductor.org/packages/release/bioc/html/MetaVolcanoR.html).
#ESTUDO NÃO ORIGINAL:
#Prada C, Lima D, Nakaya H (2020). MetaVolcanoR: Gene Expression Meta-analysis Visualization Tool. R package version 1.2.0.# 
#|--------------------------First Analyze--------------------------|
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Codigos\\LIMMA\\GSE33099" #Coloque seu diretorio do arquivo
setwd(dir.main)

#Bibliotecas:
#Bio:
library("BiocManager")
library(GEOquery)
library(ExpressionAtlas)
library(limma)

#---------------------Inicio-Expressão-Diferencial---------------------:

#Carregando o arquivo de Series contendo as samples:
gse30999 <- getGEO(filename = "GSE30999_series_matrix.txt.gz") 


#Matriz de Expressão:
eset <- exprs(gse30999)

#Dados Gerais:
samp_annot <- pData(gse30999)

#Obter anotação das Probes:
annottab <- fData(gse30999)[,c('ID', 'Gene Symbol')]    
colnames(annottab) <- c('Probe', 'Symbol')

#Filtrar probes com mais de um símbolo gênico:
annottab <- subset(annottab, !grepl('///', Symbol))      #Filtro: Probes com + de 1 gene que usam '///' para separar symbol.
annottab <- subset(annottab, Symbol != '')               #Filtro: Probes sem Anotação.    
annottab <- annottab[complete.cases(annottab),]          #Filtro: Probes com NA.
annotatb <- annottab[annottab$Probe %in% rownames(eset),] #Tira probes ausentes na matriz de expressão.
genes2probes <- split(annottab$Probe, annottab$Symbol)

#Selecionar probes com maior expressão:
mean_probes <- apply(eset, 1, mean)                      #Calcular média de Expression Set por linhas.         
selected_probes <- lapply(genes2probes, function(probes){
        expression_levels <- mean_probes[probes]
        sel_probe <- probes[which.max(expression_levels)]
        return(sel_probe)
})
selected_probes <- unlist(selected_probes)

#Matriz de Expressão Filtrada e com Gene Symbol:
exp_mat <- eset[selected_probes,]
rownames(exp_mat) <- names(selected_probes)

#Dados de Expressão e Design Matrix:
experimentalFactor <- samp_annot$`biopsy type:ch1`
experimentalFactor <- factor(experimentalFactor, levels = c('non-lesion (NL)', 'psoriasis lesion (LS)'))
design <- model.matrix(~experimentalFactor)
colnames(design) <- c('normal', 'lesional_vs_normal')

#Fit:
fit <- lmFit(exp_mat, design)
fit <- eBayes(fit)
diffexp <- topTable(fit, coef = 'lesional_vs_normal', confint = T)
diffexp <- data.frame(Symbol = rownames(diffexp),log2FC = diffexp$logFC,pvalue = diffexp$P.Value,CI.L = diffexp$CI.L,CI.R = diffexp$CI.R)
all_data <- topTable(fit, number=nrow(fit), adjust="BH", coef= 1, confint = T)

#Resultados:
head(diffexp)
head(all_data)
