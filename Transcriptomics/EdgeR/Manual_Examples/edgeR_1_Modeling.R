#Bruno Rossi Carmo
#27/10/2020
#brunorossicarmo@usp.br
#Original Experiment: https://bioconductor.org/packages/release/bioc/html/edgeR.html

#Diretórios:
dir.main <- "C:\\Users\\bruno\\Desktop\\Usp\\Laboratorio\\CRID\\Codigos\\EDGER\\Modelos" 
setwd(dir.main)

#Bibliotecas:
library("BiocManager")
library("edgeR")

#Raw Data (Experimental):
raw <- readRDS("arab.rds")
head(raw) #Temos dois fatores experimentais: 1-hrcc vc mock | 2-Tempo de cada replicada.

#Design dos Fatores Experimentais:
Trat <- factor(substring(colnames(raw),1,4)) #Substring permite retirar determinados caracteres de uma string.
Trat <- relevel(Trat, ref = "mock")
Tempo <- factor(substring(colnames(raw),5,5))

#Criar o DGEObject:
y <- DGEList(counts = raw, group = Trat)

#Filtragem e Normalização:
keep <- filterByExpr(y)   #Essa função determina quais genes tem counts suficientes para análise estatística
table(keep)
y <- y[keep,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)   #Converte raw library sizes em efetivas library sizes.

#Data Exploration:
png(file = "MDS_filter_EDGER.PNG", width = 5*300,height = 5*300,res = 300,pointsize = 8)
plotMDS(y, col=rep(1:2,each=3))
dev.off()

#Análise de Correlações Experimentais:
design <- model.matrix(~Tempo+Tempo:Trat)
logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
cor(logFC[,4:6])  #Observar as correlações entre os tempos, usando T1 como base.
#Ou seja, T2 tem menor correlação, necessitando fazer correção de batch effect.
#Para isso, é adicionado na formula de design: ~Batch_Effect + Other_Treat.

#Design Matrix:
design <- model.matrix(~Tempo+Trat)   
rownames(design) <- colnames(y)
#Intercept = Temp1. 
#Para comparar control com treat é necessário usar coef 4 (última coluna).

#Estimando Dispersão:
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion #Raiz quadrada desse valor é equivalente a BCV = Coeficient of Biological Variation.

png(file = "BCV_filter_EDGER.PNG", width = 5*350,height = 5*300,res = 300,pointsize = 8)
plotBCV(y)
dev.off()

fit <- glmQLFit(y, design, robust=TRUE)
png(file = "glmFIT_filter_EDGER.PNG", width = 5*350,height = 5*300,res = 300,pointsize = 10)
plotQLDisp(fit)
dev.off()

#DE: Differential Expression:
#Vamos utilizar o QL F-Test para significantes genes.
#Primeiro iremos comparar ControlT2 cm T1 e ControlT3 cm T1 para ver se tem diferença no logFC entre eles por conta do batch effect.
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)

#Real Differential Expression
qlf <- glmQLFTest(fit)  #Não define o coef pois por default compara último coef.
topTags(qlf)
top <- rownames(topTags(qlf))
cpm(y)[top,]
summary(decideTests(qlf))

png(file = "DE_exp1_EDGER.PNG", width = 5*350,height = 5*300,res = 300,pointsize = 10)
plotMD(qlf)
abline(h=c(-1,1), col="blue")
dev.off()

#End of the Example.
