#Bruno Rossi Carmo
#28/10/2020
#brunorossicarmo@usp.br

#Estudo oferecido em: (https://www.bioconductor.org/packages/release/bioc/html/MetaVolcanoR.html).
#ESTUDO NÃO ORIGINAL:
#Prada C, Lima D, Nakaya H (2020). MetaVolcanoR: Gene Expression Meta-analysis Visualization Tool. R package version 1.2.0.#

dir.main <- "C:\\Users\\" #Coloque seu diretorio do arquivo
setwd(dir.main)

#Bibliotecas:
#Bio:
library("BiocManager")
library(GEOquery)
library(ExpressionAtlas)
library(limma)
library(MetaVolcanoR)

#Dados de Expressão (Resultado Primeiro Code)
diffexplist <- list(diffexp)          #Transformar em Lista.
class(diffexplist)

#REM: Random Effect Model
meta_degs_rem <- rem_mv(diffexp=diffexplist,
                       pcriteria="pvalue",
                       foldchangecol='log2FC', 
                       genenamecol='Symbol',
                       geneidcol=NULL,
                       collaps=FALSE,
                       llcol='CI.L',
                       rlcol='CI.R',
                       vcol=NULL, 
                       cvar=TRUE,
                       metathr=0.01,
                       jobname="MetaVolcano",
                       outputfolder=".", 
                       draw='HTML',
                       ncores=1)

#Print REM (Variação dos Genes)
png(file= "plot_meta_degs_rem.png",width = 850, height = 642)
meta_degs_rem@MetaVolcano
dev.off()

#DrawForest (HTML/Variação do Gene)
draw_forest(remres=meta_degs_rem,
            gene="MMP9",
            genecol="Symbol", 
            foldchangecol="Log2FC",
            llcol="CI.L", 
            rlcol="CI.R",
            jobname="MetaVolcano",
            outputfolder=".",
            draw="HTML")

#Vote-counting approach 
meta_degs_vote <- votecount_mv(diffexp=diffexplist,
                               pcriteria='pvalue',
                               foldchangecol='Log2FC',
                               genenamecol='Symbol',
                               geneidcol=NULL,
                               pvalue=0.05,
                               foldchange=0, 
                               metathr=0.01,
                               collaps=FALSE,
                               jobname="MetaVolcano", 
                               outputfolder=".",
                               draw='HTML')

#Diferential Expression Plot (HTML)
meta_degs_vote@degfreq

#Combining-approach 
meta_degs_comb <- combining_mv(diffexp=diffexplist,
                               pcriteria='pvalue', 
                               foldchangecol='Log2FC',
                               genenamecol='Symbol',
                               geneidcol=NULL,
                               metafc='Mean',
                               metathr=0.01, 
                               collaps=TRUE,
                               jobname="MetaVolcano",
                               outputfolder=".",
                               draw='HTML')

#MetaVolcano Original
meta_degs_comb@MetaVolcano
