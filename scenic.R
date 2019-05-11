#scenic eg. E14_mRPC data

setwd('/data2/hehan/Rstudio/mRPC/E14_RPC/E14_seurat_hh/data_190220/scenic/')
dir.create("int")
dir.create('data')

options(stringsAsFactors = F)

#
load(file = 'Rstudio/mRPC/E14_RPC/E14_seurat_hh/data_190220/E14_umap_190220.RData')
data <- as.matrix(pbmc_E14@data)

#save anno and colvars
anno <- as.matrix(pbmc_E14@meta.data[,'res.0.6'])
rownames(anno) <- rownames(pbmc_E14@meta.data)

cellInfo<-anno
col_hh <- c('0'='#1f77b4','1'='#ff7f0e','2'='#2ca02c','3'='#d62728','4'='#9467bd',
            '5'='#8c564b','6'='#e377c2','7'='#bcbd22','8'='#17becf','9'= '#aec7e8')
colVars<-list(cluster=setNames(as.character(col_hh),names(col_hh)))
saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="int/colVars.Rds")

## filter gene and save
library(SCENIC)
library(methods)
org="mgi" # or hgnc, or dmel
dbDir="databases" # RcisTarget databases location
myDatasetTitle=analysis_name # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=ncore) 
scenicOptions@settings$dbs<-"/data2/hehan/scripts/scenic/database/mm9-tss-centered-5kb-7species.mc9nr.feather"
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#expression(>0)at least 5 cells
genesLeft_minCells<-rownames(data)[rowSums(data>0)>5]

library(RcisTarget)
motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]

saveRDS(genesLeft_minCells_inDatabases, file=getIntName(scenicOptions, "genesKept"))

exprMat_filtered <- data[genesLeft_minCells_inDatabases, ]
write.csv(file="data/filter_data.csv",t(exprMat_filtered))##用于pyscenic

##jupyter & terminal : Inference of co-expression modules

# source activate python36
# nohup pyscenic grn -o /data2/hehan//co_expression.csv --num_workers 10 /data2/hehan/Rstudio//data/filter_data.csv /data2/hehan/scripts/scenic/TF_data/tf_mus.txt > scenic.log &
# source deactivate
 
#save the co_expression
coexpression<-read.csv("data/co_expression.csv",sep="\t")
colnames(coexpression)<-c("TF","Target","weight")
coexpression$weight<-coexpression$weight/100
saveRDS(file="int/1.4_GENIE3_linkList.Rds",coexpression)
#
corrMat <- cor(t(exprMat_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
#
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, data,skipTsne=T)
runSCENIC_4_aucell_binarize(scenicOptions,skipTsne=T,skipHeatmaps = T)

date()
#sessionInfo()

