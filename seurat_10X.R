rm(list = ls())
options(stringsAsFactors = F)

##library
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(readxl)

##colors
pal_rainbow <- c("red","orange","yellow","green","cyan","blue","purple")
bertie.color <- c("#D0B38A", "#A3D171", "#533E88", "#7957A3", "#000000", 
                  "#E63325", "#0A8041", "#C5208E", "#3CBDED", "#3B55A3", 
                  "#D691BE", "#D23E28", "#6474B6", "#4288C8", "#80A469", 
                  "#FFCF3F", "#FBCC9F", "#F06360", "#437BB2", "#A43E40", 
                  "#206767", "#779E43", "#258950", "#F7D059", "#ED803B") ##25
pal_dolphin <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", 
                 "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", 
                 "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF", "grey")##15
pal_location <- brewer.pal(12, name = "Paired")
pal_stage <- c("ffd2a5", "d38cad", "8a79af")

col_hh_20 <- c('#9DC8C8','#D1B6E1','#519D9E','#E53A40',
            '#30A9DE','#F17F42','#566270','#8CD790',
            '#2EC4B6', '#44633F','#D81159','#4F86C6',
            '#AACD6E','#F16B6F','#E3E36A','#D8E6E7',
            '#D09E88','#7200da','#f100e5','#005f6b')##20
col_gyd_19 <- c('dodgerblue2','red3','green3','slateblue','darkorange',
                'skyblue1','violetred4','forestgreen','steelblue4','slategrey',
                'brown','darkseagreen','darkgoldenrod3','olivedrab','royalblue',
                'tomato4','cyan2','springgreen2','lightyellow')
col_gyd_12 <- c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
                '#8c564b','#e377c2','#bcbd22','#17becf','#aec7e8',
                '#ffbb78','#98df8a')
col_hh_10 <- c('#9DC8C8','#D1B6E1',"#30A9DE","#F17F42",'#566270',
               "#8CD790","#44633F","#F16B6F","#E3E36A","#f100e5") 
col_hj_55 <- c('#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', 
               '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6', 
               '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87', 
               '#5A0007', '#809693', '#FEFFE6', '#1B4400', '#4FC601', 
               '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', 
               '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA', 
               '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', 
               '#300018', '#0AA6D8', '#013349', '#00846F', '#372101', 
               '#FFB500', '#C2FFED', '#A079BF', '#CC0744', '#C0B9B2',
               '#C2FF99', '#001E09', '#00489C', '#6F0062', '#0CBD66', 
               '#EEC3FF', '#456D75', '#B77B68', '#7A87A1', '#788D66')
col_16 <- c('0' = '#9DC8C8', '1' = '#D1B6E1', '2' = '#519D9E', '3' = '#E53A40',
            '4' = '#30A9DE', '5' = '#F17F42', '6' = '#73d2f3', '7' = '#8CD790',
            '8' = '#2EC4B6', '9' = '#837dff','10' = '#c54fa7','11' = '#4F86C6',
            '12' = '#AACD6E','13'= '#F16B6F','14' = '#E3E36A','15' = '#D09E88')##16
col_heatmap_hh <- c('#ff7761','#ecfafb','#274555') 
col_heatmap_hj <- c("#2E294E","#9055A2","#D499B9","#FFFFF2","#FDD692","#EC7357")

col_stage_feature <- c("grey", "#f26d5b","#492540") 
col_stage_feature <- c('#E0E3DA','#e4406f')

#引入细胞周期 marker（g1$s,g2$m）
#
cell_cycle <- as.data.frame(read_excel("~/Rstudio/cell cycle/cell cycle list.xlsx",sheet = 1))
#segregate this list into markers of G1/S phase, G2/M phase and markers of S phase
g1s.genes <- as.character(na.omit(cell_cycle$`G1/S`),na.omit(cell_cycle$S))
g2m.genes <- as.character(na.omit(cell_cycle$`G2/M`,na.omit(cell_cycle$M)))
#更改大小写(human不用此步)
g1s.genes <- tolower(g1s.genes)
g2m.genes <- tolower(g2m.genes)
#首字母大写(human不用此步)
library(Hmisc)
g1s.genes <- capitalize(g1s.genes)
g2m.genes <- capitalize(g2m.genes)

##mouse
geneset <- list()
geneset$g1s <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
                 "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2",
                 "Rpa2","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7",
                 "Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1",
                 "Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1",
                 "Chaf1b","Brip1","E2f8")

geneset$g2m <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80",
                 "Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a",
                 "Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
                 "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk",
                 "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
                 "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
                 "Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")

geneset$g1s <- toupper(geneset$g1s)
geneset$g2m <- toupper(geneset$g2m)

#读入10X/mm10
data <- Read10X(data.dir="~/Rstudio/...")
data <- as.matrix(data)
#过滤 gene 表达为0的 gene
data <- data[rowMeans(data)>0,]
pbmc@raw.data <- pbmc@raw.data[rowMeans(as.matrix(pbmc@raw.data))>0,]##根据情况决定是否在创建 seurat 对象之前过滤
pbmc@data <- pbmc@data[rowMeans(as.matrix(pbmc@data))>0,]
# filt_cells
filt_cells <- read.csv(file = '',header = F)
#
pbmc <- CreateSeuratObject(raw.data = data,min.cells = 0,min.genes = -1)
# DoubletDetect 标记 doublet
pbmc@meta.data$dblets_1 <- filt_cells$V1
pbmc@meta.data$dblets_1 <- ifelse(pbmc@meta.data$dblets_1 == 1,'doublet','singlet')
doublets_1 <- rownames(pbmc@meta.data[pbmc@meta.data$dblets_1 == 'doublet',])
doublets_2 <- colnames(pbmc@raw.data)[apply(pbmc@raw.data[c('Kdm5d','Eif2s3y','Gm29650','Uty','Ddx3y'),],2,function(x) any(x>0)) & pbmc@raw.data["Xist",]>0]
pbmc@meta.data$dblets_2 <- ifelse(rownames(pbmc@meta.data) %in% doublets_2,'doublet','singlet')
doublets <- union(doublets_1,doublets_2)
pbmc@meta.data$doublets <- ifelse(rownames(pbmc@meta.data) %in% doublets,'doublet','singlet')
                                            
# 标记样本类型
pbmc@meta.data$type <- 'type'

#par(mfrow = c(1,1))#设置窗口#（针对 Error:figure margins too large）
plot(pbmc@meta.data$nUMI, pbmc@meta.data$nGene, 
     xlab = "nUMI - the number of transcripts", 
     ylab = "nGene - the number of genes")

#判断可能会过滤细胞类型，是否为 needed 细胞             
pbmc@meta.data$filt_cells_nGene <- ifelse(pbmc@meta.data$nGene < 2500,"filt_cells_nGene","others")
pbmc@meta.data$filt_cells_nUMI <- ifelse(pbmc@meta.data$nUMI < 10000,"filt_cells_nUMI","others")

## evaluate threshold for percent.mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,group.by = 'dblets_1')
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,group.by = 'dblets_2')
      
#unfilt cells_PCA 查看细胞判断其低质量 or 独立群体细胞             
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- ScaleData(pbmc)
#(对于 scale.factor 的选择：应选择最接近其 nUMI 中位数的数量级)，normalize 后的文件覆盖原 data 表达矩阵；
#seurat中Normolization（值为 log（2）与手动取 log，是不一样的，手动取 log 可能会在后续寻找 df markers 中其 logFC 值可能会被放大；）
#vars.to.regress = c("nUMI" , "percent.mito")#(vars.to.regress,对 UMI 范围不是特别大，不用执行回归) 
             
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 30, do.print = TRUE, 
               pcs.print = 1:10, genes.print = 20)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "filt_cells_nGene")
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "filt_cells_nUMI")
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "doublets")       
FeaturePlot(pbmc, features.plot = c("GYPA", "PTPRE", "GATA1", "KLF1"), 
            reduction.use = "pca", cols.use = c("grey","red"))###红细胞(GYPA:CD235a)
pbmc <- RunUMAP(pbmc,min_dist = 0.2, dims.use = 1:10)
# filter cells
##确定主成分明确不需要的群体：如根据 PC1 >20过滤细胞，及mito >5%
##PC
#pbmc@meta.data$PC1 <- pbmc@dr$pca@cell.embeddings[,1]
#View(pbmc_raw@meta.data)             
#pbmc <- FilterCells(object = pbmc, subset.names = "PC1", high.thresholds = 20)
                                            
pbmc <- FilterCells(object = pbmc, subset.names = 'nGene',low.thresholds = 1000)             
pbmc <- FilterCells(object = pbmc, subset.names = "percent.mito", high.thresholds = 0.03)
pbmc <- SubsetData(pbmc,cells.use = setdiff(colnames(pbmc@data),doublets))###去掉 doublets
ncol(pbmc@data)

## save filt_cells
filt_cells <- colnames(pbmc@data)
filt_cells_raw <- pbmc@raw.data[,filt_cells]
write.csv(as.matrix(filt_cells_raw), file = "~/Rstudio/.../filt_cells_raw.csv")  #此为保存细胞过滤后未normalise表达矩阵

## filt_cells, Detection of variable genes across the single cells
#Assign Cell-Cycle Scores
pbmc <- CellCycleScoring(object = pbmc,  s.genes = geneset$g1s, 
                         g2m.genes = geneset$g2m, set.ident = TRUE)
head(x = pbmc@meta.data)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "Phase")                                           
pbmc <- ScaleData(pbmc,vars.to.regress = c('nUMI','S.Score','G2M.Score'))##根据情况选择是否需要回归周期
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)##10X 数据建议选择2000Gene左右
length(x = pbmc@var.genes)

#hvg <- rownames(pbmc@hvg.info[order(pbmc@hvg.info[,3],decreasing = TRUE),])[1:2000]
library(org.Mm.eg.db)
gogenes <- select(org.Mm.eg.db, keys = c("GO:0007049"), columns = c("SYMBOL"), keytype = "GOALL")
hvg <- setdiff(pbmc@var.genes,unique(gogenes$SYMBOL))##2018
hvg <- setdiff(hvg,c('Ung','Hmmr'))##2016

##                                            
pbmc <- RunPCA(pbmc, pc.genes = hvg, pcs.compute = 50, do.print = TRUE, 
                        pcs.print = 1:20, genes.print = 10) 
PCElbowPlot(object = pbmc, num.pc = 30)
                                            
#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 50, do.print = TRUE, 
               pcs.print = 1:20, genes.print = 10)  #（需查看主成分中是否有细胞周期相关基因）
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2) 


# Visualize the distribution of cell cycle markers across
#RidgePlot(object = pbmc, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), nCol = 2)

#Determine statistically significant principal components
                                        
#pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = pbmc, PCs = 1:50)

PCElbowPlot(object = pbmc, num.pc = 50)

#save PCA    
pdf(file = "~/Rstudio/.../pca/pc1:6.pdf", width = 9.81, height = 6.49)
PCHeatmap(object = pbmc, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

#Cluster the cells 
##根据 pca 选择1:10,or 1:15,分别分群及 runtsne
                                        
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = , 
                     resolution = seq(0.2,2,0.2), print.output = 0, save.SNN = TRUE)
       
##RunUMAP
pbmc <- RunUMAP(pbmc, dims.use = 1:10, min_dist = 0.2, n_neighbors = )
DimPlot(pbmc, reduction.use = 'umap', do.return = T,cols.use = , group.by = 'res.0.6', pt.size = ,
        plot.title = 'pc1:10_umap_E14')
DimPlot(pbmc, reduction.use = 'umap', do.return = T,cols.use = , group.by = 'Phase', pt.size = 1.2,
        plot.title = 'pc1:15_umap_E14')
## FeaturePlot
# S-RPC
FeaturePlot(pbmc, features.plot = c('Neurog2','Olig2','Atoh7','Ascl1'), 
            reduction.use = "umap", cols.use = ,no.legend = F,pt.size = 0.8)
# Photoreceptor
FeaturePlot(pbmc, features.plot = c('Otx2','Crx','Nrl','Thrb'),
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)
FeaturePlot(pbmc, features.plot = c('Rxrg','Opn1sw','Rho','Rcvrn'), 
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)
# Bipolar
FeaturePlot(pbmc, features.plot = c('Prkca','Lhx4','Prdm8','Stx1a'),
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)
# Horizontal & Amacrine
FeaturePlot(pbmc, features.plot = c('Ptf1a','Tfap2b','Tfap2a','Calb1'),
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)
FeaturePlot(pbmc, features.plot = c('Calb2','Gad2','Slc6a9'), 
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)

# Calbindin--Calb1:horizontal and amacrine cells
# Calb2—-Carentinin：amacrine subclass and ganglion cells
# Gad65--Gad2:GABAergic amacrine cells
# glycine transporter 1 (Glyt1, glycinergic amacrine cells)(Slc6a9)

# ganglion
FeaturePlot(pbmc, features.plot = c('Pou4f2','Pou4f1','Sncg','Isl1'), 
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)
# muller
FeaturePlot(pbmc, features.plot = c('Rlbp1','Gfap','S100b'), 
            reduction.use = "umap", cols.use = c("grey","red"),nCol = 3,no.legend = F)
## 放胶
FeaturePlot(pbmc, features.plot = c('Nes','Slc1a3','Fabp7'), 
            reduction.use = "umap", cols.use = c("grey","red"),no.legend = F)

#name
pbmc@meta.data$Cluster <- plyr::mapvalues(pbmc@meta.data$res.0.4, 
                                          c(2,3,5,6), c("S-RPC", "ganglion",'photoreceptor',"H & A"))

DimPlot(pbmc, reduction.use = 'umap', do.return = T,cols.use = col_hh, group.by = 'Cluster',
        plot.title = 'pc1:15_umap_E14')
# Findmarkers
pbmc <- SetAllIdent(pbmc,id = 'res.0.4')
degs <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
cluster10.markers <- FindMarkers(object = pbmc, ident.1 = , ident.2 = c(2,1,4,6,9), 
                                 min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
       
##clustprofiler
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
c3 <- subset(degs, cluster == "3")$gene
c3.fc <- subset(degs, cluster == "3")$avg_logFC

c3.idmap <- bitr(geneID = c3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
c3.id <- c3.idmap$ENTREZID[match(c3, c3.idmap$SYMBOL)]

ggo <- groupGO(gene = c3, 
               OrgDb = org.Mm.eg.db, 
               keytype = "SYMBOL", 
               ont = "CC", 
               level = 5, 
               readable = F)
ego <- enrichGO(gene = c3, 
                OrgDb = org.Mm.eg.db, 
                keyType = "SYMBOL", 
                ont = "BP", # "BP", "CC", "MF", "ALL" 
                pvalueCutoff = 0.05, 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                minGSSize = 10, 
                maxGSSize = 500, 
                readable = F, 
                pool = F)
# ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
# dropGO(x = ego, level = 5, term = NULL)
ego1 <- simplify(x = ego)
ego2 <- gofilter(x = ego, level = 4)
barplot(ego2, 
        colorBy = "p.adjust", showCategory = 10, font.size = 10, title = "clusterProfiler results") # see help
dotplot(ego, 
        colorBy = "p.adjust", showCategory = 10, font.size = 10, title = "clusterProfiler results") # see help
enrichMap(ego)
cnetplot(ego, categorySize = "pvalue", foldChange = setNames(c3.fc, c3), font.size = 1)
plotGOgraph(ego)

                                        
###TSNEPlot PC15 及PC10
pbmc_ <- SetAllIdent(pbmc_,id = "res.0.6")
pbmc_ <-  RunTSNE(pbmc_,do.fast = T,dims.use = )
TSNEPlot(pbmc_,do.lable = T,do.return = T,colors.use = ) + labs(title = "tsne_pc1:,perplexity = ")
TSNEPlot(pbmc_hGW11_,do.lable = T,do.return = T,group.by = 'type',pt.size = 0.8) + labs(title = "Sample")
TSNEPlot(pbmc_hGW11_,do.lable = F,do.return = T,group.by = 'Phase',
         colors.use = ) + labs(title = "Cellcycle")
##确定 pca
pbmc <- pbmc_

#设置离散度（perplexity）,选择较好的呈现图，设置resolution = 0.2,0.4,0.6,0.8,1.0,1.2

table(pbmc@meta.data$res)##注意颜色选择
             
for(i in seq(20,200,10)){
  pbmc <- SetAllIdent(pbmc, id = "res.0.6")
  pbmc<- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = , perplexity = i)
  png(file = paste("~/Rstudio/.../perplexity/per",i,".png",sep = ""), 
      width = 981,height = 694）
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return=T, color.use = , pt.size = ) + labs(title=i))
  dev.off()
}

pbmc <- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = , perplexity = )
             
##Resolution = 0.2,0.4, 0.6, 0.8, 1.0, 1.2
for(i in seq(0.2,0.8,0.2)){
  pbmc <- SetAllIdent(pbmc,id = paste("res.",i,sep = ""))
  pdf(file = paste("~/Rstudio/.../tsne/res.",i,".pdf",sep = ""), 
      width = 9.81,height = 6.94)
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return = T, color.use = , pt.size = ) + labs(title = i))
  dev.off()
}

for(i in seq(1,1.8,0.2)){
  pbmc <- SetAllIdent(pbmc,id = paste("res.",i,sep = "")
  pdf(file = paste("~/Rstudio/.../tsne/res.",i,".pdf",sep = ""), 
      width = 9.81,height = 6.94)
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return = T,color.use = , pt.size = ) + labs(title = i))
  dev.off()
}

#确定 resolution
pbmc<- SetAllIdent(pbmc,id = "")
TSNEPlot(object = pbmc, do.label = TRUE, do.return = T, color.use = , pt.size = ) + labs(title = "tsne_PC1:,perplexity = ")
TSNEPlot(pbmc_hGW11_,do.lable = T,do.return = T,group.by = 'type',pt.size = 0.8) + labs(title = "Sample")
TSNEPlot(pbmc_hGW11_,do.lable = F,do.return = T,group.by = 'Phase',colors.use = ) + labs(title = "Cellcycle")

##
pbmc_   <- pbmc  ###保存 pbmc_pc_sample
#
write.csv(pbmc_@meta.data, file = '/data2/hehan/Rstudio/.../日期_anno_sample名.csv')
write.csv(pbmc_@dr$tsne@cell.embeddings,file = "/data2/hehan/Rstudio/.../日期_tsne_cycle_regressout_sample名.csv")
save(pbmc_,file = "/data2/hehan/Rstudio/.../日期_seurat_cycle_regressout_sample名.RData")

## FeaturePlot
# S-RPC
FeaturePlot(pbmc_pc5, features.plot = c('Neurog2','Olig2','Atoh7','Ascl1'), 
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F,pt.size = )
FeatureHeatmap(pbmc_pc5,features.plot = c('Neurog2','Olig2'),group.by = 'type',
               pt.size = ,do.return = TRUE)
FeatureHeatmap(pbmc_pc5,features.plot = c('Atoh7','Ascl1'),group.by = 'type',
               pt.size = ,do.return = TRUE)
                      
# Photoreceptor
FeaturePlot(pbmc_pc5, features.plot = c('Otx2','Crx','Nrl','Thrb'), pt.size = ,
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F)
FeaturePlot(pbmc_pc5, features.plot = c('Rxrg','Opn1sw','Rho','Rcvrn'), pt.size = ,
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F)
FeatureHeatmap(pbmc_pc5,features.plot = c('Rcvrn'),group.by = 'type',
               pt.size = ,do.return = TRUE)
# Bipolar
#FeaturePlot(pbmc_pc5, features.plot = c('Prkca','Lhx4','Prdm8','Stx1a'), pt.size = 0.8,
#           reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F)

# Horizontal & Amacrine
FeaturePlot(pbmc_pc5, features.plot = c('Ptf1a','Tfap2b','Tfap2a','Calb1'),pt.size = 0.8,
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = T)
FeaturePlot(pbmc_pc5, features.plot = c('Calb2','Gad2','Slc6a9'), pt.size = 0.8,
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F)

FeatureHeatmap(pbmc_pc5,features.plot = c('Calb1'),group.by = 'type',
               pt.size = 0.8,do.return = TRUE)
# Calbindin--Calb1:horizontal and amacrine cells
# Calb2—-Carentinin：amacrine subclass and ganglion cells
# Gad65--Gad2:GABAergic amacrine cells
# glycine transporter 1 (Glyt1, glycinergic amacrine cells)(Slc6a9)

# ganglion
FeaturePlot(pbmc_pc5, features.plot = c('Pou4f2','Pou4f1','Sncg','Isl1'), 
            reduction.use = "tsne", cols.use = c("grey","red"),pt.size = 0.8,no.legend = T)
FeatureHeatmap(pbmc_pc5,features.plot = c('Pou4f2','Pou4f1'),group.by = 'type',
               pt.size = 0.8,do.return = TRUE)
FeatureHeatmap(pbmc_pc5,features.plot = c('Sncg','Isl1'),group.by = 'type',
               pt.size = 0.8,do.return = TRUE)

# muller
FeaturePlot(pbmc_pc5, features.plot = c('Rlbp1','Gfap','S100b'), pt.size = 0.8,
            reduction.use = "tsne", cols.use = c("grey","red"),nCol = 3,no.legend = F)
## 放胶
FeaturePlot(pbmc_pc5, features.plot = c('Ca2','Slc1a3','Fabp7'), 
            reduction.use = "tsne", cols.use = c("grey","red"),no.legend = F)
## 
FeatureHeatmap(pbmc_pc5,features.plot = c('Slc1a2','Fabp7'),group.by = 'type',
               pt.size = 0.8,do.return = TRUE)

#name
pbmc_pc5@meta.data$Cluster <- plyr::mapvalues(pbmc_pc5@meta.data$res.0.6, c(3,6,5,7,10,11,12,14,15,16), 
                                              c("Photoreceptor", "Photoreceptor", "S-RPC", 
                                                "Horizontal & Amrcrine",'ganglion-related',
                                                'ganglion-related','out-segment','ganglion',"H & A",'muller'))
pbmc_pc5 <- SetAllIdent(pbmc_pc5, "Cluster")
TSNEPlot(pbmc_pc5,do.label = TRUE,colors.use = bertie.color,pt.size = 0.8)

## 间充质干细胞相关 marker
#Sca1(Atxn1)、CD105 (Eng)、CD140b (Pdgfrb)、Gpc3、 Tagln、 Dcn 、Col1a2 、Acta2 
FeaturePlot(pbmc_pc5, features.plot = c('Atxn1',"Eng", "Pdgfrb",'Gpc3'), 
            reduction.use = "tsne", cols.use = c("grey","red"))
FeaturePlot(pbmc_pc5, features.plot = c("Tagln", "Dcn",'Col1a2','Myl9'), 
            reduction.use = "tsne", cols.use = c("grey","red"))                      

# Findmarkers
res0.6_markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
cluster10.markers <- FindMarkers(object = pbmc, ident.1 = , ident.2 = c(2,1,4,6,9), 
                                 min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
write.csv(res0.6_markers, file = ".csv")
write.csv(cluster10.markers,file = ".csv")

#vinplot
VlnPlot(object = pbmc, features.plot = c(""))
 
 #展示 cluster
 for(i in 0:~){
   pbmc@meta.data$cluster <- ifelse(pbmc@meta.data$res.0.6 == i,i,"others")
   pbmc <- SetAllIdent(pbmc, id = "cluster")
   pdf(file = paste("~/Rstudio/.../clusterplot/",i,"_cluster.pdf",sep = ""),
       width = 7.10, height = 6.60)
   print(TSNEPlot(pbmc, do.return=T)+scale_color_manual(values = c("red","grey")))
   dev.off()
 }
 
 #确定群体名字，if needed
pbmc <- SetAllIdent(pbmc, id = "")#已设置 ident 不需要
                      
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = ...)
                      
myFeaturePlot <- function(pbmc, features.plot,pt.size, nrow = NULL, ncol = NULL, dr = c("tsne","pca","ccscore","avscore","umap"), cc.args = list(th.g1s = 2, th.g2m = 2),...){
  require(ggplot2)
  require(gridExtra)
  dr <- dr[1]
  ggData <- as.data.frame(cbind(pbmc@dr[[dr]]@cell.embeddings,FetchData(pbmc, features.plot)))
  colnames(ggData) <- c(colnames(pbmc@dr[[dr]]@cell.embeddings),gsub("-",".",features.plot))
  # print(feature.tmp)
  # ggData[,feature.tmp] <- t(pbmc@data[feature.tmp,])
  if(dr == 'tsne') {
    xx <- "tSNE_1"
    yy <- "tSNE_2"
  }
  if(dr == 'pca'){
    xx <- "PC1"
    yy <- "PC2"
  }
  if(dr == 'dpt'){
    xx <- "DC1"
    yy <- "DC2"
  }
  if(dr == 'ccscore'){
    xx <- "G1S.Score"
    yy <- "G2M.Score"
  }
  if(dr == 'avscore'){
    xx <- "Artery.Score"
    yy <- "Venous.Score"
  }
  if(dr == 'umap'){
    xx <- 'UMAP1'
    yy <- 'UMAP2'
  }
  col_stage_feature <- c('grey',"#fff0bc","#f05a28")
  ggl <- lapply(features.plot, function(feature){
    p <- ggplot(ggData) + geom_point(mapping = aes_string(x = xx, y = yy, color = gsub("-",".",feature)), size = pt.size) + 
      scale_color_gradientn(colours = col_stage_feature) + 
      xlab(label = xx) + ylab(label = yy) +
      theme(legend.title = element_blank()) + ggtitle(feature) 
    if(dr == "ccscore"){
      th.g1s <- cc.args$th.g1s
      th.g2m <- cc.args$th.g2m
      ccx <- ceiling(max(pbmc@meta.data$G1S.Score))
      ccy <- ceiling(max(pbmc@meta.data$G2M.Score))
      p <- p + geom_linerange(mapping = aes_(x = th.g1s, ymin = 0, ymax = th.g1s)) + 
        geom_segment(mapping = aes_(x = 0, y = th.g2m, xend = ccx, yend = th.g2m)) + 
        geom_segment(mapping = aes_(x = th.g1s, y = th.g2m, xend = ccx, yend = ccy)) +
        geom_text(mapping = aes_(th.g1s, quote(0.3), label = quote("Quiescent"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, quote(0.3), label = quote("G1"), hjust = 1.1)) +
        geom_text(mapping = aes_(ccx, ccy-2, label = quote("S"), hjust = 1.1)) +
        geom_text(mapping = aes_(th.g1s/2, ccy-1, label = quote("G2M"), hjust = 1.1))
    }
    p
  })
  grid.arrange(grobs = ggl, nrow= nrow,ncol = ncol)
}


myFeaturePlot(pbmc_E14,features.plot = c('Rax','Vsx2','Sox2','Pax6'), dr = 'umap',nrow = 1)
myFeaturePlot(pbmc_E14,features.plot = c('Neurog2','Olig2','Atoh7'),dr = 'umap',nrow = 1)

# Photoreceptor
myFeaturePlot(pbmc_E14,features.plot = c('Otx2','Crx','Nrl','Thrb'),dr = 'umap',nrow = 1)
myFeaturePlot(pbmc_E14,features.plot = c('Rxrg','Opn1sw','Rho','Rcvrn'),dr = 'umap')

# Horizontal & Amacrine
myFeaturePlot(pbmc_E14,features.plot = c('Ptf1a','Tfap2b','Tfap2a'),dr = 'umap',nrow = 1)
myFeaturePlot(pbmc,features.plot = c('Calb2','Gad2','Slc6a9' ),dr = 'umap')

# Calbindin--Calb1:horizontal and amacrine cells
# Calb2—-Carentinin：amacrine subclass and ganglion cells
# Gad65--Gad2:GABAergic amacrine cells
# glycine transporter 1 (Glyt1, glycinergic amacrine cells)(Slc6a9)

# ganglion
#myFeaturePlot(pbmc_E14,features.plot = c('Pou4f2','Pou4f1','Rxrg'), dr = 'umap',nrow = 1)

# muller
myFeaturePlot(pbmc, features.plot = c('S100b','Gfap'), dr = 'umap',ncol = 2)

## 放胶
myFeaturePlot(pbmc_E14,features.plot = c('Nes','Slc1a3','Fabp7'),dr = 'umap',nrow = 1)              

 
 



