setwd()
options(stringsAsFactors = F)
rm(list = ls(())
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
pal_cluster <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", 
                 "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", 
                 "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF",
                 "#D0B38A", "#533E88", "#D23E28", "#80A469", "#F06360")##20
pal_location <- brewer.pal(12, name = "Paired")
pal_stage <- c("indianred", "steelblue", "blue3")
col_hh <- c('#9DC8C8','#D1B6E1','#519D9E','#E53A40',
            '#30A9DE','#F17F42','#566270','#8CD790',
            '#2EC4B6', '#44633F','#D81159','#4F86C6',
            '#AACD6E','#F16B6F','#E3E36A','#D8E6E7',
            '#D09E88','#7200da','#f100e5','#005f6b')##20

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

##
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
doublets_2 <- colnames(pbmc@raw.data)
[apply(pbmc@raw.data[c('Kdm5d','Eif2s3y','Gm29650','Uty','Ddx3y'),],2,function(x) any(x>0)) & pbmc@raw.data["Xist",]>0]
pbmc@meta.data$dblets_2 <- ifelse(rownames(pbmc@meta.data) %in% doublets_2,'doublet','singlet')
doublets <- union(doublets_1,doublets_2)
                                            
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

#unfilt cells_PCA 查看细胞判断其低质量 or 独立群体细胞             
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- ScaleData(pbmc)
#(对于 scale.factor 的选择：应选择最接近其 nUMI 中位数的数量级)，normalize 后的文件覆盖原 data 表达矩阵；
#seurat中Normolization（值为 log（2）与手动取 log，是不一样的，手动取 log 可能会在后续寻找 df markers 中其 logFC 值可能会被放大；）
#vars.to.regress = c("nUMI" , "percent.mito")#(vars.to.regress,对 UMI 范围不是特别大，不用执行回归) 
             
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 30, do.print = TRUE, 
               pcs.print = 1:10, genes.print = 20)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "filt_cells_nGene")
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "filt_cells_nUMI")
FeaturePlot(pbmc, features.plot = c("GYPA", "HBB", "GATA1", "KLF1"), 
            reduction.use = "pca", cols.use = c("grey","red"))###红细胞(GYPA:CD235a)

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

## 仅根据细胞 nUMI，nGene，mito过滤
pbmc <- FilterCells(object = pbmc, 
                    subset.names = c("nGene",'nUMI', "percent.mito"), 
                    low.thresholds = c(1000, 4000,-Inf), 
                    high.thresholds = c(40000, 60000,0.1))
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
pbmc <- ScaleData(pbmc,vars.to.regress = c('nUMI','S.Score','G2M.Score'))
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
hvg <- rownames(pbmc@hvg.info[order(pbmc@hvg.info[,3],decreasing = TRUE),])[1:2000]

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
pbmc_pc15 <- pbmc
pbmc_pc10 <- pbmc
                                        
pbmc_ <- FindClusters(object = pbmc_, reduction.type = "pca", dims.use = , 
                     resolution = seq(0.2,2,0.2), print.output = 0, save.SNN = TRUE)
                                        
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

 
 



