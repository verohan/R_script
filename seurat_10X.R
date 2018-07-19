
setwd()
options(stringsAsFactors = F)
rm（list = ls(()）

##library
library(Seurat)
library(scater)
library(RColorBrewer)
library(readxl)

##colors
pal_rainbow <- c("red","orange","yellow","green","cyan","blue","purple")
bertie.color <- c("#D0B38A", "#A3D171", "#533E88", "#7957A3", "#000000", "#E63325", "#0A8041", "#C5208E", "#3CBDED", "#3B55A3", "#D691BE", "#D23E28", "#6474B6", "#4288C8", "#80A469", "#FFCF3F", "#FBCC9F", "#F06360", "#437BB2", "#A43E40", "#206767", "#779E43", "#258950", "#F7D059", "#ED803B") 
pal_dolphin <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF", "grey")
pal_cluster <- c("#FF00AE", "#A020F1", "#000000", "#0403E5", "#FF8C01", 
                 "#8B0101", "#007502", "#FE0000", "#FFFF01", "#FF99CB", 
                 "#4A95FB", "#61FE69", "#9A7A01", "#017F8B", "#05FDFF",
                 "#D0B38A", "#533E88", "#D23E28", "#80A469", "#F06360")
pal_location <- brewer.pal(12, name = "Paired")
pal_stage <- c("indianred", "steelblue", "blue3")

#引入细胞周期 marker（g1$s,g2$m）
cell_cycle <- as.data.frame(read_excel("~/Rstudio/cell cycle/cell cycle list.xlsx",
                                       sheet = 1))

#segregate this list into markers of G1/S phase, G2/M phase and markers of S phase
g1s.genes <- as.character(na.omit(cell_cycle$G1/S))
g2m.genes <- as.character(na.omit(cell_cycle$G2/M))

#更改大小写(human不用此步)
g1s.genes <- tolower(g1s.genes)
g2m.genes <- tolower(g2m.genes)

#首字母大写(human不用此步)
library(Hmisc)
g1s.genes <- capitalize(g1s.genes)
g2m.genes <- capitalize(g2m.genes)

#读入10X/mm10
data <- Read10X(data.dir="~/Rstudio/...")
ncol(as.matrix(data))

pbmc <- CreateSeuratObject(raw.data = data,min.cells = 3,min.genes = 200)#()

plot(pbmc@meta.data$nUMI, pbmc@meta.data$nGene, 
     xlab = "nUMI - the number of transcripts", 
     ylab = "nGene - the number of genes")

#判断可能会过滤细胞类型，是否为 needed 细胞
pbmc@meta.data$filt_cells <- ifelse(pbmc@meta.data$nGene < 2500,"filt_cells","others")
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 30, do.print = TRUE, 
               pcs.print = 1:10, genes.print = 20)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2, group.by = "filt_cells")
FeaturePlot(pbmc, features.plot = c("GYPA", "HBB", "GATA1", "KLF1"), 
            reduction.use = "pca", cols.use = c("grey","red"))

# evaluate threshold for percent.mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), 
                   value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / 
  Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, 
                    col.name = "percent.mito")
VlnPlot(object = pbmc, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)

# filter cells
pbmc <- FilterCells(object = pbmc, 
                    subset.names = c("nGene",'nUMI', "percent.mito"), 
                    low.thresholds = c(1000, 4000,-Inf), 
                    high.thresholds = c(40000, 60000,0.1))
ncol(pbmc@data)

write.csv(pbmc@data, file = "~/Rstudio/... .csv")  #此为保存细胞过滤后表达矩阵
pbmc_copy <- pbmc #备份

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
#(对于 scale.factor 的选择：应选择最接近其 nUMI 中位数的数量级)，normalize 后的文件覆盖原 data 表达矩阵；
#seurat中Normolization（值为 log（2）与手动取 log，是不一样的，手动取 log 可能会在后续寻找 df markers 中其 logFC 值可能会被放大；）

#Detection of variable genes across the single cells
pbmc <- FindVariableGenes(object = pbmc , mean.function = ExpMean , dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

#Scaling the data and removing unwanted sources of variation
pbmc <- ScaleData(object=pbmc)#vars.to.regress = c("nUMI" , "percent.mito")#(vars.to.regress,对 UMI 范围不是特别大，不用执行回归)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 50, do.print = TRUE, 
               pcs.print = 1:10, genes.print = 20)  #（需查看主成分中是否有细胞周期相关基因）
#Assign Cell-Cycle Scores
pbmc <- CellCycleScoring(object = pbmc,  s.genes = g1s.genes, 
                         g2m.genes = g2m.genes, set.ident = TRUE)
head(x = pbmc@meta.data)

#Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
pbmc <- RunPCA(object = pbmc, pc.genes = c(g1s.genes, g2m.genes), do.print = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2)#(主成分展示)

#Regress out cell cycle scores during data scaling，if needed;
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("S.Score", "G2M.Score"), display.progress = FALSE)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, pcs.compute = 50, do.print = TRUE, pcs.print = 1:5, genes.print = 30)

pbmc <- RunPCA(object = pbmc, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)

VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(pbmc, dim.1 = 1, dim.2 = 2)


#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, pcs.compute = 30 ,do.print = TRUE, pcs.print = 1:10, genes.print = 20) 


# Visualize the distribution of cell cycle markers across
RidgePlot(object = pbmc, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), 
          nCol = 2)

#Determine statistically significant principal components
pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = pbmc, PCs = 1:50)

PCElbowPlot(object = pbmc, num.pc = 50)

#save PCA
pbmc <- ProjectPCA(pbmc,do.print = FALSE)
pdf(file = "~/Rstudio/.../pca/pc1:6.pdf", width = 9.81, height = 6.49)
PCHeatmap(object = pbmc, pc.use = 1:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

#Cluster the cells 
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = __, 
                     resolution = seq(0.2,2,0.2), print.output = 0, save.SNN = TRUE)
write.csv(pbmc@meta_data,file = "~/Rstudio/.../anno.csv")

#设置离散度（perplexity）,选择较好的呈现图，设置resolution = 0.2,0.4,0.6,0.8,1.0,1.2
for(i in seq(20,200,10)){
  pbmc <- SetAllIdent(pbmc, id = "res.0.4")
  pbmc<- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = __, perplexity = i)
  pdf(file = paste("~/Rstudio/.../perplexity/per",i,".pdf",sep = ""), 
      width = 9.81,height = 6.94,onefile = F)
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return=T) + labs(title=i))
  dev.off()
}

#Resolution=0.2, 0.4, 0.6, 0.8, 1.0, 1.2
for(i in seq(0.2,0.8,0.2)){
  pbmc <- SetAllIdent(pbmc,id = paste("res.",i,sep = ""))
  pbmc <- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = __, perplexity = __)
  pdf(file = paste("~/Rstudio/.../tsne/res.",i,".pdf",sep = ""), 
      width = 9.81,height = 6.94)
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return = T) + labs(title = i))
  dev.off()
}

for(i in seq(1,1.8,0.2)){
  pbmc <- SetAllIdent(pbmc,id = paste("res.",i,sep = ""))
  pbmc <- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = __, perplexity = __)
  pdf(file = paste("~/Rstudio/.../tsne/res.",i,".pdf",sep = ""), 
      width = 9.81,height = 6.94)
  print(TSNEPlot(object = pbmc, do.label = TRUE, do.return = T) + labs(title = i))
  dev.off()
}

#确定 resolution
pbmc<- SetAllIdent(pbmc,id = "__")
pbmc <- RunTSNE(object = pbmc, do.fast = TRUE, dims.use = __, perplexity = __)
TSNEPlot(object = pbmc, do.label = TRUE, do.return = T, color.use = __) + labs(title = "__")
DimPlot(object, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
  cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE,
  cols.use = NULL, group.by = "ident", pt.shape = NULL,
  do.hover = FALSE, data.hover = "ident", do.identify = FALSE,
  do.label = FALSE, label.size = 4, no.legend = FALSE, no.axes = FALSE,
  dark.theme = FALSE, plot.order = NULL, cells.highlight = NULL,
  plot.title = NULL, vector.friendly = FALSE, png.file = NULL,
  png.arguments = c(10, 10, 100), ...)
  
# cell cycle_tsne
#pbmc <- AddMetaData(pbmc, metadata = pbmc@meta.data)#(???????有疑问)(待完善)
pbmc <- SetAllIdent(pbmc, id = "Phase")
pbmc <- RunTSNE(object = pbmc, seed.use = 1, dims.use = __, do.fast = TRUE, perplexity = __) 
pdf(file = "~/Rstudio/.../cellcycle_tsne.pdf", width = 9.81, height = 6.40)
TSNEPlot(object = pbmc, do.label = TRUE)
dev.off()


#
write.csv(pbmc@dr$tsne@cell.embeddings,file = "")
save(pbmc,file = "~")

# Findmarkers
res0.6_markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
cluster10.markers <- FindMarkers(object = pbmc, ident.1 = , ident.2 = c(2,1,4,6,9), 
                                 min.pct = 0.25, logfc.threshold = 1, only.pos = TRUE)
write.csv(res0.6_markers, file = ".csv")
write.csv(cluster10.markers,file = ".csv")

#vinplot
VlnPlot(object = pbmc, features.plot = c(""))
 
#featureplot
FeaturePlot(pbmc, features.plot = c(''),  pt.size = 1, cols.use = c("grey", "red"), 
             reduction.use = "tsne", nCol = NULL, do.return = FALSE)
FeatureHeatmap(pbmc,features.plot = c(""),group.by = "", pt.size = 1, key.position = "top")#cca展示不同批次需要 
 
 #展示 cluster
 for(i in 0:~){
   pbmc@meta.data$cluster <- ifelse(pbmc@meta.data$res.~ == i,i,"others")
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

 
 



