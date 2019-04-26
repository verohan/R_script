##----heatmap----
degs_use <- read.csv(file = 'Rstudio/mRPC/merge/data_190308*/degs/degs_filtraw.csv')

DoHeatmap(pbmc_merge,use.scaled = F, genes.use = degs_use$gene,
          group.by = 'Cluster_raw',disp.min = 0,disp.max = 2.5,
          group.order = c('RPC','Restricted-RPC','Precursors','PR','RGC','AC/HC'),
          col.low = 'darkblue',col.mid = 'white',col.high = 'red',
          group.cex = 10,cex.col = 0.1,group.label.loc = 'top',
          slim.col.label = T,remove.key = F,do.plot = T)
range(as.matrix(pbmc_merge@scale.data))[degs_use$gene,]
quantile(pbmc_merge@data)

## from hj
degs_raw <- read.csv('Rstudio/mRPC/merge/data_190424/degs/degs_raw1.csv',header = T,row.names = 1)

pheatmap_cluster <- function(pbmc, feature_genes, cluster_mean,  feature_fontsize = 12){
  data_mean <- data.frame(row.names = feature_genes)
  pbmc@meta.data$cluster_heatmap <- pbmc@meta.data[,cluster_mean]
  for (i in unique(pbmc@meta.data$cluster_heatmap)) {
    data_tmp <- pbmc@data[feature_genes,rownames(subset(pbmc@meta.data, cluster_heatmap == i))]
    mean_tmp <- apply(data_tmp, 1, mean)
    data_mean <- data.frame(data_mean, mean_tmp)
  }
  colnames(data_mean) <- unique(pbmc@meta.data$cluster_heatmap)
  library(pheatmap)
  pheatmap(data_mean, cluster_cols = F, cluster_rows = F, fontsize_row = feature_fontsize, legend = F,
           color = colorRampPalette(colors =c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"))(100))
}
##
## heatmap
degs_use <- read.csv('Rstudio/mRPC/ACHC/acah.csv',header = T,row.names = 1)
degs_matrix <- as.matrix(pbmc_ACHC@data)
degs_matrix <- degs_matrix[as.character(degs_use$gene),]# #degs 里为factor

range(degs_matrix)#0,3.17

phmat <- t(scale(t(degs_matrix)))
phmat[phmat>2]<- 2
phmat[phmat < -2]<- -2



dim(phmat)
dim(pbmc_ACHC@meta.data)

table(pbmc_ACHC@meta.data$res.0.4)
pheatmap::pheatmap(phmat[,order(pbmc_ACHC@meta.data$res.0.4)],show_colnames = F,kmeans_k = NA,
                   annotation_col = pbmc_ACHC@meta.data[,'res.0.4', drop =F],
                   color = colorRampPalette(c('midnightblue','dodgerblue3','white','goldenrod1','darkorange2'))(100),
                   cluster_rows = F,cluster_cols = F,border_color = NA,scale = 'none')
