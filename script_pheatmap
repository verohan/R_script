## 准备文件：data,data_expMat,anno(colnames=cluster),dfmarkers,gene,col(颜色)
data <- as.matrix(rawdata)
data_expMat <- apply(data,2,function(x) log2(1e4*x/sum(x) +1))
anno <- data$anno
colnames(anno) <- "cluster"

### 附值颜色
col <- setNames(c("red3", "green3", "blue3", "purple3", "pink", "navy","grey"),unique(anno$cluster))

#对差异基因进行排序（Cluster1-。。。，logfc降序）
library(dplyr)
order_marker_by_fc<-function(marker,order){
  cluster<-order
  gene<-c()
  deg<-list()
  for(i in 1:length(cluster)){
    tmp<-marker[marker$cluster==cluster[i],]
    tmp<-tmp[order(tmp$avg_logFC,decreasing = T),]
    top12<-as.character(na.omit(tmp$gene[1:12]))
    gene<-c(gene,top12)
    deg[[i]]<-tmp
  }
  all_deg<-do.call(rbind,deg)
  result<-list(deg=all_deg,gene=gene)
  return(result)
}
dfmarkers <- order_marker_by_fc(marker = dfmarkers,order = c('0','1','3','5','2','6','4'))##(order:anno中的 cluster)
gene1 <- dfmarkers$gene

#heatmap
do_heatmap<-function(data_exprMat,scale=F,gene1,anno,order=NULL,col){
  library(pheatmap)
  library(cowplot)
  library(ggplotify)
  if(is.null(order)==TRUE){
    if(scale==TRUE){
      data<-t(scale(t(data_exprMat)))
      data[data>2]<-2
      data[data< -2]<- -2
      pheatmap(data[gene1,rownames(anno)[order(anno$cluster)]],
               cluster_cols = F,annotation_col = anno,show_colnames = F,
               cluster_rows = F,legend = T,color = colorRampPalette(c("purple","black","yellow"))(100),border_color = NA)
      
    }else{
      data<-data_exprMat
      pheatmap(data[gene1,rownames(anno)[order(anno$cluster)]],
               cluster_cols = F,annotation_col = anno,show_colnames = F,
               cluster_rows = F,legend = T,color = colorRampPalette(c("navy","white","red"))(100),border_color = NA)}
  }else{
    cell<-c()
    for(i in order){
      cell1<-rownames(anno)[anno$cluster==i]
      cell<-c(cell,cell1)}
    if(scale==TRUE){
      data<-t(scale(t(data_exprMat)))
      data[data>2]<-2
      data[data< -2]<- -2
      pheatmap(data[gene1,cell],cluster_cols = F,annotation_col = anno,show_colnames = F,cluster_rows = F,annotation_colors = col,legend = T,color = colorRampPalette(c("purple","black","yellow"))(100),border_color = NA)
    }else{
      data<-data_exprMat
      pheatmap(data[gene1,cell],cluster_cols = F,annotation_col = anno,show_colnames = F,cluster_rows = F,annotation_colors = col,legend = T,color = colorRampPalette(c("navy","white","red"))(100),border_color = NA)
    }
  }
}

do_heatmap(data_expMat,gene1 = gene, anno = anno,order=c("0","1","3","5","2","6","4"),col=list(cluster=col))
