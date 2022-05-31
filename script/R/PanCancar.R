### ---------------
###
### Create: Yuan.Sh
### Date: 2022-05-31 16:15:25 
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
###
### ---------------

############# Step.00 ############# 
# clean
rm(list = ls())
gc()
# packages
library(Seurat)
library(stringr)
library(scales)
library(magrittr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(eoffice)
library(dplyr)
library(reshape2)
library(ggpubr)
theme_set(theme_cowplot())
# show_col(use_colors)
# options
options(stringsAsFactors = F)
options(as.is = T)
setwd('/media/yuansh/14THHD/scMutation/')
gene.list = read.table('geneList.txt')[,1]
############# GSE116237 SKCM ############# 
# prepare data
if(F){
  # get raw-count matrix 
  expr = read.csv('/media/yuansh/14THHD/BscModel-V4/GSE116237/GSE116237.csv', row.names=1)
  cell.id = colnames(expr)
  # get pathway detection result
  for(gene in gene.list){
    file.name = paste0('Pan_cancer_Pathyway_detection/GSE116237/Pan_cancer_',gene,'/','ScMutation_',gene,'_input_sce_embeddings.csv')
    assign(paste0('GSE116273_',gene), read.csv(file.name, row.names = 1)[cell.id,])
  }
  # get meta data
  if(T){
    df = as.data.frame(cell.id)
    df$TimePoint = 'T56'
    df[grep('T0',df$cell.id),]$TimePoint = 'T0'
    df[grep('T4',df$cell.id),]$TimePoint = 'T4'
    df[grep('T28',df$cell.id),]$TimePoint = 'T28'
    rownames(df) = df$cell.id
    for(gene in gene.list){
      ids = get(paste0('GSE116273_',gene))
      df[[gene]] = ifelse(ids$predict == 'WT',0,1)
    }
    df$mut_score = apply(df[,gene.list],1,sum)
    
    summary(df$mut_score)
    
    df$mut_type = ifelse(df$mut_score<=70,'G1',
                         ifelse(df$mut_score<=85, 'G2',
                                ifelse(df$mut_score<=94, 'G3','G4')))
    table(df$mut_type)
  }
  sce = CreateSeuratObject(count = expr,
                           meta.data = df)
  saveRDS(sce,'script/R/GSE116237/GSE116237-PathWay-Detection-seur-obj.rds')
}
# load data
if(T){
  GSE116237 = readRDS('script/R/GSE116237/GSE116237-PathWay-Detection-seur-obj.rds')
  GSE116237$TimePoint = factor(GSE116237$TimePoint, levels = c('T0','T4','T28','T56'))
  plot.path = 'script/R/GSE116237/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
  timepoint.cols = c(T0='#F8766D',T4='#ECA86E',T28='#FFEF9C',T56='#7FA4D1')
}
df = GSE116237@meta.data
df = df[,c('TimePoint','mut_type','mut_score',gene.list)]
write.csv(df,'~/Desktop/GSE116237_SKCM.csv')
ggplot(GSE116237@meta.data, aes(TimePoint,mut_score,fill =TimePoint)) +geom_boxplot()
ggplot(GSE116237@meta.data, aes(mut_type,nFeature_RNA,fill =mut_type)) +geom_boxplot()
ggplot(GSE116237@meta.data, aes(TimePoint,nFeature_RNA,fill =TimePoint)) +geom_boxplot()
ggplot(GSE116237@meta.data, aes(TimePoint,nFeature_RNA,fill =TimePoint)) +geom_boxplot()

# Dual
if(F){
  sce = read.csv('script/R/GSE116237/GSE116237_Dual_cell_plot.csv',row.names = 1)
  sce = cbind(sce,GSE116237@meta.data)
  df = sce[,c('immune_score','epi_score','driver_gene_score','mut_score')]
  colnames(df) = c('immscore','episcore','driverscore','burd')
  
  T0 = df %>% rownames() %>% grep('T0',.,value = T)
  T4 = df %>% rownames() %>% grep('T4',.,value = T)
  T28 = df %>% rownames() %>% grep('T28',.,value = T)
  df$group = ifelse(rownames(df) %in% T0,'T0',ifelse(rownames(df) %in% T4,'T4',ifelse(rownames(df) %in% T28,'T28','T56')))
  
  df$group = factor(df$group,levels = c('T0','T4','T28','T56'))
  x= ggplot(df,aes(immscore,driverscore))+
    stat_density2d(geom ="polygon",aes(fill = ..level..),bins=30 )+
    scale_fill_gradientn(colours=c("#1B7EAB", "#ffffff" ,"#EF7B66"))+#, trans="log"
    theme_bw()+labs(title='')+
    geom_hline(aes(yintercept=0),color = '#990000',linetype='dashed')+
    geom_vline(aes(xintercept=0),color = '#990000',linetype='dashed')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+
    facet_wrap(~group)
  x
  topptx(x,file="script/R/GSE116237/GSE116237_Dual-four_grid-plot.pptx")
  
  library(ggpubr) 
  my_comparisons <- list( c("T0", "T4"), c("T0", "T28"), c("T0", "T56"))
  df$group = factor(df$group,levels = c('T0','T4','T28','T56'))
  x = ggplot(df,aes(x=group,immscore,fill= group))+
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons,method = 't.test')+ # Add pairwise comparisons p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = ".all.", hide.ns = TRUE,label.y = 6) +
    labs(title='')+    theme(legend.position = "top") +scale_fill_manual(values = timepoint.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank())
  x
  topptx(x,file="script/R/GSE116237/GSE116237_Dual-immune_score_box-plot.pptx")
}
# plot mut burd
if(F){
  df = GSE116237@meta.data
  my_comparisons <- list( c("T0", "T4"), c("T0", "T28"), c("T0", "T56"))
  x=ggplot(df,aes(x = TimePoint, y = mut_score,fill=TimePoint)) + 
    geom_boxplot() + 
    scale_fill_manual(values = c("#DC0000B2","#7E6148B2",
                                 "#E64B35B2",'#698EC3'))+
    stat_compare_means(comparisons = my_comparisons,method = 't.test')+ # Add pairwise comparisons p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = ".all.", hide.ns = TRUE,label.y = 130) +
    labs(title='')+
    scale_fill_manual(values = timepoint.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank())  
  x
  topptx(x,file=paste0(plot.path,"GSE116237-Mut_score_T0-T56.pptx"))
}
# plot pathway mut freq
if(F){
  df = GSE116237@meta.data
  planes <- group_by(df, TimePoint)
  ids = 'summarise(planes,'
  for(gene in gene.list){
    if(gene == gene.list[100]){
      ids = paste0(ids,gene,'=','mean(',gene,'))')
    }else(ids = paste0(ids,gene,'=','mean(',gene,'),'))
  }
  df.gene.freq=eval(parse(text = ids)) %>% as.data.frame()
  rownames(df.gene.freq) = df.gene.freq[,1]
  
  df.gene.freq = df.gene.freq[,-1] %>% t()
  std=apply(df.gene.freq,1,var)
  std =std[ order(std,decreasing = T)] %>% names()
  pheatmap::pheatmap(df.gene.freq[std[1:20],],scale = 'row',
                     cellwidth = 10,cellheight =  10,cluster_cols = F,
                     filename = paste0(plot.path,"GSE116237-Mut_Freq_T0-T56.pdf"))
  
  ids = abs(df.gene.freq[,4] - df.gene.freq[,2])
  ids = ids[order(ids,decreasing = T)] %>% names()
  #ids = std[1:20]
  df = df.gene.freq[ids[1:15],] %>% t %>% as.data.frame()
  df$group = factor(rownames(df), levels = c('T0','T4','T28','T56'))
  mydata<-melt(
    df,                       
    id.vars=c("group"),  
    variable.name="Gene",         
    value.name="Freq"           
  )
  x = ggplot(mydata,aes(x = group, y = Freq,fill=group)) + facet_wrap(~Gene)+
    geom_bar(stat = 'identity')+ scale_y_continuous(limits = c(0,1))+scale_fill_manual(values = timepoint.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank())  
  x
  topptx(x,file=paste0(plot.path,"GSE116237-mut-freq-bar.pptx"))
}
#gene diversity
if(F){
  sce = GSE116237
  df = sce@meta.data
  df$group = df$mut_type
  df$group = factor(df$group,levels = c('G4','G3','G2','G1'))
  x= ggplot(df, aes(nFeature_RNA,group, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())+facet_wrap(~TimePoint)
  x
  topptx(x,file=paste0(plot.path,"GSE116237_gene_diversity.pptx")) 
}

# run tsne
if(F){
  # DEGs
  sce = GSE116237
  # 1.log
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize",  scale.factor = 1e6)
  sce = FindVariableFeatures(object = sce,selection.method = "vst", nfeatures = 2000)
  sce = ScaleData(object = sce)
  # 4. PCA
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:15)
  sce = FindClusters(sce) 
  # 7.tsne
  sce=RunUMAP(sce,dims = 1:15)
  DimPlot(sce,group.by = 'mut_type')
}
# plot G1-G4 Freq
if(F){
  library(reshape2)
  table(sce$mut_type)
  table(sce$mut_type, sce$TimePoint)
  
  df = as.data.frame.array(table(sce$mut_type, sce$TimePoint))
  
  df.freq = df / (matrix( rep(table(sce$TimePoint),4), nrow=4) %>% t)
  df.freq$mut_type = rownames(df.freq)
  df.freq<-melt(
    df.freq,                    
    id.vars=c('mut_type'),  
    variable.name="TimePoint",        
    value.name="Freq"
  )
  x = ggplot(df.freq,aes(TimePoint, Freq, fill=mut_type)) + 
    geom_bar(stat = 'identity') +scale_fill_manual(values = mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank())  
  x
  topptx(x,file=paste0(plot.path,"GSE116237-g1-g4-in-t0-t56.pptx"),height = 3,width = 7.5)
}
# plot gsea1 T56_G4 vs T0_G4
if(F){
  sce = GSE116237
  T0_G4 =names( which(sce$mut_type == 'G4' & sce$TimePoint == 'T0') )
  T56_G4 =names( which(sce$mut_type == 'G4' & sce$TimePoint == 'T56') )
  length(T56_G4)
  length(T0_G4)
  degs = FindMarkers(sce, ident.1=T56_G4, ident.2=T0_G4)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0_G4-T56_G4.csv"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  
  gsea_result = egmt@result
  grep('EPITHELIAL', rownames(gsea_result),value = T)
  library(enrichplot)
  x = gseaplot2(egmt,color='#00A087',
                geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                pvalue_table = TRUE)
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T56_G4-T0_G4-ETM.pptx"))
  
  x = gseaplot2(egmt,color='#c10534',
                geneSetID = 'HALLMARK_ALLOGRAFT_REJECTION',
                pvalue_table = TRUE)
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-ALLOGRAFT_REJECTION.pptx"))

  df = gsea_result#[which(gsea_result$pvalue < 0.05),]
  
  df$Description = str_replace(df$Description,'HALLMARK_','')
  
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-GSEA-bar.pptx"), width = 8,height = 5)
}
# plot gsea1 T0_G4 vs T0_G1
if(F){
  T0_G4 =names( which(sce$mut_type == 'G4' & sce$TimePoint == 'T0') )
  T0_G1 =names( which(sce$mut_type == 'G1' & sce$TimePoint == 'T0') )
  length(T0_G1)
  length(T0_G4)
  degs = FindMarkers(sce, ident.1=T0_G4, ident.2=T0_G1)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0_G4-T0_G1.csv"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  
  gsea_result = egmt@result
  grep('EPITHELIAL', rownames(gsea_result),value = T)
  library(enrichplot)
  df = gsea_result[which(gsea_result$pvalue < 0.05),]
  df$Description = str_replace(df$Description,'HALLMARK_','')
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T0_G1-GSEA-bar.pptx"), width = 8,height = 5)
}
# plot gsea1 T56_G4 vs T0_G1
if(F){
  sce = GSE116237
  T0_G1 =names( which(sce$mut_type == 'G1' & sce$TimePoint == 'T0') )
  T56_G4 =names( which(sce$mut_type == 'G4' & sce$TimePoint == 'T56') )
  length(T56_G4)
  length(T0_G1)
  degs = FindMarkers(sce, ident.1=T56_G4, ident.2=T0_G1)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0_G1-T56_G4.csv"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = clusterProfiler::read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)

  
  
  gsea_result = egmt@result
  
  df = gsea_result[which(gsea_result$pvalue < 0.05),]
  
  df$Description = str_replace(df$Description,'HALLMARK_','')
  
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-GSEA-bar.pptx"), width = 8,height = 5)
}
# plot gsea1 T28_noG4 vs T0_G1
if(F){
  sce = GSE116237
  T0 =rownames( sce@meta.data[which(sce$TimePoint == 'T0'),])
  T4 = rownames( sce@meta.data[which(sce$TimePoint == 'T4'),])
  length(T4)
  length(T0)
  degs = FindMarkers(sce, ident.1=T4, ident.2=T0)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0-T4"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = clusterProfiler::read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  
  
  gsea_result = egmt@result
  
  df = gsea_result
  
  df$Description = str_replace(df$Description,'HALLMARK_','')
  
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-GSEA-bar.pptx"), width = 8,height = 5)
}
if(F){
  sce = GSE116237
  T0 =rownames( sce@meta.data[which(sce$TimePoint == 'T0'),])
  T28 = rownames( sce@meta.data[which(sce$TimePoint == 'T28'),])
  length(T28)
  length(T0)
  degs = FindMarkers(sce, ident.1=T28, ident.2=T0)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0-T28"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = clusterProfiler::read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  
  
  gsea_result = egmt@result
  
  df = gsea_result
  
  df$Description = str_replace(df$Description,'HALLMARK_','')
  
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-GSEA-bar.pptx"), width = 8,height = 5)
}
if(F){
  sce = GSE116237
  T0 =rownames( sce@meta.data[which(sce$TimePoint == 'T0'),])
  T56 = rownames( sce@meta.data[which(sce$TimePoint == 'T56'),])
  length(T56)
  length(T0)
  degs = FindMarkers(sce, ident.1=T56, ident.2=T0)
  write.csv(degs,  paste0(plot.path,"GSE116237-T0-T56"))
  deg = degs
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  geneset = clusterProfiler::read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  
  
  gsea_result = egmt@result
  
  df = gsea_result[which(gsea_result$pvalue < 0.05),]
  
  df$Description = str_replace(df$Description,'HALLMARK_','')
  
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T56_G4-GSEA-bar.pptx"), width = 8,height = 5)
}
# plot monocle
if(F){
  library(monocle)
  Mono_tj = sce
  Mono_matrix<-as(as.matrix(GetAssayData(Mono_tj,slot = "counts")), 'sparseMatrix')
  feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
  rownames(feature_ann)<-rownames(Mono_matrix)
  Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
  sample_ann<-Mono_tj@meta.data
  Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
  Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())
  Mono.cds <- estimateSizeFactors(Mono.cds)
  Mono.cds <- estimateDispersions(Mono.cds)
  Mono.cds <- reduceDimension(
    Mono.cds,
    max_components = 2,
    method = 'DDRTree')
  Mono.cds <- orderCells(Mono.cds,num_paths=16)
  df = read.csv('script/R/GSE116237/GSE116237_Dual_cell_plot.csv',row.names = 1)
  identical(rownames(df),colnames(Mono.cds))
  Mono.cds$dual = ifelse(df$immune_score>0 & df$driver_gene_score>0,'Dual','Non_Dual')
  plot_cell_trajectory(Mono.cds,cell_size = 1,color_by='dual')+ 
    theme(legend.position = "right")+scale_color_manual(values = c(Non_Dual='#d9e6eb',Dual='#c10534'))
}
# tsne
if(F){
  sce = GSE116237
  sce = CreateSeuratObject(counts = sce@assays$RNA@counts, 
                           meta.data = sce@meta.data,
                           min.cells = 150)
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  sce$MITF_Expr = sce@assays$RNA@data['MITF',]
  sce$AXL_Expr = sce@assays$RNA@data['AXL',]
  df = sce@meta.data
  if(F){
    x1=DotPlot(sce[,which(sce$TimePoint=='T0')],group.by = 'mut_type',feature = c('MITF_Expr','AXL_Expr'),col.min = 0)+scale_size(limits=c(0,100))+
      scale_color_gradient(  low = "#FFFFFF",
                             high = "blue",limits=c(0,3))+ theme(legend.title=element_blank(),
                                                                   axis.title.x  = element_blank(),
                                                                   axis.title.y  = element_blank())  
    x2=DotPlot(sce[,which(sce$TimePoint=='T4')],group.by = 'mut_type',feature = c('MITF_Expr','AXL_Expr'),col.min = 0)+scale_size(limits=c(0,100))+
      scale_color_gradient(  low = "#FFFFFF",
                             high = "blue",limits=c(0,3))+ theme(legend.title=element_blank(),
                                                                   axis.title.x  = element_blank(),
                                                                   axis.title.y  = element_blank())  
    x3=DotPlot(sce[,which(sce$TimePoint=='T28')],group.by = 'mut_type',feature = c('MITF_Expr','AXL_Expr'),col.min = 0)+scale_size(limits=c(0,100))+
      scale_color_gradient(  low = "#FFFFFF",
                             high = "blue",limits=c(0,3))+ theme(legend.title=element_blank(),
                                                                   axis.title.x  = element_blank(),
                                                                   axis.title.y  = element_blank())  
    x4=DotPlot(sce[,which(sce$TimePoint=='T56')],group.by = 'mut_type',feature = c('MITF_Expr','AXL_Expr'),col.min = 0)+scale_size(limits=c(0,100))+
      scale_color_gradient(  low = "#FFFFFF",
                             high = "blue",limits=c(0,3))+ theme(legend.title=element_blank(),
                                                                   axis.title.x  = element_blank(),
                                                                   axis.title.y  = element_blank())  
    x = CombinePlots(plots = list(x1,x2,x3,x4),ncol = 1)
    x
    topptx(x,file=paste0(plot.path,"GSE116237-MITF.pptx"),width = 5,height = 5)
  }
  # genes cor to mut_score
  if(F){
    ids = apply(sce@assays$RNA@counts,1,
                function(x){ids=cor.test(x,sce$mut_score,method='spearman')
                ids$estimate} )
    df = data.frame(cor=ids)
    df$Gene = rownames(df)
    df = df[order(df$cor,decreasing = T),]
    rp = grep("^RP[SL][[:digit:]]", df$Gene)
    hsp = grep("^HSP", df$Gene)
    eef = grep("^EEF", df$Gene)
    df = df[-c(rp,hsp,eef),]
    
    
    features = df$Gene[c(1:4,6:8,10:12)]
    
    sce$groups = paste(sce$TimePoint,sce$mut_type, sep = '-')
    sce$groups = factor(sce$groups, levels = c("T0-G1","T0-G2","T0-G3","T0-G4","T4-G1","T4-G2","T4-G3","T4-G4","T28-G1","T28-G2","T28-G3","T28-G4","T56-G1","T56-G2","T56-G3","T56-G4"))
    
    x = DoHeatmap(sce, features = features,group.by = 'groups',slot='data',
                  disp.min=-2.5, disp.max=2.5)+
      scale_fill_gradientn(colors = c("#7FA4D1", "white", "red"), limits=c(-2.5, 2.5))
    topptx(x,filename = paste0(plot.path,"GSE116237-proliferative-heatmap.pptx"),width = 25,height = 12)
    
    x = DotPlot(sce, features = c("NUP35","CENPK","FANCD2","KNTC1",
                                  "MCM2","ATAD2","CENPU","SMC2","HMGB2"),
                group.by = 'groups',
                dot.scale=3)+
      theme(axis.title =  element_blank(),
            axis.text.x =   element_text(size=6,angle=90),
            axis.text.y =   element_text(size=6))
    x
    topptx(x,filename = paste0(plot.path,"GSE116237-proliferative.pptx"),width = 4,height = 4)
    
  } 
  if(F){
    s.genes <- cc.genes$s.genes 
    g2m.genes <- cc.genes$g2m.genes 
    sce <- CellCycleScoring(sce, s.features = s.genes, 
                            g2m.features = g2m.genes, 
                            set.ident = TRUE) 
    sce <- RunPCA(sce, features = c(s.genes, g2m.genes)) 
    sce$groups = paste(sce$TimePoint,sce$mut_type, sep = '-')
    sce$groups = factor(sce$groups, levels = c("T0-G1","T0-G2","T0-G3","T0-G4","T4-G1","T4-G2","T4-G3","T4-G4","T28-G1","T28-G2","T28-G3","T28-G4","T56-G1","T56-G2","T56-G3","T56-G4"))
    
    df = sce@meta.data
    
    
    planes <- group_by(df, groups) %>% summarise(.,
                                                 G2M.Score = median(G2M.Score),
                                                 mut_score = median(mut_score),
                                                 S.Score = median(S.Score),
    ) 
    x = ggplot(planes,aes(mut_score,G2M.Score))+geom_point()+
      geom_smooth(method = "lm")+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),legend.position = 'none')+
      stat_cor(method = "spearman",  label.x = 70,
               label.y = 0)
    x
    topptx(x,filename = paste0(plot.path,"GSE116237-mutation_score_cor_G2.pptx"),width = 4,height = 4)
  } # g2 score to mutscore
  }
# immune turnover
if(F){
  sce = GSE116237
  df = read.csv('script/R/GSE116237/GSE116237_Dual_cell_plot.csv',row.names = 1)
  df = df[which(df$immune_score>0&df$driver_gene_score>0),] %>% rownames()
  cells = sce@meta.data[which(sce$TimePoint %in% c('T0','T56') & colnames(sce)%in%df),] %>% rownames()
  sce = subset(sce, cells = cells)
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  Idents(sce) = sce$TimePoint
  degs = FindAllMarkers(sce)
  genes = read.csv('script/R/GeneList.csv')
  genes = genes$Symbol
  ids = intersect(rownames(degs),genes)
  deg = degs[ids,]
  deg = deg[which(deg$p_val_adj<0.05),]
  DoHeatmap(sce,features = rownames(deg))
  DotPlot(sce,features = rownames(deg))+ theme(axis.text.x = element_text(angle = 90))
}
# immune diveristy
if(F){
  genes = read.csv('script/R/GeneList.csv')
  genes = genes$Symbol
  sce = GSE116237[genes,]
  sce = CreateSeuratObject(counts = sce@assays$RNA@counts)
  df = sce@meta.data
  identical(rownames(df),rownames(GSE116237@meta.data))
  df$TimePoint = GSE116237$TimePoint
  df$group = GSE116237$mut_type
  df$group = factor(df$group,levels = c('G1','G4'))
  df = df[which(df$group %in% c('G1','G4')),]
  
  x= ggplot(df, aes(group,nFeature_RNA, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())+ 
    stat_compare_means()

  x
  topptx(x,file=paste0(plot.path,"GSE116237_immune_gene_diversity.pptx")) 
}

############# GSE132465 COAD ############# 
# prepare data
if(F){
  # read data
  expr = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt',sep = '\t', row.names = 1)
  colnames(expr) = str_replace(colnames(expr),'\\.','-')
  meta = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt',sep = '\t', row.names = 1)
  cancer = read.csv('/media/yuansh/14THHD/BscModel-V4/GSE132465/merge_data_predict_Integration_.csv',row.names = 1)
  rownames(cancer) = str_replace(rownames(cancer) ,'\\.','-')
  cell.id = read.csv('Pan_cancer_Pathyway_detection/GSE132465/Pan_cancer_CNTNAP5/ScMutation_CNTNAP5_input_sce_embeddings.csv', row.names = 1) %>% rownames()
  cell.id = intersect(cell.id, colnames(expr))
  # get pathway detection result
  for(gene in gene.list){
    file.name = paste0('Pan_cancer_Pathyway_detection/GSE132465/Pan_cancer_',gene,'/','ScMutation_',gene,'_input_sce_embeddings.csv')
    assign(paste0('GSE132465_',gene), read.csv(file.name, row.names = 1)[cell.id,])
  }
  # get meta data
  if(T){
    df = as.data.frame(cell.id)
    for(gene in gene.list){
      ids = get(paste0('GSE132465_',gene))
      df[[gene]] = ifelse(ids$predict == 'WT',0,1)
    }
    df$mut_score = apply(df[,gene.list],1,sum)
    
    summary(df$mut_score)
    
    df$mut_type = ifelse(df$mut_score<=5,'G1',
                         ifelse(df$mut_score<=20, 'G2',
                                ifelse(df$mut_score<=35, 'G3','G4')))
    table(df$mut_type)
  }
  col.names = c(gene.list,'mut_score','mut_type')
  ids = intersect(cell.id, rownames(meta))
  for(i in col.names){
    meta[[i]] = NA
    meta[ids,i] =df[,i]
  }
  
  meta[['casee']] = NA
  cancer = cancer[intersect(cell.id, rownames(cancer)),]
  meta[intersect(cell.id, rownames(cancer)),'casee'] = cancer$scale_predict
  
  dual = read.csv('script/R/GSE132465/GSE132465_Dual_cell_plot.csv',row.names = 1)
  dual = dual[ids,]
  meta$im_score = NA
  meta$dr_score = NA
  meta[ids,'im_score'] = dual$immune_score
  meta[ids,'dr_score'] = dual$driver_gene_score
  
  meta$dual = NA
  meta[ids,'dual'] = ifelse(meta[ids,'im_score'] > 0 & meta[ids,'dr_score'] > 0, 'Dual','Non-Dual')
  identical(rownames(meta),colnames(expr))
  sce = CreateSeuratObject(count = expr,
                           meta.data = meta)
  saveRDS(sce,'/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465_seurat.obj.rds')
}
# load data
if(T){
  GSE132465 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465_seurat.obj.rds')
  table(GSE132465$casee,GSE132465$Class)
  Normal = which(GSE132465$casee=='Normal' & GSE132465$Class == 'Normal') %>% names()
  Cancer = which(GSE132465$casee=='Cancer' & GSE132465$Class == 'Tumor') %>% names()
  length(Normal)
  length(Cancer)
  plot.path = 'script/R/GSE132465/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
}
# plot dual
if(F){
  df = GSE132465@meta.data[c(Normal, Cancer),]
  x= ggplot(df,aes(im_score,dr_score))+
    stat_density2d(geom ="polygon",aes(fill = ..level..),bins=25 )+
    scale_fill_gradientn(colours=c("#1B7EAB", "#ffffff" ,"#EF7B66"))+#, trans="log"
    theme_bw()+labs(title='')+ 
    geom_hline(aes(yintercept=0),color = '#990000',linetype='dashed')+
    geom_vline(aes(xintercept=0),color = '#990000',linetype='dashed')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+
    facet_wrap(~mut_type)
  x
  topptx(x,file="script/R/GSE132465/GSE132465_Dual-G1-4_grid-plot.pptx")
  
  library(ggpubr) 
  my_comparisons <- list( c("G1", "G2"), c("G1", "G3"), c("G1", "G4"))
  df$group = factor(df$mut_type,levels = c('G1','G2','G3','G4'))
  x = ggplot(df,aes(x=group,im_score,fill= group))+
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons,method = 't.test')+ # Add pairwise comparisons p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = ".all.", hide.ns = TRUE,label.y = 6) +
    labs(title='')+    theme(legend.position = "top") +scale_fill_manual(values = mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank())
  x
  topptx(x,file="script/R/GSE132465/GSE132465_Dual-immune_score_box-plot.pptx")
}
# check mut freq 
if(F){
  sce = subset(GSE132465, cells = Cancer)
  sce$group = sce$mut_type
  
  df = sce@meta.data
  df = df[,c('Sample','mut_type','mut_score',gene.list)]
  write.csv(df,'~/Desktop/GSE132465_COAD.csv')
  
  sce$sample_name
  result = NULL
  for(gene in gene.list){
    freq = as.matrix(table(sce@meta.data[['mut_type']],sce@meta.data[[gene]]))
    ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
    df = data.frame(MUT = freq[,1],
                    WT = freq[,2],
                    Freq = ids,
                    Gene = gene)
    df$group = rownames(df)
    result = rbind(result, df)
  }
  result$group = paste0(result$group,' (N=',result$MUT+result$WT,')')
  df = result[which(result$group=='G4 (N=4674)'),]
  ids = df[order(df$Freq),'Gene']
  result$Gene = factor(result$Gene, levels = ids)
  x = ggplot(result,aes(Gene, Freq,fill=group)) + scale_fill_manual(values = c("#F8766D","#ECA86E","#FFEF9C","#7FA4D1")) +
    geom_bar(stat = 'identity')+facet_wrap(~group,nrow = 4) +
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank()) +
    theme(axis.title =  element_text(size=6,face = "bold"),
          axis.text.x =   element_text(size=6,angle=90,  
                                       hjust = 1),
          axis.text.y = element_text(size=6)) 
  x
  topptx(x,file=paste0(plot.path,"GSE132465_Freq_bar-in-Cancer.pptx"),width = 10,height = 6)

}
# freq bar in dif group
if(F){
  sce = subset(GSE132465, cells = Cancer)
  mut.df = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465-mutation.csv',row.names = 1)
  df = as.data.frame.array(table(sce$Patient, sce$mut_type))
  ids = apply(df,1, sum)
  df$Freq_G1 = df$G1 / ids
  df$Freq_G2 = df$G2 / ids
  df$Freq_G3 = df$G3 / ids
  df$Freq_G4 = df$G4 / ids
  df = df[,c("Freq_G1","Freq_G2" ,"Freq_G3", 'Freq_G4')]
  colnames(df) = c('G1','G2','G3','G4')
  df$Cell_type = paste0(rownames(df),' Stage=',mut.df$Stage)
  df$Stage = mut.df$Stage
  df = df[order(df$Stage),]
  df$Cell_type = factor(df$Cell_type,levels = df$Cell_type)
  df = df[,1:5]
  ids = df
  df = reshape2::melt(df)
  x = ggplot(df,aes( value, Cell_type,fill=variable)) + 
    geom_bar(stat = 'identity') + scale_fill_manual(values =mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank(),
                                                                                   legend.position = "top")  
  x
  topptx(x,file=paste0(plot.path,"GSE132465_Freq_bar-in-patient.pptx"))
} # samples
if(F){
  sce = subset(GSE132465, cells = c(Cancer,Normal))
  result = NULL
  for(gene in gene.list){
    freq = as.data.frame.array(table(sce@meta.data[['casee']],sce@meta.data[[gene]]))
    ids1  = freq$`0`/ (freq[,1] + freq[,2])
    ids2  = freq$`1` / (freq[,1] + freq[,2])
    df = data.frame(WT = freq$`0`,
                    MUT = freq$`1`,
                    WT_Freq = ids1,
                    MUT_Freq = ids2,
                    Gene = gene) %>% set_rownames(rownames(freq))
    df$group = rownames(df)
    result = rbind(result, df)
  }
  #write.csv(result,'1.csv')
  ids1 = result[which(result$group == 'Cancer'),]
  ids2 = result[which(result$group == 'Normal'),]
  identical(ids1$Gene,ids2$Gene)
  df = data.frame(logFc = log(ids1$MUT_Freq/ids2$MUT_Freq,2),
                  Gene = ids1$Gene)
  df$dir = ifelse(df$logFc>0,'Up','Dn')
  df = df[order(df$logFc),]
  df$Gene = factor(df$Gene, levels = df$Gene)
  df = df[c(1:10,90:100),]
  x = ggplot(df,aes(logFc,Gene, fill=logFc)) + 
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-2.2,2.2))+
    geom_bar(stat = 'identity')
  x
  topptx(x,file=paste0(plot.path,"GSE132465_Freq_Ration_in_Cancer_vs_normal.pptx"))
  
}
# tsne & programs scores
if(F){
  sce = subset(GSE132465, cells = Cancer)
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  programs = read.csv('Cell-Human_program_file.csv')
  cols = colnames(programs)
  ids = cols
  for(program in cols){
    print(program)
    gene = as.list(programs[[program]])
    sce = AddModuleScore(
      object = sce,
      features = gene,
      name = program
    )
  }
  DimPlot(sce,group.by = 'mut_type')
  cols = paste0(colnames(programs),'1')
  df = sce@meta.data[,cols]
  th=1
  bk <- c(seq(-th,-0.1,by=0.01),seq(0,th,by=0.01))
  colnames(df) = ids
  pheatmap::pheatmap(cor(df),
                     color = c(colorRampPalette(colors = c("#3262A4","white"))(length(bk)/2),colorRampPalette(colors = c("white","#C91923"))(length(bk)/2)), # 设置颜色
                     legend_breaks=seq(-th,th,2), cellwidth = 25, cellheight = 25,
                     breaks=bk,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,legend = F,
                     filename = paste0(plot.path,"GSE132465_Cor_HeatMap.pdf"))
  
  df = sce@meta.data[which(sce$mut_type == 'G1'),cols]
  colnames(df) = ids
  pheatmap::pheatmap(cor(df),
                     color = c(colorRampPalette(colors = c("#3262A4","white"))(length(bk)/2),colorRampPalette(colors = c("white","#C91923"))(length(bk)/2)), # 设置颜色
                     legend_breaks=seq(-th,th,2), cellwidth = 25, cellheight = 25,
                     breaks=bk,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,legend = F,
                     filename = paste0(plot.path,"GSE132465_Cor_HeatMap_G1.pdf"))
  df = sce@meta.data[which(sce$mut_type == 'G4'),cols]
  colnames(df) = ids
  pheatmap::pheatmap(cor(df),
                     color = c(colorRampPalette(colors = c("#3262A4","white"))(length(bk)/2),colorRampPalette(colors = c("white","#C91923"))(length(bk)/2)), # 设置颜色
                     legend_breaks=seq(-th,th,2), cellwidth = 25, cellheight = 25,
                     breaks=bk,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,legend = F,
                     filename = paste0(plot.path,"GSE132465_Cor_HeatMap_G4.pdf"))
}
#  whole-genome sequencing vs RNA
if(F){
  sce = subset(GSE132465, cells = c(Cancer))
  sce = CreateSeuratObject(counts = sce@assays$RNA@counts, 
                           meta.data = sce@meta.data,
                           min.cells = 50)
  sce = NormalizeData(sce)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  
  mut.df = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465-mutation.csv',row.names = 1)
  ids = mut.df[,c('Stage','No.of.mutations')]
  ids$Patient = rownames(ids)
  sce@meta.data = merge(sce@meta.data,ids)
  
  df = sce@meta.data
  planes <- group_by(df, Patient) %>% summarise(.,
                                              No.of.mutations = median(No.of.mutations),
                                               mut_score = mean(mut_score)) 
  
  planes$No.of.mutations_log = log(planes$No.of.mutations,2)
  planes$mut_score_log = log(planes$mut_score,2)
  
  x = ggplot(planes,aes(mut_score,No.of.mutations_log))+geom_point()+
    geom_smooth(method = "glm",formula = y~I(x^5))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),legend.position = 'none')+
    stat_cor(method = "pearson",  label.x = 30,
             label.y =12)
  x
  topptx(x,file=paste0(plot.path,"GSE132465_WGS_RNA.pptx"))

  
}
# plot gene diversity
if(F){
  sce = subset(GSE132465, cells = c(Cancer,Normal))
  sce$group = ifelse(sce$casee == 'Normal',sce$casee,sce$mut_type)
  df = sce@meta.data
  df$group = factor(df$group,levels = c('G4','G3','G2','G1','Normal'))
  x= ggplot(df, aes(nFeature_RNA,group, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"GSE132465_gene_diversity.pptx")) 
  
  sce = subset(GSE132465)
  sce$group = ifelse(sce$Cell_type!='Epithelial cells',sce$Cell_type,ifelse(sce$casee == 'Normal',sce$casee,sce$mut_type))
  df = sce@meta.data
  df$group = factor(df$group,levels = c('G4','G3','G2','G1','Normal','Stromal cells','Myeloids','T cells',"B cells","Mast cells" ))
  df = df[which(!is.na(df$group)),]
  x= ggplot(df, aes(nFeature_RNA,group, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    #scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"GSE132465_allcells_diversity.pptx")) 
  
  x= ggplot(df, aes(nCount_RNA,group, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"GSE132465_count_diversity.pptx")) 
  
  x= ggplot(df, aes(nCount_RNA/nFeature_RNA,group, fill=group)) +
    geom_boxplot()+
    labs(title='')+
    theme(legend.position = "top")+
    scale_fill_manual(values = c(mut_type.cols,Normal='#d9e6eb'))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),
          axis.title.y  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"GSE132465_feature_count_diversity.pptx")) 
  
}
# CheckPoint expression in Cacner
if(F){
  sce = subset(GSE132465, cells = c(Cancer))
  sce$PD1_expr = sce@assays$RNA@counts['CD274',]
  sce$PD1_expr = ifelse(sce$PD1_expr>0,1,0)
  df = table(sce$mut_type,sce$PD1_expr) %>% as.data.frame.array()
  ids = apply(df,1,sum)
  df$Percent_CD274_Expression = df$`1`/ids *100
  df$group = c(1,2,3,4)
  x = ggplot(df,aes(group,Percent_CD274_Expression, size=Percent_CD274_Expression))+ geom_point() + 
    geom_smooth(method = "lm", formula = y ~ I(5^x))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
    stat_cor(method = "spearman")
  x
  topptx(x,file=paste0(plot.path,"GSE132465_PD1_expression_in_Cacner.pptx"))
  
  
  sce = subset(GSE132465, cells = c(Cancer))
  sce$CTLA4_expr = sce@assays$RNA@counts['CTLA4',]
  sce$CTLA4_expr = ifelse(sce$CTLA4_expr>0,1,0)
  df = table(sce$mut_type,sce$CTLA4_expr) %>% as.data.frame.array()
  ids = apply(df,1,sum)
  df$Percent_CTLA4_Expression = df$`1`/ids *100
  df$group = c(1,2,3,4)
  x = ggplot(df,aes(group,Percent_CTLA4_Expression, size=Percent_CTLA4_Expression))+ geom_point() + 
    geom_smooth(method = "lm")+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
    stat_cor(method = "pearson")
  x
  topptx(x,file=paste0(plot.path,"GSE132465_CTLA4_expression_in_Cacnery.pptx"))
}
# CheckPoint expression in TME
if(F){
  sce = GSE132465[,which(GSE132465$Cell_type %in% c("Myeloids","T cells") | colnames(GSE132465) %in% Cancer)]
  sce$PD1_expr = sce@assays$RNA@counts['CD274',]
  sce$PD1_expr = ifelse(sce$PD1_expr>0,1,0)
  df1 = table(sce$Patient,sce$PD1_expr) %>% as.data.frame.array()
  df2 = table(sce$Patient,sce$mut_type) %>% as.data.frame.array()
  df2$Freq_G4 = df2$G4 / apply(df2, 1, sum)
  df1$Percent_PD1_Expression = df1$`1`/apply(df1,1,sum) *100
  
  df = cbind(df1,df2)
  x = ggplot(df,aes(Freq_G4,Percent_PD1_Expression))+ geom_point() + 
    geom_smooth( method = 'lm')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
    stat_cor(method = "pearson")
  x
  topptx(x,file=paste0(plot.path,"GSE132465_PD1_expression_in_TME.pptx"))
  
}
# inferCNV
if(F){
  sce = subset(GSE132465,cells = Cancer)
  cells = sce[,which(sce$mut_type == 'G1' | sce$mut_type=='G4')] %>% colnames()
  sce = subset(sce, cells = cells)
  dfcount = as.data.frame(sce@assays$RNA@counts)
  groupinfo= data.frame(cellId = colnames(dfcount))
  groupinfo$cellType = sce@meta.data$mut_type
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
  
  expFile='script/R/GSE132465/GSE132465_G4_vs_G1_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='script/R/GSE132465/GSE132465_G4_vs_G1_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='script/R/GSE132465/GSE132465_G4_vs_G1_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  if(T){
    rm(list=ls())
    options(stringsAsFactors = F)
    library(Seurat)
    library(ggplot2)
    library(infercnv)
    library(future)
    plan("multicore", workers = 16)
    expFile='script/R/GSE132465/GSE132465_G4_vs_G1_expFile.txt'
    groupFiles='script/R/GSE132465/GSE132465_G4_vs_G1_groupFiles.txt'
    geneFile='script/R/GSE132465/GSE132465_G4_vs_G1_geneFile.txt'
    library(infercnv)
    gc()
    # infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
    #                                     annotations_file=groupFiles,
    #                                     delim="\t",
    #                                     gene_order_file= geneFile,
    #                                     ref_group_names=c('G1')) 
    # # Run infer CNV
    # infercnv_all = infercnv::run(infercnv_obj,
    #                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
    #                              out_dir= "InferCNV/GSE132465",  # dir is auto-created for storing outputs
    #                              cluster_by_groups=TRUE,   # cluster
    #                              num_threads=32,
    #                              denoise=T,
    #                              HMM=F)
    #saveRDS(infercnv_obj,'script/R/GSE132465/infercnv_obj.rds')
    #saveRDS(infercnv_all,'script/R/GSE132465/infercnv_all.rds')
    library(RColorBrewer)
    infercnv_all=readRDS('script/R/GSE132465/infercnv_all.rds')
    library(future)
    plan("multicore", workers = 16)
    infercnv::plot_cnv(infercnv_all, 
                       obs_title = "Observations (G4)",
                       ref_title = "References (G1)",
                       plot_chr_scale = T, contig_cex=2,
                       output_filename = "script/R/GSE132465/GSE132465_better_plot",output_format = "png")
  }
  
}
# monocle2
if(F){
  library(monocle)
  Mono_tj = subset(GSE132465, cells = Cancer) 
  immune.gene = read.csv('script/R/GeneList.csv')
  genes = unique(immune.gene$Symbol)
  genes = intersect(genes, rownames(Mono_tj))
  Mono_tj = Mono_tj[genes,]
  
  Mono_matrix<-as(as.matrix(GetAssayData(Mono_tj,slot = "counts")), 'sparseMatrix')
  feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
  rownames(feature_ann)<-rownames(Mono_matrix)
  Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
  sample_ann<-Mono_tj@meta.data
  Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
  Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())
  Mono.cds <- estimateSizeFactors(Mono.cds)
  Mono.cds <- estimateDispersions(Mono.cds)
  Mono.cds <- reduceDimension(
    Mono.cds,
    max_components = 2,
    method = 'DDRTree')
  Mono.cds <- orderCells(Mono.cds,num_paths=32)
  x = plot_cell_trajectory(Mono.cds,cell_size = 1,color_by='mut_type')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),
          legend.position = c(0.75,0.3)) +scale_color_manual(values = c(G1 = "#1a476f",G2 = "#698EC3",G3 = "#E64B35B2",G4 = "#c10534"))
  x
  topptx(x,file=paste0(plot.path,"GSE132465-immune-pseduo.pptx"))
}
#plot gsea G4-dual vs G1
if(F){
  sce = subset(GSE132465,cells = Cancer)
  g1 = which(sce$mut_type == 'G1') %>% names()
  g4 = which(sce$dual =='Dual' & sce$mut_type == 'G4') %>% names()
  sce = subset(sce,cells = c(g1,g4))
  table(sce$mut_type)
  
  sce = NormalizeData(sce,scale.factor = 1000000)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)

  degs = FindMarkers(sce, ident.1=g4, ident.2=g1)
  write.csv(degs,  paste0(plot.path,"GSE132465-dual-nodual.csv"))
  deg = degs
  library(EnhancedVolcano)
  x = EnhancedVolcano(deg,#xlim=c(-15,15),
                  lab = rownames(deg),
                  pCutoff = 0.01,
                  FCcutoff = 1.5,
                  x = 'avg_log2FC',
                  y = 'p_val_adj')
  topptx(x,file=paste0(plot.path,"GSE132465-G4-G1-Volcano.pptx"), width = 8,height = 5)
  
  geneList = deg$avg_log2FC
  names(geneList) = toupper(rownames(deg))
  geneList = sort(geneList, decreasing = T)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gmtfile = 'ref/h.all.v7.5.1.symbols.gmt'
  library(GSEABase)
  library(enrichplot)
  
  geneset = read.gmt(gmtfile=gmtfile)
  egmt = GSEA(geneList, TERM2GENE = geneset,
              minGSSize = 1,
              pvalueCutoff = 0.95)
  
  x = gseaplot2(egmt,color='#c10534',
                geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                subplots = 1:2)
  x
  topptx(x,file=paste0(plot.path,"GSE132465-G4-G1-GSEA-EMT.pptx"), width = 8,height = 5)
  
  gsea_result = egmt@result
  df = gsea_result[which(gsea_result$p.adjust < 0.05),]
  dim(df)
  df$Description = str_replace(df$Description,'HALLMARK_','')
  df$log_Pvalue = -log10(df$pvalue)
  df$group = ifelse(df$enrichmentScore>0,1,-1)
  df$log_Pvalue = df$log_Pvalue * df$group
  df = df[order(df$log_Pvalue),]
  df$Description = factor(df$Description, levels = df$Description)
  x = ggplot(df, aes(log_Pvalue,Description,fill = enrichmentScore)) + geom_bar(stat = 'identity')+
    scale_fill_gradientn(colours = c("#1a476f", "white","#c10534"),limits = c(-0.8,0.8),breaks=c(-0.8,-0.4,0,0.4,0.8))
  x
  topptx(x,file=paste0(plot.path,"GSE116237-T0_G4-T0_G1-GSEA-bar.pptx"), width = 8,height = 5)
}
# FindMarker genes
if(F){
  sce = subset(GSE132465,cells = Cancer)
  sce = NormalizeData(sce,scale.factor = 1000000)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  Idents(sce) = sce$mut_type
  markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  DoHeatmap(sce, features = top10$gene) + NoLegend()
  # IG
  if(F){
    gene = grep('^IG[^F]',rownames(sce@assays$RNA@scale.data),value = T) %>% list()
    sce = AddModuleScore(
      object = sce,
      features = gene,
      name = 'IG'
    )
    df = sce@meta.data
    planes <- group_by(df, Patient) %>% summarise(.,
                                                  IG = median(IG1),
                                                  mut_score = median(mut_score)) 
    x = ggplot(planes,aes(mut_score,IG))+geom_point()+
      geom_smooth(method = "lm")+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),legend.position = 'none')+
      stat_cor(method = "spearman",  label.x = 40,
               label.y = 0)
    x
    topptx(x,file=paste0(plot.path,"GSE132465-IG-score.pptx"))
    
  }
  gene = grep('^IG[^F]',rownames(sce@assays$RNA@scale.data),value = T)
  gene = c("IGHA1","IGHG1","IGHG3","IGHG4","IGKC","IGLC2","IGLC3")
  DotPlot(sce, features = gene,group.by = 'mut_type')
}
# Survival
if(F){
  library(survminer)
  library(survival)
  df = read.table('~/Desktop/GDC TCGA Colon Cancer (COAD)-PMP22.tsv',header = T,sep = '\t')
  df = df[which(df[,5]!='middle'),]
  df$surdate = df$OS.time/30
  df$event = df$OS
  plot_title_name='COAD-PMP22'
  annotate_x=40
  df$label = df$gene.expression.RNAseq...HTSeq...FPKM.UQ..PMP22
  df = na.omit(df)
  if(T){
    survival<-df[,c('surdate','event','label')]
    survival = na.omit(survival)
    label = survival$label
    label=as.matrix(label)
    survival = as.matrix(survival[,1:2])
    summary = summary(coxph(Surv(survival[,1], survival[,2]) ~ label))
    HR = round(summary$coefficients[2],2)
    CI95 = round(summary$conf.int[3:4],2)
    LogRankP = signif(summary$sctest[3],digits = 3)
    pp = paste("LogRank p = ",LogRankP)
    HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
    label[label==1,1] <- 'High-risk'
    label[label==0,1] <- 'Low-risk'
    xlab<-"Progression-Free Interval(months)"
    ylab<-"Survival probability probability"
    surv_TTP<-survfit(Surv(survival[,1], survival[,2]) ~ label)
    res <- ggsurvplot(fit = surv_TTP,data = as.data.frame(label), 
                      legend = "bottom",legend.labs=c("High-risk","Low-risk"),legend.title="",
                      size = 2.1,fontsize=7,palette = c("#99151B","#305596"),
                      ggtheme = theme_minimal(),
                      risk.table = T,risk.table.title = "Number at risk",
                      risk.table.y.text = F,risk.table.col="strata", cumevents = F,
                      #break.time.by = deg,
                      conf.int = F,surv.median.line = 'none',combine = T,newpage=T,
                      ncensor.plot=F,xlab = xlab , ylab = ylab
    )
    res
    res$table <- res$table + 
      theme(axis.line = element_blank(),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.x = element_text(family = "serif",face ="bold.italic",size=12 )
      )
    res$plot <- res$plot +  
      labs(title=plot_title_name) + 
      theme(axis.title.x = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.y = element_text(family = "serif",face ="bold.italic",size=12),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            axis.text.y = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(hjust = 0.5,size=12),
            legend.text = element_text(size = 12)
            
      )
    res$plot <- res$plot + 
      annotate("text", x = annotate_x,y= 0.15,
               label = pp,size = 4,family = "serif",
               colour = "black", fontface = "bold" #,fontface = "italic"
      ) + 
      annotate("text", x = annotate_x,y= 0.08,
               label = HHRR,size = 4,family = "serif",
               colour = "black",fontface = "bold" #,fontface = "italic",
      )
    print(res)
  }
  x = res$plot
  topptx(x,file=paste0(plot.path,"TCGA-PMP22_Survival-curve.pptx"))
  x = res$table
  topptx(x,file=paste0(plot.path,"TCGA-PMP22_Survival-table.pptx"))
  
  
  library(survminer)
  library(survival)
  df = read.table('~/Desktop/GDC TCGA Colon Cancer (COAD)-RGS2.tsv',header = T,sep = '\t')
  df = df[which(df[,5]!='middle'),]
  df$surdate = df$OS.time/30
  df$event = df$OS
  plot_title_name='COAD-RGS2'
  annotate_x=40
  df$label = df$gene.expression.RNAseq...HTSeq...FPKM.UQ..RGS2
  if(T){
    survival<-df[,c('surdate','event','label')]
    survival = na.omit(survival)
    label = survival$label
    label=as.matrix(label)
    survival = as.matrix(survival[,1:2])
    summary = summary(coxph(Surv(survival[,1], survival[,2]) ~ label))
    HR = round(summary$coefficients[2],2)
    CI95 = round(summary$conf.int[3:4],2)
    LogRankP = signif(summary$sctest[3],digits = 3)
    pp = paste("LogRank p = ",LogRankP)
    HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
    label[label==1,1] <- 'High-risk'
    label[label==0,1] <- 'Low-risk'
    xlab<-"Progression-Free Interval(months)"
    ylab<-"Survival probability probability"
    surv_TTP<-survfit(Surv(survival[,1], survival[,2]) ~ label)
    res <- ggsurvplot(fit = surv_TTP,data = as.data.frame(label), 
                      legend = "bottom",legend.labs=c("High-risk","Low-risk"),legend.title="",
                      size = 2.1,fontsize=7,palette = c("#99151B","#305596"),
                      ggtheme = theme_minimal(),
                      risk.table = T,risk.table.title = "Number at risk",
                      risk.table.y.text = F,risk.table.col="strata", cumevents = F,
                      #break.time.by = deg,
                      conf.int = F,surv.median.line = 'none',combine = T,newpage=T,
                      ncensor.plot=F,xlab = xlab , ylab = ylab
    )
    res
    res$table <- res$table + 
      theme(axis.line = element_blank(),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.x = element_text(family = "serif",face ="bold.italic",size=12 )
      )
    res$plot <- res$plot +  
      labs(title=plot_title_name) + 
      theme(axis.title.x = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.y = element_text(family = "serif",face ="bold.italic",size=12),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            axis.text.y = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(hjust = 0.5,size=12),
            legend.text = element_text(size = 12)
            
      )
    res$plot <- res$plot + 
      annotate("text", x = annotate_x,y= 0.15,
               label = pp,size = 4,family = "serif",
               colour = "black", fontface = "bold" #,fontface = "italic"
      ) + 
      annotate("text", x = annotate_x,y= 0.08,
               label = HHRR,size = 4,family = "serif",
               colour = "black",fontface = "bold" #,fontface = "italic",
      )
    print(res)
  }
  x = res$plot
  topptx(x,file=paste0(plot.path,"TCGA-RGS2_Survival-curve.pptx"))
  x = res$table
  topptx(x,file=paste0(plot.path,"TCGA-RGS2_Survival-table.pptx"))
  
}
# grade
if(F){
  sce = subset(GSE132465, cells = Cancer) 
  lung = sce@meta.data
  mut.df = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465-mutation.csv',row.names = 1)
  ids = mut.df[,c('Stage','No.of.mutations')]
  ids$Patient = rownames(ids)
  lung = merge(lung,ids)
  
  
  planes <- group_by(lung, Patient) %>% summarise(.,
                                                    grade = mean(Stage),
                                                    mut_score = median(mut_score)) 
  x = ggplot(planes,aes(grade,mut_score))+geom_point()+
    geom_smooth(method = "lm")+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),legend.position = 'none')+
    stat_cor(method = "spearman",  label.x = 1,
             label.y =1)
  x
  
}

############# PRJNA591860 LUAD ############# 
# prepare data
if(F){
  if(F){
    raw.data <- read.csv("/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_datafinal.csv", header=T, row.names = 1)
    dim(raw.data)
    head(colnames(raw.data))
    metadata <- read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_metacells.csv', row.names=1, header=T)
    head(metadata)
    erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
    percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
    ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
    raw.data <- raw.data[-ercc.index,]
    dim(raw.data)
    main_tiss <- CreateSeuratObject(counts = raw.data)
    # add rownames to metadta 
    row.names(metadata) <- metadata$cell_id
    # add metadata to Seurat object 
    main_tiss <- AddMetaData(object = main_tiss, metadata = metadata)
    main_tiss <- AddMetaData(object = main_tiss, percent.ercc, col.name = "percent.ercc")
    # Head to check
    head(main_tiss@meta.data)
    ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = main_tiss@assays$RNA@data), value = TRUE)
    percent.ribo <- Matrix::colSums(main_tiss@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(main_tiss@assays$RNA@data)
    main_tiss <- AddMetaData(object = main_tiss, metadata = percent.ribo, col.name = "percent.ribo")
    main_tiss
    sce <- subset(x=main_tiss, subset = nCount_RNA > 50000 & nFeature_RNA > 500)
    sce@meta.data$sample_name <- as.character(sce@meta.data$sample_name)
    sample_name <- as.character(sce@meta.data$sample_name)
    # Make table 
    tab.1 <- table(sce@meta.data$sample_name) 
    # Which samples have less than 10 cells 
    samples.keep <- names(which(tab.1 > 10))
    metadata_keep <- filter(sce@meta.data, sample_name %in% samples.keep)
    # Subset Seurat object 
    sce <- subset(sce, cells=as.character(metadata_keep$cell_id))
    sce
    table(sce@meta.data$sample_name)
    table(sce@meta.data$patient_id)
    saveRDS(sce,'/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_datafinal.rds')
    
  }
  sce = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_datafinal.rds')
  meta = sce@meta.data
  cancer = read.csv('/media/yuansh/14THHD/BscModel-V4/PRJNA591860/Candidate_cancer_predict_Integration_.csv',row.names = 1)
  cell.id = read.csv('Pan_cancer_Pathyway_detection/PRJNA591860/Pan_cancer_CNTNAP5/ScMutation_CNTNAP5_input_sce_embeddings.csv', row.names = 1) %>% rownames()
  cell.id = intersect(cell.id, colnames(sce))
  # get pathway detection result
  for(gene in gene.list){
    file.name = paste0('Pan_cancer_Pathyway_detection/PRJNA591860/Pan_cancer_',gene,'/','ScMutation_',gene,'_input_sce_embeddings.csv')
    assign(paste0('PRJNA591860_',gene), read.csv(file.name, row.names = 1)[cell.id,])
  }
  if(T){
    df = as.data.frame(cell.id)
    for(gene in gene.list){
      ids = get(paste0('PRJNA591860_',gene))
      df[[gene]] = ifelse(ids$predict == 'WT',0,1)
    }
    df$mut_score = apply(df[,gene.list],1,sum)
    summary(df$mut_score)
    df$mut_type = ifelse(df$mut_score<=30,'G1',
                         ifelse(df$mut_score<=55, 'G2',
                                ifelse(df$mut_score<=70, 'G3','G4')))
    table(df$mut_type)
    
  }
  
  col.names = c(gene.list,'mut_score','mut_type')
  ids = intersect(cell.id, rownames(meta))
  for(i in col.names){
    meta[[i]] = NA
    meta[ids,i] =df[,i]
  }
  
  meta[['casee']] = NA
  cancer = cancer[intersect(cell.id, rownames(cancer)),]
  meta[intersect(cell.id, rownames(cancer)),'casee'] = cancer$scale_predict
  
  dual = read.csv('script/R/PRJNA591860/PRJNA591860_Dual_cell_plot.csv',row.names = 1)
  dual = dual[ids,]
  meta$im_score = NA
  meta$dr_score = NA
  meta[ids,'im_score'] = dual$immune_score
  meta[ids,'dr_score'] = dual$driver_gene_score
  
  meta$dual = NA
  meta[ids,'dual'] = ifelse(meta[ids,'im_score'] > 0 & meta[ids,'dr_score'] > 0, 'Dual','Non-Dual')
  identical(rownames(meta),colnames(sce))
  sce@meta.data = meta
  saveRDS(sce,'/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_full_info.rds')
}
# load data
if(T){
  PRJNA591860 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_full_info.rds')
  table(PRJNA591860$casee)
  Normal = which(PRJNA591860$casee=='Normal') %>% names()
  Cancer = which(PRJNA591860$casee=='Cancer') %>% names()
  length(Normal)
  length(Cancer)
  plot.path = 'script/R/PRJNA591860/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
}
# dual  
if(F){
  df = PRJNA591860@meta.data
  df = df[which(!is.na(df$casee)),]
  df$group = df$mut_type
  colors = c("#1B7EAB", "#ffffff" ,"#EF7B66")
  x= ggplot(df,aes(im_score,dr_score))+
    stat_density2d(geom ="polygon",aes(fill = ..level..),bins=30 )+
    scale_fill_gradientn(colours=colors)+#, trans="log"
    theme_bw()+labs(title='')+ 
    geom_hline(aes(yintercept=0),color = '#990000',linetype='dashed')+
    geom_vline(aes(xintercept=0),color = '#990000',linetype='dashed')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+
    facet_wrap(~group)
  x
  topptx(x,file=paste0(plot.path,"PRJNA591860-dual.pptx"))
  
  my_comparisons <- list( c("G1", "G2"),
                          c("G1", "G3"), 
                          c("G1", "G4"))
  library(ggpubr)
  x = ggplot(df,aes(x=group,im_score,fill= group))+
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons,method = 't.test')+ # Add pairwise comparisons p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = ".all.", hide.ns = TRUE,label.y = 6) +
    labs(title='')+    theme(legend.position = "top") +scale_fill_manual(values = mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"PRJNA591860-box.pptx"))
}
# freq bar in dif group
if(F){
  sce = subset(PRJNA591860, cells = Cancer)
  df = as.data.frame.array(table(sce$sample_name, sce$mut_type))
  ids = apply(df,1, sum)
  df$Freq_G1 = df$G1 / ids
  df$Freq_G2 = df$G2 / ids
  df$Freq_G3 = df$G3 / ids
  df$Freq_G4 = df$G4 / ids
  df = df[,c("Freq_G1","Freq_G2" ,"Freq_G3", 'Freq_G4')]
  colnames(df) = c('G1','G2','G3','G4')
  df$Cell_type = rownames(df)
  df$Cell_type = factor(df$Cell_type,levels = df$Cell_type)
  df = df[,1:5]
  df = reshape2::melt(df)
  x = ggplot(df,aes( value, Cell_type,fill=variable)) + 
    geom_bar(stat = 'identity') + scale_fill_manual(values =mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank(),
                                                                                   legend.position = "top")  
  x
  topptx(x,file=paste0(plot.path,"PRJNA591860_Freq_bar-in-patient.pptx"))
} # samples
# Survival
if(F){
  library(survminer)
  library(survival)
  sce = subset(PRJNA591860, cells = Cancer)
  df.sur = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/Survival.csv',row.names = 1)
  df = as.data.frame.array(table(sce$patient_id, sce$mut_type))
  df$label = ifelse((df$G3+df$G4) /(df$G1+df$G2+df$G3+df$G4) >0.7,'High_Risk','Low_Risk')
  table(df$label)
  ids = intersect(rownames(df),rownames(df.sur))
  df = df[ids,]
  df.sur = df.sur[ids,]
  df$event = df.sur[,c('best_response_status')]
  df$surdate = df.sur[,c('pfs_month')]
  df = na.omit(df)
  plot_title_name = 'Survival'
  annotate_x <- 10
  df = df[which(df$G1+df$G2+df$G3+df$G4 > 10),]
  df$event = ifelse(df$event == 'PD',1,0)
  if(T){
    survival<-df[,c('surdate','event','label')]
    survival = na.omit(survival)
    label = survival$label
    label=as.matrix(label)
    survival = as.matrix(survival[,1:2])
    summary = summary(coxph(Surv(survival[,1], survival[,2]) ~ label))
    HR = round(summary$coefficients[2],2)
    CI95 = round(summary$conf.int[3:4],2)
    LogRankP = signif(summary$sctest[3],digits = 3)
    pp = paste("LogRank p = ",LogRankP)
    HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
    label[label==1,1] <- 'High-risk'
    label[label==0,1] <- 'Low-risk'
    xlab<-"Progression-Free Interval(months)"
    ylab<-"Survival probability probability"
    surv_TTP<-survfit(Surv(survival[,1], survival[,2]) ~ label)
    res <- ggsurvplot(fit = surv_TTP,data = as.data.frame(label), 
                      legend = "bottom",legend.labs=c("High-risk","Low-risk"),legend.title="",
                      size = 2.1,fontsize=7,palette = c("#99151B","#305596"),
                      ggtheme = theme_minimal(),
                      risk.table = T,risk.table.title = "Number at risk",
                      risk.table.y.text = F,risk.table.col="strata", cumevents = F,
                      #break.time.by = deg,
                      conf.int = F,surv.median.line = 'none',combine = T,newpage=T,
                      ncensor.plot=F,xlab = xlab , ylab = ylab
    )
    res$table <- res$table + 
      theme(axis.line = element_blank(),
            axis.text.x = element_text(family = "serif",face ="bold",size=15),
            plot.title = element_text(family = "serif",face ="bold.italic",size=20),
            axis.title.x = element_text(family = "serif",face ="bold.italic",size=20 )
      )
    res$plot <- res$plot +  
      labs(title=plot_title_name) + 
      theme(axis.title.x = element_text(family = "serif",face ="bold.italic",size=20),
            axis.title.y = element_text(family = "serif",face ="bold.italic",size=20),
            axis.text.x = element_text(family = "serif",face ="bold",size=15),
            axis.text.y = element_text(family = "serif",face ="bold",size=15),
            plot.title = element_text(hjust = 0.5,size=20),
            legend.text = element_text(size = 20)
            
      )
    res$plot <- res$plot + 
      annotate("text", x = annotate_x,y= 0.15,
               label = pp,size = 7,family = "serif",
               colour = "black", fontface = "bold" #,fontface = "italic"
      ) + 
      annotate("text", x = annotate_x,y= 0.08,
               label = HHRR,size = 7,family = "serif",
               colour = "black",fontface = "bold" #,fontface = "italic",
      )+ 
      annotate("text", x = annotate_x,y= 0.01,
               label = " C-index = 0.6933",size = 7,family = "serif",
               colour = "black",fontface = "bold" #,fontface = "italic",
      )
    print(res)
  }
  x = CombinePlots(plots = list(res$plot,res$table),ncol = 1)
  topptx(x,file=paste0(plot.path,"PRJNA591860_Survival.pptx"))
  
}
# CheckPoint expression in Cacner
if(F){
  sce = subset(PRJNA591860, cells = c(Cancer))
  sce$PD1_expr = sce@assays$RNA@counts['CD274',]
  sce$PD1_expr = ifelse(sce$PD1_expr>0,1,0)
  df = table(sce$mut_type,sce$PD1_expr) %>% as.data.frame.array()
  ids = apply(df,1,sum)
  df$Percent_PD1_Expression = df$`1`/ids *100
  df$group = c(1,2,3,4)
  x = ggplot(df,aes(group,Percent_PD1_Expression, size=Percent_PD1_Expression))+ geom_point() + 
    geom_smooth(method = "lm")+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
    stat_cor(method = "pearson")
  x
  topptx(x,file=paste0(plot.path,"PRJNA591860_PD1_expression_in_Cacner.pptx"))
}
# inferCNV
if(F){
  sce = subset(PRJNA591860,cells = Cancer)
  cells = sce[,which(sce$mut_type == 'G1' | sce$mut_type=='G4')] %>% colnames()
  sce = subset(sce, cells = cells)
  dfcount = as.data.frame(sce@assays$RNA@counts)
  groupinfo= data.frame(cellId = colnames(dfcount))
  groupinfo$cellType = sce@meta.data$mut_type
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
  
  expFile='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  if(T){
    rm(list=ls())
    options(stringsAsFactors = F)
    library(Seurat)
    library(ggplot2)
    library(infercnv)
    library(future)
    plan("multiprocess", workers = 16)
    expFile='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_expFile.txt'
    groupFiles='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_groupFiles.txt'
    geneFile='script/R/PRJNA591860/PRJNA591860_G4_vs_G1_geneFile.txt'
    library(infercnv)
    gc()
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                        annotations_file=groupFiles,
                                        delim="\t",
                                        gene_order_file= geneFile,
                                        ref_group_names=c('G1')) 
    # Run infer CNV
    infercnv_all = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir= "InferCNV/PRJNA591860",  # dir is auto-created for storing outputs
                                 cluster_by_groups=TRUE,   # cluster
                                 num_threads=32,
                                 denoise=T,
                                 HMM=F)
  }
}
# FindMarker genes
if(F){
  sce = subset(PRJNA591860,cells = Cancer)
  sce = NormalizeData(sce,scale.factor = 1000000)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  Idents(sce) = sce$mut_type
  markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  DoHeatmap(sce, features = top10$gene) + NoLegend()
  # IG
  if(F){
    gene = grep('^IG',rownames(sce@assays$RNA@scale.data),value = T) %>% list()
    gene
    sce = AddModuleScore(
      object = sce,
      features = gene,
      name = 'IG'
    )
    df = sce@meta.data
    planes <- group_by(df, patient_id) %>% summarise(.,
                                                  IG = median(IG1),
                                                  mut_score = median(mut_score)) 
    x = ggplot(planes,aes(mut_score,IG))+geom_point()+
      geom_smooth(method = "lm")+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),legend.position = 'none')+
      stat_cor(method = "spearman",  label.x = 50,
               label.y = 0)
    x
    topptx(x,file=paste0(plot.path,"GSE132465-IG-score.pptx"))
    
  }
  
}
# gender
if(F){
  sce = subset(PRJNA591860, cells = Cancer)
  df = sce@meta.data
  sce$CD274_expr.count = sce@assays$RNA@data['CD274',]
  sce$CD274_expr = ifelse(sce$CD274_expr.count>0,'Positive','Negative')
  df$CD274_expr = sce$CD274_expr
  df$CD274_expr.count = sce$CD274_expr.count
  planes <- group_by(df, sample_name) %>% summarise(.,
                                                  CD274_expr = mean(CD274_expr.count),
                                                  mut_score = median(mut_score)) 
  ids = df[,c('sample_name','gender')]
  male = ids[which(ids$gender=='Male'),'sample_name'] %>% unique()
  female = ids[which(ids$gender=='Female'),'sample_name'] %>% unique()
  planes$gender = ifelse(planes$sample_name %in%male,'male',ifelse(planes$sample_name %in%female,'Female',NA))
  x = ggplot(planes,aes(mut_score,CD274_expr))+geom_point()+
    geom_smooth(method = "lm")+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),legend.position = 'none')+
    stat_cor(method = "spearman",  label.x = 1,
             label.y =30) + facet_wrap(~gender)
  x
  ggplot(planes,aes(gender,CD274_expr,fill=gender))+geom_boxplot()
  topptx(x,file=paste0(plot.path,"GSE132465-Gender-PD-L1.pptx"))
  
}
############# GSE176078 BRCA ############# 
# prepare data
if(F){
  if(F){
    # expr = Matrix::readMM("../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx")
    # genes = read.table('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv')
    # barcode = read.table('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv')
    # 
    # colnames(expr) = barcode$V1
    # rownames(expr) = genes$V1
    # saveRDS(expr,'../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/raw_sce_count_matrix.rds')
    #sce = CreateSeuratObject(counts=expr, project = "GSE176078_BRCA")
    #saveRDS(sce,'../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/seurat.obj.rds')
    #sce = readRDS('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/seurat.obj.rds')
    #meta = read.csv('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/metadata.csv',row.names = 1)
    #sce = AddMetaData(sce, metadata = meta)
    #table(sce$celltype_major)
    
    #cell.id = names(which(sce$celltype_major == 'Cancer Epithelial'))
    #GSE176078.epi = subset(sce, cells = cell.id)
    
    #GSE176078.epi = CreateSeuratObject(counts=GSE176078.epi@assays$RNA@counts, 
    #                                   meta.data = GSE176078.epi@meta.data,
    #                                   project = "GSE176078_BRCA_Cancer",
    #                                   min.features = 150,
    #                                   min.cells = 10)
    #
    #write.csv(GSE176078.epi@assays$RNA@counts,'../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/GSE176078_BRCA_Cancer')
  }
  sce = readRDS('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/seurat.obj.rds')
  meta = read.csv('../Pan_Cancer_landscape/data/GSE176078_BRCA/Wu_etal_2021_BRCA_scRNASeq/metadata.csv',row.names = 1)
  sce = AddMetaData(sce, metadata = meta)
  
  cell.id = read.csv('Pan_cancer_Pathyway_detection/GSE176078/Pan_cancer_CNTNAP5/ScMutation_CNTNAP5_input_sce_embeddings.csv', row.names = 1) %>% rownames()
  cell.id = intersect(cell.id, colnames(sce))
  
  for(gene in gene.list){
    file.name = paste0('Pan_cancer_Pathyway_detection/GSE176078/Pan_cancer_',gene,'/','ScMutation_',gene,'_input_sce_embeddings.csv')
    assign(paste0('GSE176078_',gene), read.csv(file.name, row.names = 1)[cell.id,])
  }
  
  if(T){
    df = as.data.frame(cell.id)
    for(gene in gene.list){
      ids = get(paste0('GSE176078_',gene))
      df[[gene]] = ifelse(ids$predict == 'WT',0,1)
    }
    df$mut_score = apply(df[,gene.list],1,sum)
    summary(df$mut_score)
    df$mut_type = ifelse(df$mut_score<=22,'G1',
                         ifelse(df$mut_score<=29, 'G2',
                                ifelse(df$mut_score<=36.00, 'G3','G4')))
    table(df$mut_type)
  }
  
  col.names = c(gene.list,'mut_score','mut_type')
  ids = intersect(cell.id, rownames(meta))
  for(i in col.names){
    meta[[i]] = NA
    meta[ids,i] =df[,i]
  }
  dual = read.csv('script/R/GSE176078/GSE176078_Dual_cell_plot.csv',row.names = 1)
  dual = dual[ids,]
  meta$im_score = NA
  meta$dr_score = NA
  meta[ids,'im_score'] = dual$immune_score
  meta[ids,'dr_score'] = dual$driver_gene_score
  
  meta$dual = NA
  meta[ids,'dual'] = ifelse(meta[ids,'im_score'] > 0 & meta[ids,'dr_score'] > 0, 'Dual','Non-Dual')
  identical(rownames(meta),colnames(sce))
  sce@meta.data = meta
  saveRDS(sce,'/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE176078_BRCA/GSE176078_seurat_obj.rds')
}
# load data
if(F){
  GSE176078 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE176078_BRCA/GSE176078_seurat_obj.rds')
  table(GSE176078$celltype_major)
  Normal = which(GSE176078$celltype_major=='Normal Epithelial') %>% names()
  Cancer = which(GSE176078$celltype_major=='Cancer Epithelial') %>% names()
  length(Normal)
  length(Cancer)
  plot.path = 'script/R/GSE176078/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
}
# inferCNV
if(F){
  sce = subset(GSE176078,cells = Cancer)
  cells = sce[,which(sce$mut_type == 'G1' | sce$mut_type=='G4')] %>% colnames()
  sce = subset(sce, cells = cells)
  dfcount = as.data.frame(sce@assays$RNA@counts)
  groupinfo= data.frame(cellId = colnames(dfcount))
  groupinfo$cellType = sce@meta.data$mut_type
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
  
  expFile='script/R/GSE176078/GSE176078_G4_vs_G1_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='script/R/GSE176078/GSE176078_G4_vs_G1_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='script/R/GSE176078/GSE176078_G4_vs_G1_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  if(T){
    rm(list=ls())
    options(stringsAsFactors = F)
    library(Seurat)
    library(ggplot2)
    library(infercnv)
    library(future)
    plan("multiprocess", workers = 16)
    expFile='script/R/GSE176078/GSE176078_G4_vs_G1_expFile.txt'
    groupFiles='script/R/GSE176078/GSE176078_G4_vs_G1_groupFiles.txt'
    geneFile='script/R/GSE176078/GSE176078_G4_vs_G1_geneFile.txt'
    library(infercnv)
    gc()
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                        annotations_file=groupFiles,
                                        delim="\t",
                                        gene_order_file= geneFile,
                                        ref_group_names=c('G1')) 
    # Run infer CNV
    infercnv_all = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir= "InferCNV/GSE176078",  # dir is auto-created for storing outputs
                                 cluster_by_groups=TRUE,   # cluster
                                 num_threads=32,
                                 denoise=T,
                                 HMM=F)
  }
}
# FindMarker genes
if(F){
  sce = subset(GSE176078,cells = Cancer)
  sce = sce[,which(sce$subtype == 'TNBC')]
  sce = NormalizeData(sce,scale.factor = 1000000)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  Idents(sce) = sce$mut_type
  markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  DoHeatmap(sce, features = top10$gene) + NoLegend()
  # IG
  if(F){
    gene = grep('^IL|^IG',rownames(sce),value = T) %>% list()
    gene
    sce = AddModuleScore(
      object = sce,
      features = gene,
      name = 'CCL'
    )
    df = sce@meta.data
    planes <- group_by(df, orig.ident) %>% summarise(.,
                                                     CCL = median(CCL1),
                                                     mut_score = median(mut_score)) 
    x = ggplot(planes,aes(mut_score,CCL))+geom_point()+
      geom_smooth(method = "lm")+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),legend.position = 'none')+
      stat_cor(method = "spearman",  label.x = 30,
               label.y = 0)
    x
    topptx(x,file=paste0(plot.path,"GSE132465-IG-score.pptx"))
    
  }
  
}
# dual  
if(F){
  df = GSE176078@meta.data
  df = df[which(df$celltype_major=='Cancer Epithelial'),]
  df$group = df$mut_type
  colors = c("#1B7EAB", "#ffffff" ,"#EF7B66")
  x= ggplot(df,aes(im_score,dr_score))+
    stat_density2d(geom ="polygon",aes(fill = ..level..),bins=30 )+
    scale_fill_gradientn(colours=colors)+#, trans="log"
    theme_bw()+labs(title='')+ 
    geom_hline(aes(yintercept=0),color = '#990000',linetype='dashed')+
    geom_vline(aes(xintercept=0),color = '#990000',linetype='dashed')+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+
    facet_wrap(~group)
  x
  topptx(x,file=paste0(plot.path,"GSE176078-dual.pptx"))
  
  my_comparisons <- list( c("G1", "G2"),
                          c("G1", "G3"), 
                          c("G1", "G4"))
  library(ggpubr)
  x = ggplot(df,aes(x=group,im_score,fill= group))+
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons,method = 't.test')+ # Add pairwise comparisons p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = ".all.", hide.ns = TRUE,label.y = 6) +
    labs(title='')+    theme(legend.position = "top") +scale_fill_manual(values = mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank())
  x
  topptx(x,file=paste0(plot.path,"PRJNA591860-box.pptx"))
}
############# GSE186344 Brain Metastases ############# 
# prepare data
if(F){
  sce = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE186344_BRAINMetas/seurat_object/raw_merge.RDS')
  meta = sce@meta.data
  cell.id = read.csv('Pan_cancer_Pathyway_detection/GSE186344/Pan_cancer_CNTNAP5/ScMutation_CNTNAP5_input_sce_embeddings.csv', row.names = 1) %>% rownames()
  cell.id = gsub('\\.','-',cell.id)
  cell.id = intersect(cell.id, colnames(sce))
  for(gene in gene.list){
    file.name = paste0('Pan_cancer_Pathyway_detection/GSE186344/Pan_cancer_',gene,'/','ScMutation_',gene,'_input_sce_embeddings.csv')
    df = read.csv(file.name, row.names = 1)
    rownames(df) = gsub('\\.','-',rownames(df))
    assign(paste0('GSE186344_',gene), df[cell.id,])
  }
  
  if(T){
    df = as.data.frame(cell.id)
    for(gene in gene.list){
      ids = get(paste0('GSE186344_',gene))
      df[[gene]] = ifelse(ids$predict == 'WT',0,1)
    }
    df$mut_score = apply(df[,gene.list],1,sum)
    summary(df$mut_score)
    df$mut_type = ifelse(df$mut_score<=13,'G1',
                         ifelse(df$mut_score<=25, 'G2',
                                ifelse(df$mut_score<=70, 'G3','G4')))
    table(df$mut_type)
    
  }
  
  col.names = c(gene.list,'mut_score','mut_type')
  ids = intersect(cell.id, rownames(meta))
  for(i in col.names){
    meta[[i]] = NA
    meta[ids,i] =df[,i]
  }
  identical(rownames(meta),colnames(sce))
  sce@meta.data = meta
  saveRDS(sce,'/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE186344_BRAINMetas/seurat_object/GSE186344.ScMutation.rds')
}
# load data
if(F){
  GSE186344 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE186344_BRAINMetas/seurat_object/GSE186344.ScMutation.rds')
  table(GSE186344$Cell_Type)
  table(GSE186344$orig.ident)
  Cancer_type = data.frame(SampleID = c('GSM5645888', 'GSM5645889', 'GSM5645890', 
                                        'GSM5645891', 'GSM5645892', 'GSM5645893','GSM5645908',
                                        'GSM5645894', 'GSM5645895', 'GSM5645896', 
                                        'GSM5645897', 'GSM5645898',
                                        'GSM5645900', 
                                        'GSM5645902', 
                                        'GSM5645904'),
                           Groups = c('Melan','Melan','Melan',
                                      'Breast','Breast','Breast','Breast',
                                      'Lung','Lung','Lung',
                                      'Ovarian','Ovarian',
                                      'Colorectal',
                                      'Renal',
                                      'Unknown'))
  GSE186344$Cancer_type = ifelse(GSE186344$orig.ident %in% c('GSM5645888', 'GSM5645889', 'GSM5645890'),'Melan',
                                 ifelse(GSE186344$orig.ident %in% c( 'GSM5645891', 'GSM5645892', 'GSM5645893','GSM5645908'),'Breast',
                                        ifelse(GSE186344$orig.ident %in% c('GSM5645894', 'GSM5645895', 'GSM5645896'), 'Lung',
                                               ifelse(GSE186344$orig.ident %in% c('GSM5645897', 'GSM5645898'),'Ovarian',
                                                      ifelse(GSE186344$orig.ident == 'GSM5645900','Colorectal',
                                                             ifelse(GSE186344$orig.ident == 'GSM5645902','Renal',
                                                                    ifelse(GSE186344$orig.ident == 'GSM5645904','Unknown',-1))))))) 
  
  table(GSE186344$Cancer_type)
  Cancer = which(GSE186344$Cell_Type=='MTC') %>% names()
  length(Cancer)
  plot.path = 'script/R/GSE186344/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
}
# inferCNV
if(F){
  sce = subset(GSE186344,cells = Cancer)
  cells = sce[,which(sce$mut_type == 'G1' | sce$mut_type=='G4')] %>% colnames()
  sce = subset(sce, cells = cells)
  dfcount = as.data.frame(sce@assays$RNA@counts)
  groupinfo= data.frame(cellId = colnames(dfcount))
  groupinfo$cellType = sce@meta.data$mut_type
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
  
  expFile='script/R/GSE186344/GSE186344_G4_vs_G1_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='script/R/GSE186344/GSE186344_G4_vs_G1_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='script/R/GSE186344/GSE186344_G4_vs_G1_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
  if(T){
    rm(list=ls())
    options(stringsAsFactors = F)
    library(Seurat)
    library(ggplot2)
    library(infercnv)
    library(future)
    plan("multicore", workers = 16)
    expFile='script/R/GSE186344/GSE186344_G4_vs_G1_expFile.txt'
    groupFiles='script/R/GSE186344/GSE186344_G4_vs_G1_groupFiles.txt'
    geneFile='script/R/GSE186344/GSE186344_G4_vs_G1_geneFile.txt'
    library(infercnv)
    gc()
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                        annotations_file=groupFiles,
                                        delim="\t",
                                        gene_order_file= geneFile,
                                        ref_group_names=c('G1')) 
    # Run infer CNV
    infercnv_all = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir= "InferCNV/GSE186344",  # dir is auto-created for storing outputs
                                 cluster_by_groups=TRUE,   # cluster
                                 num_threads=32,
                                 denoise=T,
                                 HMM=F)
    
    #saveRDS(infercnv_obj,'script/R/GSE186344/infercnv_obj.rds')
    #saveRDS(infercnv_all,'script/R/GSE186344/infercnv_all.rds')
    library(RColorBrewer)
    infercnv_all=readRDS('script/R/GSE186344/infercnv_all.rds')
    infercnv::plot_cnv(infercnv_all, 
                       obs_title = "Observations (G4)",
                       ref_title = "References (G1)",
                       plot_chr_scale = T, 
                       output_filename = "script/R/GSE186344/GSE186344_better_plot",output_format = "pdf",
                       custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD")))
    
  }
}
# freq in cancer type 
if(F){
  sce = subset(GSE132465, cells = Cancer)
  df = as.data.frame.array(table(sce$Cancer_type, sce$mut_type))
  ids = apply(df,1, sum)
  df$Freq_G1 = df$G1 / ids
  df$Freq_G2 = df$G2 / ids
  df$Freq_G3 = df$G3 / ids
  df$Freq_G4 = df$G4 / ids
  df = df[,c("Freq_G1","Freq_G2" ,"Freq_G3", 'Freq_G4')]
  colnames(df) = c('G1','G2','G3','G4')
  df$Cancer_type = rownames(df)
  df = df[order(df$G4,decreasing = T),]
  df$Cancer_type = factor(df$Cancer_type,levels = df$Cancer_type)
  df = df[,1:5]
  ids = df
  df = reshape2::melt(df)
  x = ggplot(df,aes( value, Cancer_type,fill=variable)) + 
    geom_bar(stat = 'identity') + scale_fill_manual(values =mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank(),
                                                                                   legend.position = "top")

  x
  topptx(x,file=paste0(plot.path,"GSE186344_Freq_bar-in-Cancer_type.pptx"),width = 15,height = 10)
}
############# Pancancer ############# 
# prepare data
if(F){
  GSE132465 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE132465_COAD/GSE132465_seurat.obj.rds')
  Cancer = which(GSE132465$casee=='Cancer' & GSE132465$Class == 'Tumor') %>% names()
  GSE132465 = subset(GSE132465, cells = Cancer)
  
  PRJNA591860 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/S01_full_info.rds')
  Cancer = which(PRJNA591860$casee=='Cancer') %>% names()
  PRJNA591860 = subset(PRJNA591860, cells = Cancer)
  
  GSE176078 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE176078_BRCA/GSE176078_seurat_obj.rds')
  Cancer = which(GSE176078$celltype_major=='Cancer Epithelial') %>% names()
  GSE176078 = subset(GSE176078, cells = Cancer)
  
  GSE186344 = readRDS('/media/yuansh/14THHD/Pan_Cancer_landscape/data/GSE186344_BRAINMetas/seurat_object/GSE186344.ScMutation.rds')
  Cancer_type = data.frame(SampleID = c('GSM5645888', 'GSM5645889', 'GSM5645890', 'GSM5645891', 'GSM5645892', 'GSM5645893','GSM5645908','GSM5645894', 'GSM5645895', 'GSM5645896', 'GSM5645897', 'GSM5645898','GSM5645900', 'GSM5645902', 'GSM5645904'),
                           Groups = c('Melan','Melan','Melan','Breast','Breast','Breast','Breast','Lung','Lung','Lung','Ovarian','Ovarian','Colorectal','Renal','Unknown'))
  GSE186344$Cancer_type = ifelse(GSE186344$orig.ident %in% c('GSM5645888', 'GSM5645889', 'GSM5645890'),'Melan',ifelse(GSE186344$orig.ident %in% c( 'GSM5645891', 'GSM5645892', 'GSM5645893','GSM5645908'),'Breast',ifelse(GSE186344$orig.ident %in% c('GSM5645894', 'GSM5645895', 'GSM5645896'), 'Lung',ifelse(GSE186344$orig.ident %in% c('GSM5645897', 'GSM5645898'),'Ovarian',ifelse(GSE186344$orig.ident == 'GSM5645900','Colorectal',ifelse(GSE186344$orig.ident == 'GSM5645902','Renal',ifelse(GSE186344$orig.ident == 'GSM5645904','Unknown',-1))))))) 
  
  Cancer = which(GSE186344$Cell_Type=='MTC') %>% names()
  GSE186344 = subset(GSE186344, cells = Cancer)
  
  plot.path = 'script/R/Pancancer/'
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
  
  GSE132465$Dataset = 'GSE132465'
  PRJNA591860$Dataset = 'PRJNA591860'
  GSE176078$Dataset = 'GSE176078'
  GSE186344$Dataset = 'GSE186344'
  sce =merge(GSE132465,y=c(PRJNA591860,GSE176078,GSE186344))
  saveRDS(sce,'script/R/Pancancer/pancancer.rds')
  
}
# load data
if(T){
  pancancer = readRDS('script/R/Pancancer/pancancer.rds')
  plot.path = 'script/R/Pancancer/'
  
  Cancer_type = data.frame(SampleID = c('GSM5645888', 'GSM5645889', 'GSM5645890', 'GSM5645891', 'GSM5645892', 'GSM5645893','GSM5645908','GSM5645894', 'GSM5645895', 'GSM5645896', 'GSM5645897', 'GSM5645898','GSM5645900', 'GSM5645902', 'GSM5645904'),
                           Groups = c('Melan','Melan','Melan','Breast','Breast','Breast','Breast','Lung','Lung','Lung','Ovarian','Ovarian','Colorectal','Renal','Unknown'))
  data_info = data.frame(DataSet = c('GSE132465','PRJNA591860','GSE176078'),
                         Groups = c('Colorectal','Lung','Breast'))
  
  pancancer$Cancer_type = ifelse(!is.na(pancancer$Cancer_type),pancancer$Cancer_type,
                                 ifelse(pancancer$Dataset == 'GSE132465','Colorectal',
                                        ifelse(pancancer$Dataset == 'PRJNA591860','Lung',
                                               ifelse(pancancer$Dataset == 'GSE176078','Breast',NA))))
  
  pancancer = pancancer[,which(pancancer$Cancer_type%in%c('Breast', 'Colorectal',  'Lung', 'Melan', 'Ovarian'))]
  pancancer$pan.type = ifelse(pancancer$mut_score<=15,'P1',
                        ifelse(pancancer$mut_score<=27, 'P2',
                               ifelse(pancancer$mut_score<=47, 'P3','P4')))
  pancancer$SampleIds = ifelse(pancancer$orig.ident == 'SeuratProject',pancancer$sample_name,pancancer$orig.ident)
}
# TSNE
if(T){
  data.cols = c(GSE132465='#a0522d',PRJNA591860='#ECA86E',GSE176078='#DC0000B2',GSE186344='#7FA4D1')
  cancer.cls = c(Breast='#DC0000B2',Ovarian='#E64B35B2',Colorectal='#a0522d',Melan='#1a476f',Lung='#ECA86E')
  sce = CreateSeuratObject(counts = pancancer@assays$RNA@counts,
                           meta.data = pancancer@meta.data,
                           min.cells = 200)
  
  sce = NormalizeData(sce,scale.factor = 1000000)
  sce = FindVariableFeatures(sce)
  sce = ScaleData(object = sce)
  sce = RunPCA(object = sce, do.print = FALSE)
  sce= FindNeighbors(sce, dims = 1:20)
  sce = FindClusters(sce) 
  sce=RunUMAP(sce,dims = 1:20)
  DimPlot(sce,group.by = 'Dataset',cols = data.cols)+theme(legend.position = c(0.8,0.15),
                                                                plot.title = element_blank(),
                                                                axis.title = element_blank(),
                                                                axis.text = element_blank(),
                                                                axis.ticks = element_blank())
  DimPlot(sce,group.by = 'SampleIds')+ theme(legend.position="bottom")
}
# Redefine groups
if(T){
  df = sce@meta.data
  ggplot(df,aes(nFeature_SCT,Cancer_type,fill= Cancer_type))+
    geom_boxplot()+
    labs(title='')+
    theme()+
    scale_fill_manual(values = cancer.cls)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
          legend.title=element_blank(),
          axis.title.x  = element_blank(),
          axis.title.y  = element_blank(),
          legend.position = "top")
  
  df = sce@meta.data
  
  df[which(df$Cancer_type=='Colorectal'),'mut_score'] %>% summary()
  df[which(df$Cancer_type=='Lung'),'mut_score'] %>% summary()
  df[which(df$Cancer_type=='Breast'),'mut_score'] %>% summary()
  df[which(df$Cancer_type=='Melan'),'mut_score'] %>% summary()
  df[which(df$Cancer_type=='Ovarian'),'mut_score'] %>% summary()
  
  df[which(df$Cancer_type=='Colorectal'),]$mut_type = ifelse(df[which(df$Cancer_type=='Colorectal'),]$mut_score<=10,'G1',
                                                             ifelse(df[which(df$Cancer_type=='Colorectal'),]$mut_score<=25, 'G2',
                                                                    ifelse(df[which(df$Cancer_type=='Colorectal'),]$mut_score<=45, 'G3','G4')))
  
  df[which(df$Cancer_type=='Lung'),]$mut_type = ifelse(df[which(df$Cancer_type=='Lung'),]$mut_score<=35,'G1',
                                                       ifelse(df[which(df$Cancer_type=='Lung'),]$mut_score<=55, 'G2',
                                                              ifelse(df[which(df$Cancer_type=='Lung'),]$mut_score<=75, 'G3','G4')))
  
  df[which(df$Cancer_type=='Breast'),]$mut_type = ifelse(df[which(df$Cancer_type=='Breast'),]$mut_score<=13,'G1',
                                                         ifelse(df[which(df$Cancer_type=='Breast'),]$mut_score<=22, 'G2',
                                                                ifelse(df[which(df$Cancer_type=='Breast'),]$mut_score<=33, 'G3','G4')))
  
  df[which(df$Cancer_type=='Melan'),] $mut_type = ifelse(df[which(df$Cancer_type=='Melan'),]$mut_score<=10,'G1',
                                                         ifelse(df[which(df$Cancer_type=='Melan'),]$mut_score<=13, 'G2',
                                                                ifelse(df[which(df$Cancer_type=='Melan'),]$mut_score<=22, 'G3','G4')))
  
  df[which(df$Cancer_type=='Ovarian'),] $mut_type = ifelse(df[which(df$Cancer_type=='Ovarian'),]$mut_score<=15,'G1',
                                                           ifelse(df[which(df$Cancer_type=='Ovarian'),]$mut_score<=25, 'G2',
                                                                  ifelse(df[which(df$Cancer_type=='Ovarian'),]$mut_score<=43, 'G3','G4')))
  
  table(df$Cancer_type,df$mut_type)
  
  identical(rownames(df),rownames(sce@meta.data))
}
meta = df
library(ggbump)

x = ggplot(meta,aes(mut_type,mut_score,fill = mut_type))+
  geom_boxplot()+facet_wrap(~Cancer_type)+ scale_fill_manual(values =mut_type.cols)+
  theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                 axis.title.x  = element_blank(),
                                                                                 axis.title.y  = element_blank(),
                                                                                 legend.position = "top")

topptx(x,file=paste0(plot.path,"Pancancer_mut_type.pptx"))

#ggplot(meta[which(meta$Dataset!='GSE186344'),],aes(mut_type,mut_score,fill = mut_type))+geom_boxplot()+facet_wrap(~Cancer_type)
#ggplot(meta[which(meta$Dataset=='GSE186344'),],aes(mut_type,mut_score,fill = mut_type))+geom_boxplot()+facet_wrap(~Cancer_type)

sce@meta.data = meta
# check mut freq & TCGA
if(F){
  tcga = read.csv('TCGA-MUT-BURD.csv') %>% na.omit()
  
  tcga$cancer = ifelse(tcga$DISEASE == 'BRCA','Breast',ifelse(tcga$DISEASE == 'COAD','Colorectal',ifelse(tcga$DISEASE == 'OV','Ovarian',ifelse(tcga$DISEASE == 'SKCM','Melan',ifelse(tcga$DISEASE == "LUAD" | tcga$DISEASE == "LUSC" ,'Lung',NA)))))
  tcga$cancer = factor(tcga$cancer, levels = c("Breast","Ovarian","Colorectal","Lung","Melan"))
  pancancer$Cancer_type = factor(pancancer$Cancer_type, levels = c("Breast","Ovarian","Colorectal","Lung","Melan"))
  
  cls = c(Breast='#DC0000B2',Ovarian='#E64B35B2',Colorectal='#a0522d',Melan='#1a476f',Lung='#3C5488B2')
  # check mut freq 
  if(F){
    
    result = NULL
    for(gene in gene.list){
      freq = as.matrix(table(sce@meta.data[['Cancer_type']],sce@meta.data[[gene]]))
      ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
      df = data.frame(MUT = freq[,1],
                      WT = freq[,2],
                      Freq = ids,
                      Gene = gene)
      df$group = rownames(df)
      result = rbind(result, df)
    }
    result$group.no = paste0(result$group,' (N=',result$MUT+result$WT,')')
    result[which(result$group=='Breast'),'Freq'] %>% summary()
    result[which(result$group=='Colorectal'),'Freq'] %>% summary()
    result[which(result$group=='Lung'),'Freq'] %>% summary()
    result[which(result$group=='Melan'),'Freq'] %>% summary()
    result[which(result$group=='Ovarian'),'Freq'] %>% summary()
    
    result[which(result$group=='Breast'),'Freq'] %>% var()
    result[which(result$group=='Colorectal'),'Freq'] %>% var()
    result[which(result$group=='Lung'),'Freq'] %>% var()
    result[which(result$group=='Melan'),'Freq'] %>% var()
    result[which(result$group=='Ovarian'),'Freq'] %>% var()
    
    #result$Gene = factor(result$Gene, levels = ids)
    x = ggplot(result,aes(Gene, Freq,fill=group)) + scale_fill_manual(values = cls) +
      geom_bar(stat = 'identity')+facet_wrap(~group,nrow = 5) +
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                     axis.title.x  = element_blank(),
                                                                                     axis.title.y  = element_blank()) +
      theme(axis.title =  element_text(size=6,face = "bold"),
            axis.text.x =   element_text(size=6,angle=90,  
                                         hjust = 1),
            axis.text.y = element_text(size=6)) 
    x
    topptx(x,file=paste0(plot.path,"Pancancer_Freq_bar-in-Cancer.pptx"),width = 10,height = 6)
    
  }
  ggplot(tcga,aes(cancer,log10_mut, fill=cancer))+ 
    geom_boxplot() + 
    scale_fill_manual(values = cls) +   
    theme(legend.position = "top")
  
  ggplot(pancancer@meta.data,aes(Cancer_type,mut_score, fill=Cancer_type))+ 
    geom_boxplot() + 
    scale_fill_manual(values = cls) +   
    theme(legend.position = "top")
}
# CheckPoint expression in Cacner
if(F){
  if(F){
    df = meta
    sce$CD274_expr.count = sce@assays$RNA@data['CD274',]
    sce$CD274_expr = ifelse(sce$CD274_expr.count>0,'Positive','Negative')
    df$CD274_expr = sce$CD274_expr
    df$CD274_expr.count = sce$CD274_expr.count
    planes <- group_by(df, SampleIds) %>% summarise(.,
                                                    CD274_expr = mean(CD274_expr.count),
                                                    mut_score = median(mut_score)) 
    x = ggplot(planes,aes(mut_score,CD274_expr))+geom_point()+
      geom_smooth(method = "lm")+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),legend.position = 'none')+
      stat_cor(method = "spearman",  label.x = 1,
               label.y =1)
    x
    topptx(x,file=paste0(plot.path,"Pancancer_CD274_expression.pptx"))
    
    
    
    df1 = table(df$SampleIds,df$CD274_expr) %>% as.data.frame.array()
    df2 = table(df$SampleIds,df$mut_type) %>% as.data.frame.array()
    df2$Freq_G4 = df2$G4 / apply(df2, 1, sum)
    df1$Percent_CD274_Expression = df1$Positive/apply(df1,1,sum) *100
    ids = cbind(df1,df2)
    ids$SampleIds = rownames(ids)
    samples = unique(df[,c('SampleIds','Cancer_type','Dataset')])
    
    ids = merge(ids,samples)
    
    x1 = ggplot(ids[which(ids$Cancer_type == 'Breast'),],aes(Freq_G4,Percent_CD274_Expression))+ geom_point()+ 
      geom_smooth( method = 'lm')+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
      stat_cor(method = "pearson")
    x1

    x2 = ggplot(ids[which(ids$Cancer_type == 'Lung'),],aes(Freq_G4,Percent_CD274_Expression))+ geom_point()+ 
      geom_smooth( method = 'lm')+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
      stat_cor(method = "pearson")
    x2
    
    x3 = ggplot(ids[which(ids$Cancer_type == 'Colorectal'),],aes(Freq_G4,Percent_CD274_Expression))+ geom_point()+ 
      geom_smooth( method = 'lm')+
      theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'),
            legend.title=element_blank(),axis.title.x  = element_blank(),legend.position = 'none')+
      stat_cor(method = "pearson")
    x3
    
    topptx(x1,file=paste0(plot.path,"Breast_CD274_expression.pptx"))
    topptx(x2,file=paste0(plot.path,"Lung_CD274_expression.pptx"))
    topptx(x3,file=paste0(plot.path,"Colorectal_CD274_expression.pptx"))
  } # cd274
}
# Survival
if(F){
  library(survminer)
  library(survival)
  df.sur = read.csv('/media/yuansh/14THHD/Pan_Cancer_landscape/data/PRJNA591860_LUNG/Survival.csv',row.names = 1)
  ids = meta[which(meta$Dataset == 'PRJNA591860' ),]
  df = as.data.frame.array(table(ids$patient_id, ids$mut_type))
  df$label = ifelse((df$G3+df$G4) /(df$G1+df$G2+df$G3+df$G4) >0.7,'High_Risk','Low_Risk')
  table(df$label)
  ids = intersect(rownames(df),rownames(df.sur))
  df = df[ids,]
  df.sur = df.sur[ids,]
  df$event = df.sur[,c('best_response_status')]
  df$surdate = df.sur[,c('pfs_month')]
  df = na.omit(df)
  plot_title_name = 'Survival'
  annotate_x <- 10
  df = df[which(df$G1+df$G2+df$G3+df$G4 > 10),]
  df$event = ifelse(df$event == 'PD',1,0)
  if(T){
    survival<-df[,c('surdate','event','label')]
    survival = na.omit(survival)
    label = survival$label
    label=as.matrix(label)
    survival = as.matrix(survival[,1:2])
    summary = summary(coxph(Surv(survival[,1], survival[,2]) ~ label))
    HR = round(summary$coefficients[2],2)
    CI95 = round(summary$conf.int[3:4],2)
    LogRankP = signif(summary$sctest[3],digits = 3)
    pp = paste("LogRank p = ",LogRankP)
    HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
    label[label==1,1] <- 'High-risk'
    label[label==0,1] <- 'Low-risk'
    xlab<-"Progression-Free Interval(months)"
    ylab<-"Survival probability probability"
    surv_TTP<-survfit(Surv(survival[,1], survival[,2]) ~ label)
    res <- ggsurvplot(fit = surv_TTP,data = as.data.frame(label), 
                      legend = "bottom",legend.labs=c("High-risk","Low-risk"),legend.title="",
                      size = 2.1,fontsize=7,palette = c("#99151B","#305596"),
                      ggtheme = theme_minimal(),
                      risk.table = T,risk.table.title = "Number at risk",
                      risk.table.y.text = F,risk.table.col="strata", cumevents = F,
                      #break.time.by = deg,
                      conf.int = F,surv.median.line = 'none',combine = T,newpage=T,
                      ncensor.plot=F,xlab = xlab , ylab = ylab
    )
    res
    res$table <- res$table + 
      theme(axis.line = element_blank(),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.x = element_text(family = "serif",face ="bold.italic",size=12 )
      )
    res$plot <- res$plot +  
      labs(title=plot_title_name) + 
      theme(axis.title.x = element_text(family = "serif",face ="bold.italic",size=12),
            axis.title.y = element_text(family = "serif",face ="bold.italic",size=12),
            axis.text.x = element_text(family = "serif",face ="bold",size=12),
            axis.text.y = element_text(family = "serif",face ="bold",size=12),
            plot.title = element_text(hjust = 0.5,size=12),
            legend.text = element_text(size = 12)
            
      )
    res$plot <- res$plot + 
      annotate("text", x = annotate_x,y= 0.15,
               label = pp,size = 4,family = "serif",
               colour = "black", fontface = "bold" #,fontface = "italic"
      ) + 
      annotate("text", x = annotate_x,y= 0.08,
               label = HHRR,size = 4,family = "serif",
               colour = "black",fontface = "bold" #,fontface = "italic",
      )
    print(res)
  }
  x = res$plot
  topptx(x,file=paste0(plot.path,"Pancancer_Survival-curve.pptx"))
  
}
# Pancancer-gene-diversity-boxplot
if(F){
  df = as.data.frame.array(table(meta$Dataset,meta$mut_type))
  ids = apply(df,1, sum)
  df$Freq_G1 = df$G1 / ids
  df$Freq_G2 = df$G2 / ids
  df$Freq_G3 = df$G3 / ids
  df$Freq_G4 = df$G4 / ids
  df = df[,c("Freq_G1","Freq_G2" ,"Freq_G3", 'Freq_G4')]
  colnames(df) = c('G1','G2','G3','G4')
  df$Cell_type = rownames(df)
  df$Cell_type = factor(df$Cell_type,levels = df$Cell_type)
  df = df[,1:5]
  df = reshape2::melt(df)
  mut_type.cols = c(G1='#F8766D',G2='#ECA86E',G3='#FFEF9C',G4='#7FA4D1')
  
  x = ggplot(df,aes( value, Cell_type,fill=variable)) + 
    geom_bar(stat = 'identity') + scale_fill_manual(values =mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank(),
                                                                                   legend.position = "top")  
  x
  topptx(x,file=paste0(plot.path,"Pancancer_mut-freq-dataset.pptx"))
  
  
  
  meta$group = paste(meta$Cancer_type, meta$mut_type,sep = '-')
  x = ggplot(meta, aes(nFeature_RNA,group,fill=mut_type))+
    geom_boxplot()+facet_grid(~Cancer_type)+theme(axis.title =  element_blank()) + scale_fill_manual(values =mut_type.cols)+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank(),
                                                                                   legend.position = "top")  
  
  topptx(x,file=paste0(plot.path,"Pancancer-gene-diversity-boxplot.pptx"))
}
#TCGA survival
if(F){
  library(survminer)
  library(survival)
  files = dir('~/Desktop/TCGA/',pattern = 'BRCA|LUAD')
  for(file in files){
    ids = str_split(file,' ',simplify = T)[,5] %>% 
      gsub('\\.tsv|\\(|\\)','',.) %>% 
      gsub('-','_',.)
    assign(ids,read.table(paste0('~/Desktop/TCGA/',file),sep = '\t',header = T))
  }
  ids = str_split(files,' ',simplify = T)[,5] %>% 
    gsub('\\.tsv|\\(|\\)','',.) %>% 
    gsub('-','_',.)

  for(i in 1:length(ids)){
    df = get(ids[i])
    df = df[,c(4,3,5)]
    colnames(df) = c('surdate','event','label')
    df$surdate = df$surdate /31
    plot_title_name=ids[i]
    annotate_x=100
    if(T){
      survival<-df[,c('surdate','event','label')]
      survival = na.omit(survival)
      label = survival$label
      label=as.matrix(label)
      survival = as.matrix(survival[,1:2])
      summary = summary(coxph(Surv(survival[,1], survival[,2]) ~ label))
      HR = round(summary$coefficients[2],2)
      CI95 = round(summary$conf.int[3:4],2)
      LogRankP = signif(summary$sctest[3],digits = 3)
      pp = paste("LogRank p = ",LogRankP)
      HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
      label[label==1,1] <- 'High-risk'
      label[label==0,1] <- 'Low-risk'
      xlab<-"Progression-Free Interval(months)"
      ylab<-"Survival probability probability"
      surv_TTP<-survfit(Surv(survival[,1], survival[,2]) ~ label)
      res <- ggsurvplot(fit = surv_TTP,data = as.data.frame(label), 
                        legend = "bottom",legend.labs=c("High-risk","Low-risk"),legend.title="",
                        palette = c("#99151B","#305596"),
                        ggtheme = theme_minimal(),
                        risk.table = T,risk.table.title = "Number at risk",
                        risk.table.y.text = F,risk.table.col="strata", cumevents = F,
                        conf.int = F,surv.median.line = 'none',combine = T,newpage=T,
                        ncensor.plot=F,xlab = xlab , ylab = ylab
      )
      res$plot <- res$plot + 
        annotate("text", x = annotate_x,y= 0.15,
                 label = pp,
                 colour = "black"
        ) + 
        annotate("text", x = annotate_x,y= 0.08,
                 label = HHRR,
                 colour = "black"
        )
    }
    x = res$plot
    x
    ggsave(filename = paste0('~/Desktop/TCGA/',ids[i],'.pdf'),  width = 5,
           height = 4)
  }
}
# lung
if(F){
  df.lung = meta[which(meta$Cancer_type == 'Lung'),]
  unique(df.lung$mut_type)  
  unique(df.lung$Dataset)
  result = NULL
  for(gene in gene.list){
    freq = as.matrix(table(df.lung[['Dataset']],df.lung[[gene]]))
    ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
    df = data.frame(MUT = freq[,1],
                    WT = freq[,2],
                    Freq = ids,
                    Gene = gene)
    df$group = rownames(df)
    result = rbind(result, df)
  }
  ids1 = result[which(result$group=='PRJNA591860'),]
  ids2 = result[which(result$group=='GSE186344'),]
  
  ids1 = ids1[which(ids1$Freq<30),]
  ids2 = ids2[which(ids2$Freq>50),]
  intersect(ids1$Gene, ids2$Gene)
  
  ids = ids[order(ids$Freq),'Gene'] %>% unique()
  result$Gene = factor(result$Gene, levels = ids)
  
  ggplot(result, aes(Gene,Freq,fill=group)) + 
    geom_bar(stat='identity')+facet_wrap(~group,nrow = 2)+
    theme(axis.title =  element_blank(),
          axis.text.x =   element_text(angle=90))+
    theme(strip.background = element_rect(fill='#FFFFFF',colour='#FFFFFF'))+ theme(legend.title=element_blank(),
                                                                                   axis.title.x  = element_blank(),
                                                                                   axis.title.y  = element_blank())
  
  
}
# venn
if(F){
  df=meta
  df.lung = meta[which(meta$Cancer_type == 'Lung'),]
  df.brca = meta[which(meta$Cancer_type == 'Breast'),]
  df.coda = meta[which(meta$Cancer_type == 'Colorectal'),]
  
  result = NULL
  for(gene in gene.list){
    freq = as.matrix(table(df.coda$mut_type,df.coda[[gene]]))
    ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
    df = data.frame(WT = freq[,1],
                    MUT = freq[,2],
                    Freq = ids,
                    Gene = gene)
    df$group = rownames(df)
    result = rbind(result, df)
  }
  ids1 = result[which(result$group=='G1'),]
  ids4 = result[which(result$group=='G4'),]
  df3 = data.frame(fre.diff = ids4$Freq - ids1$Freq,
                   fre.rat = ids4$Freq / ids1$Freq,row.names = ids1$Gene)
  result = NULL
  for(gene in gene.list){
    freq = as.matrix(table(df.brca$mut_type,df.brca[[gene]]))
    ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
    df = data.frame(WT = freq[,1],
                    MUT = freq[,2],
                    Freq = ids,
                    Gene = gene)
    df$group = rownames(df)
    result = rbind(result, df)
  }
  ids1 = result[which(result$group=='G1'),]
  ids4 = result[which(result$group=='G4'),]
  df2 = data.frame(fre.diff = ids4$Freq - ids1$Freq,
                   fre.rat = ids4$Freq / ids1$Freq,row.names = ids1$Gene)
  
  result = NULL
  for(gene in gene.list){
    freq = as.matrix(table(df.lung$mut_type,df.lung[[gene]]))
    ids  = freq[,2] / (freq[,1] + freq[,2]) * 100
    df = data.frame(WT = freq[,1],
                    MUT = freq[,2],
                    Freq = ids,
                    Gene = gene)
    df$group = rownames(df)
    result = rbind(result, df)
  }
  ids1 = result[which(result$group=='G1'),]
  ids4 = result[which(result$group=='G4'),]
  df1 = data.frame(fre.diff = ids4$Freq - ids1$Freq,
                   fre.rat = ids4$Freq / ids1$Freq,row.names = ids1$Gene)
  
  
  df1$cancer = 'Lung'
  df2$cancer = 'Brca'
  df3$cancer = 'Coad'
  
  df1$genes = rownames(df1)
  df2$genes = rownames(df2)
  df3$genes = rownames(df3)
  df = rbind(df1,df2,df3)
  
  ids = df[which(df$fre.diff>50 & df$fre.rat > 30),]
  dim(ids)
  unique(ids$genes) %>% length()
  
  lung = ids[which(ids$cancer == 'Lung'),'genes']
  coad = ids[which(ids$cancer == 'Coad'),'genes']
  brca = ids[which(ids$cancer == 'Brca'),'genes']
  
  x = list(
    lung = lung, 
    coad = coad, 
    brca = brca
  )
  library('ggVennDiagram')
  x = ggVennDiagram(x, label_alpha = 0)+
    scale_fill_gradient(low="#FFFFFF",high = "#7FA4D1")
  x
  topptx(x,file=paste0(plot.path,"Pancancer_Veen.pptx"))
}
