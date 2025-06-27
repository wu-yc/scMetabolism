#' scMetabolism
#'
#' scMetabolism
#' @param obj
#' @param pathway
#' @param phenotype
#' @keywords scMetabolism
#' @examples
#' DimPlot.metabolism()
#' DotPlot.metabolism()
#' BoxPlot.metabolism()
#' @export DimPlot.metabolism
#' @export DotPlot.metabolism
#' @export BoxPlot.metabolism


library(ggplot2)
library(wesanderson)
library(data.table)

DimPlot.metabolism <- function(obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, size= 1){
  
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")

  #umap
  if (dimention.reduction.type == "umap"){

    if (dimention.reduction.run == T) obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc<-obj@reductions$umap@cell.embeddings

    row.names(umap.loc)<-colnames(obj)
    signature_exp<-obj@assays$METABOLISM$score

    input.pathway <- pathway

    signature_ggplot<-data.frame(umap.loc, t(signature_exp[input.pathway,]))

    library(wesanderson)
    pal <- wes_palette("Zissou1", 100, type = "continuous")


    library(ggplot2)
    plot <- ggplot(data=signature_ggplot, aes(x=UMAP_1, y=UMAP_2, color = signature_ggplot[,3])) +  #this plot is great
      geom_point(size = size) +
      scale_fill_gradientn(colours = pal) +
      scale_color_gradientn(colours = pal) +
      labs(color = input.pathway) +
      #xlim(0, 2)+ ylim(0, 2)+
      xlab("UMAP 1") +ylab("UMAP 2") +
      theme(aspect.ratio=1)+
      #theme_bw()
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }

  #tsne
  if (dimention.reduction.type == "tsne"){
    if (dimention.reduction.run == T) obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc<-obj@reductions$tsne@cell.embeddings

    row.names(tsne.loc)<-colnames(obj)
    signature_exp<-obj@assays$METABOLISM$score

    input.pathway <- pathway

    signature_ggplot<-data.frame(tsne.loc, t(signature_exp[input.pathway,]))

    pal <- wes_palette("Zissou1", 100, type = "continuous")


    library(ggplot2)
    plot <- ggplot(data=signature_ggplot, aes(x=tSNE_1, y=tSNE_2, color = signature_ggplot[,3])) +  #this plot is great
      geom_point(size = size) +
      scale_fill_gradientn(colours = pal) +
      scale_color_gradientn(colours = pal) +
      labs(color = input.pathway) +
      #xlim(0, 2)+ ylim(0, 2)+
      xlab("tSNE 1") +ylab("tSNE 2") +
      theme(aspect.ratio=1)+
      #theme_bw()
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))


  }
  plot
}


DotPlot.metabolism <- function(obj, pathway, phenotype, norm = "y"){
  input.norm = norm
  input.pathway <- pathway
  input.parameter<-phenotype

  metadata<-obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score

  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")


  metadata[,input.parameter]<-as.character(metadata[,input.parameter])
  metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])

  #arrange large table
  gg_table<-c()
  for (i in 1:length(input.pathway)){
    gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
  }
  gg_table<-data.frame(gg_table)

  #get median value
  gg_table_median<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))


  for (x in 1:length(input.group.x)){
    for (y in 1:length(input.group.y)){
      gg_table_sub<-subset(gg_table, gg_table[,1] == input.group.x[x] & gg_table[,2] == input.group.y[y])
      gg_table_median<-rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], median(as.numeric(as.character(gg_table_sub[,3])))))

    }
  }
  gg_table_median<-data.frame(gg_table_median)
  gg_table_median[,3]<-as.numeric(as.character(gg_table_median[,3]))


  #normalize
  gg_table_median_norm<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))


  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  if (input.norm == "y")
    for (y in 1:length(input.group.y)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,2] == input.group.y[y])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }

  if (input.norm == "x")
    for (x in 1:length(input.group.x)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,1] == input.group.x[x])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }

  if (input.norm == "na") gg_table_median_norm<-gg_table_median


  gg_table_median_norm<-data.frame(gg_table_median_norm)
  gg_table_median_norm[,3]<-as.numeric(as.character(gg_table_median_norm[,3]))




  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")

  ggplot(data=gg_table_median_norm, aes(x=gg_table_median_norm[,1], y=gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
    geom_point(data=gg_table_median_norm, aes(size = gg_table_median_norm[,3])) + #geom_line() +
    #theme_bw()+theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Metabolic Pathway")+ xlab(input.parameter)+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    scale_color_gradientn(colours = pal) +
    labs(color = "Value", size = "Value") +
    #facet_wrap(~tissueunique, ncol = 1) +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
}


BoxPlot.metabolism <- function(obj, pathway, phenotype, ncol = 1){
  input.pathway<-pathway
  input.parameter<-phenotype

  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")

  metadata<-obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score



  metadata[,input.parameter]<-as.character(metadata[,input.parameter])
  metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])

  #arrange large table
  gg_table<-c()
  for (i in 1:length(input.pathway)){
    gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
  }
  gg_table<-data.frame(gg_table)
  gg_table[,3]<-as.numeric(as.character(gg_table[,3]))


  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")

  ggplot(data=gg_table, aes(x=gg_table[,1], y=gg_table[,3], fill = gg_table[,1])) +
    geom_boxplot(outlier.shape=NA)+
    ylab("Metabolic Pathway")+
    xlab(input.parameter)+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    #scale_color_gradientn(colours = pal) +
    facet_wrap(~gg_table[,2], ncol = ncol, scales = "free") +
    labs(fill = input.parameter) +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
}



