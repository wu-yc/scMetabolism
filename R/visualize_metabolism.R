#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' plot.metabolism()
#' @export plot.metabolism

plot.metabolism <- function(obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, size= 1){
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
