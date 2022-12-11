#' Evaluate Performance of Normalization Methods by Comparing clustering based on Euclidean distance.
#' @description The function generates Euclidean-based clustering plots based on normalized DSP data. A better normalization is expected to yield a clustering prediction consistent with known biological subgroups.
#' @return a pheatmap object.
#' @export
#' @import pheatmap
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs to be compared, for instance, norm.methods.list = 1:17.
#' @param grouping.info a character vector describing grouping info of each ROI, must be arranged in the same order of the columns of DSP expression matrix. Each unique character in the vector is to be treated as a unique biological subgroup. For optimal visualization, it is strongly recommended to consider only two unique biological groups of interest in each round of comparison.
#' @param main a character string. If unspecified, then the method numbers are displayed as titles. If specified, then the specified string will be added in front of the method numbers in title.
#' @param title.size font size of title.
#' @examples
#' data(COVID.19.DSP)
#'
#' data.normed.list <- COVID.19.DSP$Expr.normed.DSP
#'
#' group.info <- COVID.19.DSP$Meta.DSP
#'
#' controls <- c("Ms IgG1", "Ms IgG2a", "Rb IgG",
#'               "S6", "GAPDH", "Histone H3")
#'
#' group.info$PanCK.Status[group.info$PanCK.Status=="PanCK-"] <- "PanCK Neg"
#' group.info$PanCK.Status[group.info$PanCK.Status=="PanCK+"] <- "PanCK Pos"
#'
#' pick.id.1 <- group.info$organ.short=="Alveolar"
#'
#' pick.id.2 <- group.info$patient=="D9"
#'
#' group.info.pick <- group.info[pick.id.1 & pick.id.2, ]
#'
#' norm.data.list.s <- list()
#' for (method_id in c(1:5, 8:16)) {
#'   data.pick <- log(data.normed.list[[method_id]])
#'   norm.data.list.s[[method_id]] <-
#'      data.pick[!row.names(data.pick)%in%controls, pick.id.1 & pick.id.2]
#' }
#'
#' DSPNorm.eval.euclidean.cluster(normed.data.list = norm.data.list.s,
#'                                norm.methods.list = c(1:5, 8:16),
#'                                grouping.info = group.info.pick$PanCK.Status,
#'                                main = "Sample D9",
#'                                title.size = 7.5)
DSPNorm.eval.euclidean.cluster <- function(normed.data.list, norm.methods.list, grouping.info, main="", title.size=7){

  plots.list = list()
  plot_id = 1

  if(length(unique(grouping.info)) > 1){

    grp.uniq = sort(unique(grouping.info))

    for (method_id in norm.methods.list) {

      data.scaled.bg_corrected.normed.pt.i = normed.data.list[[method_id]]

      annotation <- data.frame(Type=rep(NA, ncol(data.scaled.bg_corrected.normed.pt.i)))
      rownames(annotation) <- colnames(data.scaled.bg_corrected.normed.pt.i)

      for (type_id in 1:length(grp.uniq)) {

        annotation$Type[grouping.info==grp.uniq[type_id]] <- grp.uniq[type_id]

      }


      if(main==""){
        main.use <- paste("Method ", method_id, sep = "")
      } else{
        main.use <- paste(main, "\n", "Method ", method_id, sep = "")
      }
      myplot = pheatmap::pheatmap(data.scaled.bg_corrected.normed.pt.i, # the data table
                                  annotation_col = annotation, # annotation for the samples
                                  legend = F,
                                  annotation_names_col = F,
                                  annotation_legend = F,
                                  fontsize = title.size,
                                  scale = "none", # scale by row
                                  cluster_cols=T, # whether to cluster by columns
                                  cluster_rows = T, # whether to cluster by rows
                                  main = main.use,
                                  cellheight = 0,
                                  labels_row = rep("", nrow(data.scaled.bg_corrected.normed.pt.i)),
                                  labels_col = rep("", ncol(data.scaled.bg_corrected.normed.pt.i)),
                                  treeheight_row = 0)

      plots.list[[plot_id]] = myplot[[4]]

      plot_id = plot_id + 1


    }


    return(plots.list)

  } else{

    cat("   Found only one group, euclidean cluster plot is therefore ignored!\n")

    return(NULL)
  }


}


