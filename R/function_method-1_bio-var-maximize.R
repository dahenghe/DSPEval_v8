#' Evaluate Performance of Normalization Methods by Comparing Biological Variation
#' @description The function generates PCA plots of normalized DSP data of the ROIs that are treated as biological replicates by a list of normalization methods. A better normalization method is expected to yield larger biological variations between two known biological subgroups.
#' @return ggplot2 object of PCA plot.
#' @export
#' @import graphics stats ggplot2
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs of the list object of normed.data.list to be compared, for instance, norm.methods.list = 1:17.
#' @param grouping.info a character vector describing grouping info of each ROI, must be arranged in the same order of the columns of DSP expression matrix. Each unique character in the vector is to be treated as a unique biological subgroup. For optimal visualization, it is strongly recommended to consider only two unique biological groups of interest in each round of comparison.
#' @param main title to be added to the output figure.
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
#' DSPNorm.eval.bio.var(normed.data.list = norm.data.list.s,
#'                      norm.methods.list = c(1:5, 8:16),
#'                      grouping.info = group.info.pick[ ,"PanCK.Status"],
#'                      main = "Sample D9",
#'                      title.size = 6.5)
DSPNorm.eval.bio.var <- function(normed.data.list, norm.methods.list, grouping.info, main="", title.size=7){



  plots.list = list()
  plot_id = 1


  # collect PCA plots for all normalization methods on the same page:
  if(length(unique(grouping.info)) > 1){
    for (method_id in norm.methods.list) {

      data.scaled.bg_corrected.normed.pt.i = normed.data.list[[method_id]]


      for.pca = as.data.frame(t(data.scaled.bg_corrected.normed.pt.i))

      x <- prcomp(for.pca)
      z1 <- data.frame(State = rownames(x$x), x$x[, 1:2])
      x.summary = summary(x)

      myplot.keep = ggplot2::ggplot(data = z1, aes(x=PC1, y=PC2, color = grouping.info), show.legend = FALSE) +
        ggplot2::xlab(paste("PC1 (", signif(x.summary$importance["Proportion of Variance","PC1"], 3)*100,"%)", sep = "")) +
        ggplot2::ylab(paste("PC2 (", signif(x.summary$importance["Proportion of Variance","PC2"], 3)*100,"%)", sep = "")) +
        ggplot2::geom_point(shape=1) +
        ggplot2::stat_ellipse(geom = "polygon",
                              level = 0.95,
                              aes(fill = grouping.info),
                              alpha = 0.25) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       legend.position = "none",
                       axis.line = ggplot2::element_line(colour = "black"),
                       plot.title = ggplot2::element_text(size=title.size, face = "bold"),
                       axis.title=ggplot2::element_text(size=14,face="bold"),
                       axis.text=ggplot2::element_text(size=9,face="bold")
        )

      if(main==""){
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste("Method ", method_id, sep = ""))
      } else{
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste(main, "\n", "Method ", method_id, sep = ""))
      }


      plots.list[[plot_id]] = myplot.keep


      plot_id = plot_id + 1

    }

  } else{
    cat("   Found only one group, PCA plot is therefore ignored!\n")
  }


  return(plots.list)
}











