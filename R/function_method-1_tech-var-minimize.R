#' Evaluate Performance of Normalization Methods by Comparing Technical Variation among Replicates
#' @description The function generates distribution plots of normalized DSP data of the ROIs that are treated as biological replicates by a list of normalization methods. A better normalization method is expected to yield smaller technical variations among biological replicates.
#' @return ggplot2 object of distribution plot.
#' @export
#' @import graphics stats ggplot2
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs to be compared, for instance, norm.methods.list = 1:17.
#' @param grouping.info a character vector describing grouping info of each ROI, must be arranged in the same order of the columns of DSP expression matrix.
#' @param main title to be added to the output figure.
#' @param title.size font size of title.
#' @param deviation.from.median if TRUE, the distribution of expressions of proteins within each ROI is evaluated relative to its median.
#' @param group.focus a character vector specifying which group in \code{"grouping.info"} is to be collected for distribution plot, the ROIs in the selected group are to be treated as biological replicates.
#' @param range.universal a numerical vector of length two, specifying the range of distribution to be plotted. If unspecified, then the range is decided automatically by the actual range of input data.
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
#' DSPNorm.eval.tech.var(normed.data.list = norm.data.list.s,
#'                       norm.methods.list = c(1:5, 8:16),
#'                       grouping.info = group.info.pick[ ,"PanCK.Status"],
#'                       main = "Sample D9",
#'                       title.size = 6.5,
#'                       deviation.from.median=FALSE,
#'                       group.focus = "PanCK+",
#'                       range.universal = c(-2.5,9.5))
DSPNorm.eval.tech.var <- function(normed.data.list, norm.methods.list, grouping.info, main="", title.size=8, deviation.from.median=FALSE, group.focus="treatment", range.universal=c(NA, NA)){

  plots.list = list()
  plot_id = 1


  # collect distribution plots for all normalization methos on the same page:
  grp.uniq = group.focus

  for (grp_id in 1:length(grp.uniq)) {

    for (method_id in norm.methods.list) {

      data.scaled.bg_corrected.normed.pt.i = normed.data.list[[method_id]]

      for.densityplot = data.scaled.bg_corrected.normed.pt.i[ ,grouping.info==grp.uniq[grp_id]]



      for.densityplot.long = data.frame(Expr.Norm=for.densityplot[ ,1])
      if(deviation.from.median){
        for.densityplot.long = data.frame(Expr.Norm=(for.densityplot[ ,1] - median(as.numeric(for.densityplot[ ,1]))))
      }
      for.densityplot.long$ROI = colnames(for.densityplot)[1]

      for (col_id in 2:ncol(for.densityplot)) {
        for.densityplot.long.temp = data.frame(Expr.Norm=for.densityplot[ ,col_id])
        if(deviation.from.median){
          for.densityplot.long.temp = data.frame(Expr.Norm=(for.densityplot[ ,col_id] - median(as.numeric(for.densityplot[ ,col_id]))))
        }
        for.densityplot.long.temp$ROI = colnames(for.densityplot)[col_id]

        for.densityplot.long = rbind(for.densityplot.long, for.densityplot.long.temp)

      }

      myplot.keep = ggplot2::ggplot(for.densityplot.long, aes(Expr.Norm, fill = ROI))
      myplot.keep = myplot.keep + ggplot2::geom_density(alpha = 0.4)
      if(!is.na(range.universal[1])){
        myplot.keep = myplot.keep + ggplot2::xlim(range.universal[1], range.universal[2])
      }
      myplot.keep = myplot.keep + ggplot2::theme(axis.text=ggplot2::element_text(size=7, face = "bold"),
                                                 axis.title=ggplot2::element_text(size=12, face = "bold"),
                                                 plot.title = ggplot2::element_text(size=title.size, face = "bold"),
                                                 legend.position = "none",
                                                 panel.grid.minor = ggplot2::element_blank(),
                                                 panel.background = ggplot2::element_blank(),
                                                 axis.line = ggplot2::element_line(colour = "black"))
      if(main==""){
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste("Method ", method_id, sep = ""))
      } else{
        myplot.keep = myplot.keep + ggplot2::ggtitle(paste(main, "\n","Method ", method_id, sep = ""))
      }

      plots.list[[plot_id]] = myplot.keep
      plot_id = plot_id + 1

    }


  }



  return(plots.list)
}


