#' Evaluate Performance of Normalization Methods by Comparing Coefficients of Variation
#' @description The function generates box plots of coefficients of variation (COV) of protein across the ROIs of biological replicates. A better normalization is expected to yield a lower level of COV.
#' @return ggplot2 object of box plot.
#' @export
#' @import graphics stats ggplot2
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs to be compared, for instance, norm.methods.list = 1:17.
#' @param controls a character vector listing both housekeeping and negative control IDs, which are to be ignored in evaluation of COV.
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
#' pick.id.2 <- group.info$PanCK.Status=="PanCK+"
#' pick.id.3 <- group.info$X.Sample_characteristics_ch1==
#'                 "infection: SARS-CoV-2 infected"
#'
#' norm.data.list.s <- list()
#' for (method_id in c(1:5, 8:16)) {
#'   data.pick <- data.normed.list[[method_id]]
#'   norm.data.list.s[[method_id]] <-
#'      data.pick[ ,pick.id.1 & pick.id.2 & pick.id.3]
#' }
#'
#' DSPNorm.eval.cov(normed.data.list = norm.data.list.s,
#'                  norm.methods.list = c(1:5, 8:16),
#'                  controls = controls,
#'                  main = "All Patients/Alveolar/Infected/PanCK+",
#'                  title.size=9)
DSPNorm.eval.cov <- function(normed.data.list, norm.methods.list, controls, main="", title.size=7){


  tab.out.cov = matrix(NA, nrow = 1, ncol = 3)

  for (method_id in norm.methods.list) {

    normed.data.i = normed.data.list[[method_id]]
    normed.data.i <- normed.data.i[!rownames(normed.data.i)%in%controls, ]


    # ii index refers to each protein:
    for (ii in 1:nrow(normed.data.i)) {

      values.ii <- as.numeric(normed.data.i[ii, ])

      if(mean(values.ii, na.rm = T)!=0 & sd(values.ii, na.rm = T)!=0){

        cov.ii <- sd(values.ii, na.rm = T)/mean(values.ii, na.rm = T)

        tab.out.cov.temp <- c(method_id, rownames(normed.data.i)[ii], cov.ii)

        tab.out.cov <- rbind(tab.out.cov, tab.out.cov.temp)
      }

    }


  }

  colnames(tab.out.cov) = c("Method", "Antibody.Name", "COV")

  tab.out.cov = as.data.frame(tab.out.cov[-1, ])



  tab.out.cov$COV = as.numeric(tab.out.cov$COV)




  tab.out.cov$Method = factor(tab.out.cov$Method, levels = unique(tab.out.cov$Method))


  plots.list = list()


  myplot.keep = ggplot2::ggplot(tab.out.cov, aes(x=Method, y=COV, fill=Method)) +
    ggplot2::ggtitle(paste(main, sep = "") ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=title.size, face = "bold"),
                   axis.title=ggplot2::element_text(size=20,face="bold"),
                   axis.text=ggplot2::element_text(size=16,face="bold"),
                   legend.position = "none"
    ) +
    ggplot2::geom_hline(yintercept=0.5, linetype="solid",
                        color = "red", size=0.85)


  plots.list[[1]] = myplot.keep



  return(plots.list)


}

