#' Evaluate Performance of Normalization Methods by Comparing Spearman correlation between DSP protein expression and the corresponding mRNA-level expression
#' @description The function generates box plots of Spearman correlations of matched protein-mRNA expressions. A better normalization method is expected to yield stronger correlation between protein expression and the corresponding mRNA expression within the same ROI.
#' @return ggplot2 object of box plot, for both Spearman correlation strengths and Spearman correlation p values.
#' @export
#' @import graphics stats ggplot2
#' @param normed.data.list a list object containing the normalized DSP expression matrix under all normalization methods, each normalized matrix can be either a data.frame or matrix.
#' @param norm.methods.list a numerical vector of the method IDs to be compared, for instance, norm.methods.list = 1:17.
#' @param reference.data.matrix mRNA-Seq expression data of either a matrix or a data.frame class. Must be arranged in exactly the same order of the DSP protein expression data.
#' @param main title to be added to the output figure.
#' @examples
#' data(COVID.19.DSP)
#'
#' high.Expr.normed.DSP.matched.to.WTA <-
#'    COVID.19.DSP$high.Expr.normed.DSP.matched.to.WTA
#'
#' high.Expr.Q3normed.WTA <-
#'    COVID.19.DSP$high.Expr.Q3normed.WTA
#'
#' DSPNorm.eval.spearman.correlation(
#'    normed.data.list = high.Expr.normed.DSP.matched.to.WTA,
#'    norm.methods.list = c(1:5, 8:16),
#'    reference.data.matrix = high.Expr.Q3normed.WTA,
#'    main = "All Patients: across ROIs within Each Molecule")
DSPNorm.eval.spearman.correlation <- function(normed.data.list, norm.methods.list, reference.data.matrix, main=""){


  tab.out.spearman.correlation = matrix(NA, nrow = 1, ncol = 3)

  for (method_id in norm.methods.list) {

    normed.data.i = normed.data.list[[method_id]]

    if(is.list(reference.data.matrix)){
      reference.data.i = reference.data.matrix[[method_id]]
    } else{
      reference.data.i = reference.data.matrix
    }

    for (jj in 1:ncol(normed.data.i)) {

      normed.data.i.j = as.numeric(normed.data.i[ ,jj])

      reference.data.i.j = as.numeric(reference.data.i[ ,jj])

      if(sd(normed.data.i.j, na.rm = T)!=0 & sd(reference.data.i.j, na.rm = T)!=0){

        vv = stats::cor.test(normed.data.i.j, reference.data.i.j, method = "spearman", exact = F)
        #vv = cor.test(normed.data.i.j, reference.data.i.j, method = "pearson", exact = F)


        estimate.spearman = stats::cor(normed.data.i.j, reference.data.i.j, method = "spearman")
        #estimate.spearman = cor(normed.data.i.j, reference.data.i.j, method = "pearson")

        tab.out.spearman.correlation.temp = c(method_id, estimate.spearman, vv$p.value)

        tab.out.spearman.correlation = rbind(tab.out.spearman.correlation, tab.out.spearman.correlation.temp)

      } else{

        tab.out.spearman.correlation.temp = c(method_id, NA, NA)

        tab.out.spearman.correlation = rbind(tab.out.spearman.correlation, tab.out.spearman.correlation.temp)
      }



    }

  }

  colnames(tab.out.spearman.correlation) = c("Normalization.Method", "Spearman.Correlation", "Spearman.PValue")

  tab.out.spearman.correlation = as.data.frame(tab.out.spearman.correlation[-1, ])

  tab.out.spearman.correlation = tab.out.spearman.correlation[!is.na(tab.out.spearman.correlation$Spearman.PValue), ]


  tab.out.spearman.correlation$Spearman.Correlation = as.numeric(tab.out.spearman.correlation$Spearman.Correlation)

  tab.out.spearman.correlation$Spearman.PValue = as.numeric(tab.out.spearman.correlation$Spearman.PValue)

  tab.out.spearman.correlation$Negative.Log.Spearman.PValue = -log(tab.out.spearman.correlation$Spearman.PValue)

  # sort methods based in their median p-values:
  tab.out.spearman.correlation$Negative.Log.Spearman.PValue.Median.within.Method <- NA
  for (method_id in c(1:5, 8:16)) {

    tab.out.spearman.correlation$Negative.Log.Spearman.PValue.Median.within.Method[tab.out.spearman.correlation$Normalization.Method==method_id] <- median(tab.out.spearman.correlation$Negative.Log.Spearman.PValue[tab.out.spearman.correlation$Normalization.Method==method_id], na.rm = T)

  }

  tab.out.spearman.correlation <- tab.out.spearman.correlation[order(tab.out.spearman.correlation$Negative.Log.Spearman.PValue.Median.within.Method, decreasing = T), ]


  tab.out.spearman.correlation$Normalization.Method = factor(tab.out.spearman.correlation$Normalization.Method, levels = unique(tab.out.spearman.correlation$Normalization.Method))


  plots.list = list()


  #library(ggplot2)

  myplot.keep = ggplot2::ggplot(tab.out.spearman.correlation, aes(x=Normalization.Method, y=Negative.Log.Spearman.PValue)) +
    ggplot2::ggtitle(paste(main, sep = "") ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_hline(yintercept=-log(0.05), linetype="solid",
               color = "red", size=0.85)


  plots.list[[1]] = myplot.keep


  myplot.keep = ggplot2::ggplot(tab.out.spearman.correlation, aes(x=Normalization.Method, y=Spearman.Correlation)) +
    ggplot2::ggtitle(paste(main, sep = "") ) +
    ggplot2::geom_boxplot()

  plots.list[[2]] = myplot.keep


  return(plots.list)


}









