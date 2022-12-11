#' For internal usage of box plot
#' @description For internal usage of box plot
#' @import graphics stats
#' @param y vector of values to plot
#' @param x vector of grouping info
#' @param title title to show
#' @param type groups of interest
#' @param ylab ylab to show
#' @noRd
myboxplot.DSP <- function(y, x, title="", type=c("Group A", "Group B"), ylab="Protein Expr")
{

  y = y[!is.na(x)]; x = x[!is.na(x)]
  gr = unique(x)[order(unique(x),decreasing = F)]
  m = length(unique(x))
  par(mar=c(4.5,4.5,3.5,1.1))
  pval.wilcox = wilcox.test(y[x==gr[1]], y[x==gr[2]], exact = F)$p.value



  main.use <- paste("Wilcoxon P = ", signif(pval.wilcox[1],3), sep="")
  boxplot(y ~ as.factor(x), col=c("white", "white"), xlab = "", xaxt="n", range=0, outline=F, main = main.use, ylab=ylab)

  for (i in 1:m){
    vv = y[x==gr[i]]
    # add some random noise to x position of each points
    points(jitter(rep(i,length(vv)),amount=0.1), vv, pch=15, col=c("blue", "orange", "brown", "purple","red","green","black","cyan","grey")[i], cex=0.8)
  }
  axis(gr, side = 1, at = 1:m, padj=0.5, las = 1)



  box(lwd=2)


}

