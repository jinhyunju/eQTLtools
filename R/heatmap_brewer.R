#' @import fields
#' @import parallel
#' @import RColorBrewer
#' @import parallel
#' @import zoom.plot
#'
#' @export
heatmap_brewer <- function(pvals, plot.title, n.cores = 2, sig.only = FALSE, file = TRUE){
  # saving graphics environment to restore later

  default.par <- par(no.readonly = T)

  n.genotypes <- nrow(pvals)
  n.probes <- ncol(pvals)

  # dimensions for the png are set
  if(n.probes > 1200){
    heatmap.probes <- 1200
  } else {
    heatmap.probes <- n.probes
  }

  if(n.genotypes > 1920){
    heatmap.snps <- 1920
  } else {
    heatmap.snps <- n.genotypes
  }


  if(heatmap.snps == n.genotypes & heatmap.probes == n.probes){
    #cat("Skipping derez \n")
    heatmap <- pvals
    snp.idx <- 1:n.genotypes
    probe.idx <- 1:n.probes
  } else {
    derez.list <- zoom.plot::derez(pvals, heatmap.snps, heatmap.probes, 1:n.genotypes, 1:n.probes, mc.cores = n.cores)
    heatmap <- derez.list$x
    snp.idx <- derez.list$row.pos
    probe.idx <- derez.list$col.pos
  }

  genotype.pval.means <- rowMeans(-log10(pvals), na.rm = TRUE)
  phenotype.pval.means <- colMeans(-log10(pvals), na.rm = TRUE)

  pval.hits <- data.frame(which(pvals < (0.05 / length(pvals)), arr.ind=TRUE), row.names = NULL)

  # Approximate dimensions so that each entry in heatmap array ~ 1 pixel
  if(file){
    png(sprintf('%s_heatmap.png', plot.title), width=200+heatmap.snps, height=130+heatmap.probes, pointsize = 15)
  }


  par(mar=c(0, 0, 0, 0), oma=c(6, 6, 3, 1), pch='.', las=1)
  layout(matrix(c(2, 1, 4, 3), 2), heights=c(1, 8), widths=c(8, 1))
  par(cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  if(sig.only){
    image(snp.idx, probe.idx, matrix(0, nrow = nrow(heatmap), ncol = ncol(heatmap)),
          zlim = 0:1,
          xlab='SNPs', ylab='Genes',
          col = brewer.pal(9,"Greens"), xpd = NA, mgp = c(4,0.5,0) )
    with(pval.hits, points(row, col, cex = 3)) # have to fix this part so that it aligns with 

  } else{
    image(snp.idx, probe.idx, heatmap,
          zlim = 0:1,
          xlab='SNPs', ylab='Genes',
          col = brewer.pal(10,"RdYlBu"), xpd = NA, mgp = c(4,0.5,0) )
    with(pval.hits, points(row, col, cex = 3))
  }
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
  plot(genotype.pval.means, type='l', xaxt='n', xaxs='i')
  mtext(plot.title, side = 3, line = 0.5, outer = TRUE, cex = 2)
  plot(phenotype.pval.means, 1:ncol(pvals), type='l', yaxt='n',  yaxs='i')
  plot(0, ann=FALSE, bty='n', type='n', xaxt='n', yaxt='n')
  
  if(file){
    draw.scale(brewer.pal(10, "RdYlBu"), c(0,1), size = 3, cex = 1.5, width.to.height =10)
    dev.off()
  } else {
    draw.scale(brewer.pal(10, "RdYlBu"), c(0,1), size = 1, cex = 1)
    par(default.par)
  }

}
