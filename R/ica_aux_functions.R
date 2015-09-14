#' @import lrgpr
#' @import formula.tools
#' @export
ica_genotype_test <- function(ica.result, genotype.mx, n.cores = 1){
  ica.loadings <- t(ica.result$A)

  ic.vs.geno <- glmApply(ica.loadings ~ SNP,
                         features = genotype.mx,
                         nthreads = n.cores)$pValues

  colnames(ic.vs.geno) <- rownames(ica.result$A)
  sig <- which(ic.vs.geno < (0.05/length(ic.vs.geno) ), arr.ind = TRUE)

  genetic.factors <- colnames(ic.vs.geno)[unique(sig[,"col"])]
  non.genetic <- colnames(ic.vs.geno)[which(!(colnames(ic.vs.geno) %in% genetic.factors))]


  return(list("genetic" = genetic.factors, "hf" = non.genetic))
}

#' @import lrgpr
#' @import formula.tools
#' @export
pca_genotype_test <- function(pca.result, genotype.mx, n.cores = 1){
  pca.loadings <- pca.result$x

  pc.vs.geno <- glmApply(pca.loadings ~ SNP,
                         features = genotype.mx,
                         nthreads = n.cores)$pValues

  colnames(pc.vs.geno) <- rownames(pca.result$x)
  sig <- which(pc.vs.geno < (0.05/length(pc.vs.geno) ), arr.ind = TRUE)

  genetic.factors <- colnames(pc.vs.geno)[unique(sig[,"col"])]
  non.genetic <- colnames(pc.vs.geno)[which(!(colnames(pc.vs.geno) %in% genetic.factors))]


  return(list("genetic" = genetic.factors, "hf" = non.genetic))
}



#' @import lrgpr
#' @export
get_similarity_mx <- function(ica.result, hidden.factors){
  lmm.mx <- ica.result$A[hidden.factors,]
  lmm.norm <- apply(lmm.mx, 1, function(x) ((x - mean(x)) / sd(x)))


  weights <- ica.result$ica.stat.df[colnames(lmm.norm),"percent.var"] / sum(ica.result$ica.stat.df[colnames(lmm.norm),"percent.var"])

  weighted.lmm <- apply(lmm.norm, 1, function(x) x * weights)

  similarity.mx <- t(weighted.lmm) %*% weighted.lmm
  return(similarity.mx)
}
