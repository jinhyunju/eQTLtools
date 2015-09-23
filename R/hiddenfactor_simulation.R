#' Simulating hidden factor effects
#'
#' Generates a phenotype matrix based on a given input genotype matrix.
#'
#' @param n.genes Number of genes in the original phenotype matrix.
#' @param n.samples Number of samples in the original phenotype matrix
#' @param fraction.affected Fraction of genes affected by the hidden factor.
#'
#' @return hidden factor matrix with dimensions n.samples x n.genes
#' @keywords keywords
#'
#'
#' @export
hf_sim <- function(n.genes, n.samples,
                   hf.type = c("sparse", "normal", "uniform"),
                   coeff.dist = c("normal","uniform","binom"),
                   fraction.affected = 0.2,
                   binom.prob = 0.4,
                   factor.effect.size = 2){

  if( hf.type == "sparse"){
    hiddenfactors <- matrix(sparse_factor(n.genes,
                                          fraction.affected = fraction.affected,
                                          effect.mean = factor.effect.size,
                                          effect.sd = 0.5),
                            nrow = n.genes,
                            ncol = 1)
  } else if (hf.type == "normal"){
    hiddenfactors <- matrix(rnorm(n.genes, 0, factor.effect.size),
                            nrow = n.genes,
                            ncol = 1)
  } else if (hf.type == "uniform"){
    hiddenfactors <- matrix(runif(n.genes, -factor.effect.size, factor.effect.size),
                            nrow = n.genes,
                            ncol = 1)
  }


  if(coeff.dist == "normal"){
    coeff <- rnorm(n.samples, 0 , 1)
  } else if (coeff.dist == "uniform"){
    coeff <- runif(n.samples, -1, 1)
  } else {
    coeff <- rbinom(n.samples, 1, binom.prob)
  }

  hf.coeff <- matrix(coeff,
                     nrow = 1,
                     ncol = n.samples,
                     byrow = TRUE)


  hf.effect <- hiddenfactors %*% hf.coeff
  return( list("factors" = hiddenfactors, "coeff" = hf.coeff, "effect" = hf.effect))
}

#' @export
sparse_factor <- function(n.genes,
                          fraction.affected,
                          effect.mean,
                          effect.sd){
  hidden.factor <- rbinom(n.genes, 1, fraction.affected) *
    rnorm(n.genes, mean = effect.mean, sd = effect.sd) *
    sample(c(-1,1), n.genes, replace = TRUE)
  return(hidden.factor)
}

#' @export
convert2matrix <- function(vector, nr, nc){
  return(matrix(vector, nrow = nr, ncol = nc))
}
