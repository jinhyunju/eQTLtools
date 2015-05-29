#' Simulating phenotype position information
#'
#' Generates a data frame with simulated chromosome position information
#'
#' @param n.pheno Number of phenotypes to use
#' @param n.chr Number of chromosomes to simulate
#' @return Dataframe with phenotype, pheno_chr, pheno_start, pheno_end, idx as columns
#'
#'
#' @export
pheno_position_simulator <- function(n.pheno, n.chr){

  sim.chr <- sort(sample(c(1:n.chr), n.pheno, replace = TRUE))
  sim.start <- Reduce(c,lapply(table(sim.chr), function(x) 1:x))
  sim.end <- sim.start + 1000


  sim.geneinfo <- data.frame("phenotype" = paste("gene",c(1:n.pheno), sep = "_"),
                             "pheno_chr" = sim.chr,
                             "pheno_start" = sim.start,
                             "pheno_end" = sim.end,
                             "idx" = c(1:n.pheno))

  return(sim.geneinfo)

}
