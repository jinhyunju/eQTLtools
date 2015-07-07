#' Simulating eQTL phenotypes.
#'
#' Generates a phenotype matrix based on a given input genotype matrix.
#'
#' @param genotypes A genotype matrix with dimensions N x g.
#' @param n.pheno Number of phenotypes to simulate
#' @param n.eqtl Number of eQTL pairs
#' @param cis.trans.ratio Ratio between cis and trans hits
#' @param trans.impact Maximum number of genes a trans hit can influence.
#' @param trans.nerf Ratio of trans effect size compared to cis-hits.
#'        Accounts for the observation that trans hits usually have smaller
#'        effect sizes.
#' @param coeff.sample A vector from which the effect sizes (coefficiens) are
#'        going to be sampled from.
#'
#' @return Phenotype matrix with dimensions N x n.pheno
#' @keywords keywords
#'
#'
#' @export
eqtl_simulator <- function(genotypes,
                           n.pheno = 1000,
                           n.eqtl = 800,
                           cis.trans.ratio = 0.8,
                           trans.impact = 50,
                           trans.nerf = 0.6,
                           coeff.sample = NULL,
                           cis.sd = 0.5){

  if(is.null(coeff.sample)){
    coeff.sample <- sample(c(-1,1),1) * rnorm(n.cis, mean = 3, sd = 1)
  }
  # number of phenotypes

  # number of genotypes
  n.geno <- ncol(genotypes)

  # create cis effect indicator matrix
  eqtl.effect <- matrix(0, nrow = n.geno, ncol = n.pheno)


  # creat genotype phenotype pairs for eqtl
  # sampling genotypes for eqtls based on the desired numbers of eqtl
  eqtl.geno <- sample(1:n.geno, n.eqtl)

  # number of cis relationships determined by the cis.trans.ratio
  n.cis <- round(n.eqtl * cis.trans.ratio)

  # The rest of the genotypes will be trans-eqtls
  n.trans <- n.eqtl - n.cis

  # from eqtl genotypes sample cis candidates
  cis.geno <- sample(eqtl.geno, n.cis, replace = FALSE)

  # the rest will be used for trans candidates
  trans.geno <- eqtl.geno[!(eqtl.geno %in% cis.geno)]

  # create ground truth dataframe which has all the genotype / phenotype
  # relationships and effectsize saved
  eqtl.indexes <- data.frame("geno" = sort(cis.geno),
                             "pheno" = sort(sample(1:n.pheno, n.cis, replace = FALSE)),
                             "effect" = sample(coeff.sample, n.cis, replace = FALSE),
                             "label" = "cis")

  # create a matrix that will decide how many genes a single trans genotype
  # will affect
  trans.numbers <- matrix(0, nrow = 2, ncol = n.trans)
  trans.numbers[1,] <- trans.geno
  trans.numbers[2,] <- sample(1:trans.impact, n.trans, replace = TRUE)

  # loop over trans.numbers matrix to create trans relationships
  # and add it to the eqlt.indexes data frame

  for(i in 1:ncol(trans.numbers)){
    trans.size <- trans.numbers[2,i]
    trans.pheno <- sample(1:n.pheno, trans.size, replace = FALSE)
    # common perception is that trans eqtls have smaller effet sizes
    # so here we are nerfing the effect size
    trans.indexes <- data.frame("geno" = rep(trans.numbers[1,i], trans.size),
                                "pheno" = trans.pheno,
                                "effect" = sample(coeff.sample, trans.size) * trans.nerf,
                                "label" = "trans")
    eqtl.indexes <- rbind(eqtl.indexes, trans.indexes)
  }

  for(i in 1:nrow(eqtl.indexes)){
    eqtl.effect[eqtl.indexes[i,"geno"], eqtl.indexes[i,"pheno"]] <- eqtl.indexes[i,"effect"]
  }

  phenotype.mx <- genotypes %*% eqtl.effect


  return(list("phenotype" = phenotype.mx, "eqtl" = eqtl.indexes, "effectmx" = eqtl.effect))

}
