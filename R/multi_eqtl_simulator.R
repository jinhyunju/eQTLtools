#' @export
multi_eqtl_simulator <- function(input.geno = NULL,
                                 n.pheno,
                                 n.eqtl,
                                 simulation.id,
                                 working.directory,
                                 coeff.ratio = 2.5, # now it is mean of normal distribution
                                 n.sim = 5,
                                 ct.ratio = 0.7,
                                 trans.impact = 80,
                                 trans.nerf = 0.6,
                                 hidden.factor = TRUE,
                                 factors = NULL,
                                 factor.coeff = NULL,
                                 hf.frac = 0.2,
                                 effect.size = 2){

  n.sample <- nrow(input.geno)
#  coeff.dist <- rexp(100000,coeff.ratio) * sample(c(-1,1), 100000, replace = TRUE)
#  coeff.sample <- coeff.dist[which(abs(coeff.dist) > 1)]
  sim.geneinfo <- pheno_position_simulator(n.pheno, 15)

  # Genotype LD estimation (using simple correlation)
  geno.cor.mx <- cor(input.geno)
  sparse.cor.mx <- 1*(abs(geno.cor.mx) > 0.8)

  output.path <- working.directory

  dir.create(output.path, showWarnings = FALSE)

  for (k in 1:n.sim){
    ################################################################################
    #################              Simulate Phenotype              #################
    ################################################################################
    pheno.list <- list()

    eqtl.sim <- eqtl_simulator(input.geno, n.pheno = n.pheno,
                               n.eqtl = n.eqtl,
                               cis.trans.ratio = ct.ratio,
                               trans.impact = trans.impact,
                               trans.nerf = trans.nerf,
                               coeff.sample = coeff.ratio,
                               cis.sd = 0.5)


    pheno.list[["noisefree"]] <- eqtl.sim$phenotype


    effect.mx <- eqtl.sim$effectmx

    sparse.effect <- 1*(effect.mx != 0)


    # generating ground truth that accounts for LD structure (based on correlation)
    # loop over phenotype
    for( p in 1: ncol(sparse.effect)){
      # identify genotypes which are contributing
      geno.pos <- which(sparse.effect[,p] == 1)
      for (single.geno in geno.pos){
        ld.block <- which(sparse.cor.mx[single.geno,] == 1)
        sparse.effect[ld.block, p] <- 1
      }
    }

    ground.truth <- data.frame(which(sparse.effect == 1, arr.ind = TRUE))
    ground.truth$pair <- paste(ground.truth[,1], ground.truth[,2], sep = "_")

    generative.truth <- eqtl.sim$eqtl

    colnames(pheno.list[["noisefree"]]) <- sim.geneinfo$phenotype


    ################################################################################
    #################              Add Gaussian Noise              #################
    ################################################################################

    pheno.list[["noise"]] <- pheno.list[["noisefree"]] +
                                        rnorm(length(pheno.list[["noisefree"]]),
                                              mean = 0,
                                              sd = 0.5)


    ################################################################################
    #################              Add Hidden factors              #################
    ################################################################################
    if(hidden.factor){
      factor.details <- matrix(cbind(factors,factor.coeff),
                               ncol = 2,
                               byrow = FALSE)

      hf_list <- apply(factor.details, 1, function(x) hf_sim(n.genes = n.pheno,
                                                             n.samples = n.sample,
                                                             hf.type = x[1],
                                                             coeff.dist = x[2],
                                                             fraction.affected = hf.frac,
                                                             factor.effect.size = effect.size))

      # add all hidden factor effects
      hf.effect <- Reduce(`+`, lapply(hf_list, function(x) x$effect))


      pheno.list[["hf"]] <- pheno.list[["noise"]] + t(hf.effect)
      save(pheno.list, hf_list, sim.geneinfo, factor.details, ground.truth, input.geno,
           generative.truth,
           file = paste(output.path, "/",simulation.id, "_pheno", k, ".RData", sep = ""))
    } else {

      pheno.list[["hf"]] <- pheno.list[["noise"]]
      save(pheno.list, sim.geneinfo, ground.truth, input.geno, generative.truth,
           file = paste(output.path, "/",simulation.id, "_pheno", k, ".RData", sep = ""))

    }


  }
}
