#' Simulating eQTL phenotypes.
#'
#' Generates a phenotype matrix based on a given input genotype matrix.
#'
#' @param genotypes A genotype matrix with dimensions N x g.
#' @param n.pheno Number of phenotypes to simulate
#' @param n.eqtl Number of eQTL pairs
#' @param cis.trans.ratio Ratio between cis and trans hits
#' @param trans.impact Maximum percentage of total genes a trans hit can influence.
#' @param trans.nerf Ratio of trans effect size compared to cis-hits.
#'        Accounts for the observation that trans hits usually have smaller
#'        effect sizes.
#' @param coeff.sample A vector from which the effect sizes (coefficiens) are
#'        going to be sampled from.
#'
#' @return Phenotype matrix with dimensions N x n.pheno
#' @keywords keywords
#'
#' @import rhdf5
#' @export
eqtl_simulator_h5 <- function(input.geno = NULL,
                             n.pheno,
                             n.eqtl,
                             simulation.id,
                             output.path = "./simpheno",
                             coeff.mean = 2.5, # now it is mean of normal distribution
                             n.sim = 5,
                             cis.trans.ratio = 0.7,
                             trans.impact = 0.1,
                             trans.nerf = 0.7,
                             hidden.factor = TRUE,
                             factors = NULL,
                             factor.coeff = NULL,
                             hf.frac = 0.2,
                             effect.size = 2){
    output_h5 <- paste(gsub("-", "", Sys.Date()),simulation.id,".h5",sep = "")
    output_h5 <- paste(output.path, output_h5, sep = "")
    pheno.list <- list()
    n.sample <- nrow(input.geno)

    sim.geneinfo <- pheno_position_simulator(n.pheno, 15)

    # Genotype LD estimation (using simple correlation)
    geno.cor.mx <- cor(input.geno)
    sparse.cor.mx <- 1*(abs(geno.cor.mx) > 0.8)

    dir.create(output.path, showWarnings = FALSE)

    ################################################################################
    #################              Simulate Phenotype              #################
    ################################################################################
    #pheno.list <- list()


    n.geno <- ncol(input.geno)

    trans.number <- round(n.pheno * trans.impact)
    # create cis effect indicator matrix
    eqtl.effect <- matrix(0, nrow = n.geno, ncol = n.pheno)


    # creat genotype phenotype pairs for eqtl
    # sampling genotypes for eqtls based on the desired numbers of eqtl
    eqtl.geno <- sample(1:n.geno, n.eqtl)

    # number of cis relationships determined by the cis.trans.ratio
    n.cis <- round(n.eqtl * cis.trans.ratio)

    # sample effect of cis eQTL
    cis.coeff <- sample(c(-1,1),1) * rnorm(n.cis, mean = coeff.mean, sd = 1)

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
                             "effect" = cis.coeff,
                             "label" = "cis")

    # create a matrix that will decide how many genes a single trans genotype
    # will affect
    trans.mx <- matrix(0, nrow = 2, ncol = n.trans)
    trans.mx[1,] <- trans.geno
    trans.mx[2,] <- sample(1:trans.number, n.trans, replace = TRUE)

    # loop over trans.mx matrix to create trans relationships
    # and add it to the eqlt.indexes data frame

    for(i in 1:ncol(trans.mx)){
        trans.size <- trans.mx[2,i]
        trans.pheno <- sample(1:n.pheno, trans.size, replace = FALSE)
        # common perception is that trans eqtls have smaller effet sizes
        # so here we are nerfing the effect size
        trans.indexes <- data.frame("geno" = rep(trans.mx[1,i], trans.size),
                                    "pheno" = trans.pheno,
                                    "effect" = sample(cis.coeff, trans.size, replace = FALSE) * trans.nerf,
                                    "label" = "trans")
        eqtl.indexes <- rbind(eqtl.indexes, trans.indexes)
    }

    for(i in 1:nrow(eqtl.indexes)){
        eqtl.effect[eqtl.indexes[i,"geno"], eqtl.indexes[i,"pheno"]] <- eqtl.indexes[i,"effect"]
    }

    pheno.list[["noisefree"]] <- input.geno %*% eqtl.effect

    # key results = phenotype.mx, eqtl.indexes, eqtl.effect

    sparse.effect <- 1*(eqtl.effect != 0)

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

    #eqtl.indexes = generative.truth

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


    simdetails <- data.frame("N_pheno" = n.pheno, "N_eqtl" = n.eqtl, "eqtl_coeff" = coeff.mean,
                     "cis_trans_ratio" = cis.trans.ratio, "trans_impact" = trans.impact,
                     "trans_nerf" = trans.nerf, "N_hf" = length(factors), "HF_frac" = hf.frac,
                     "HF_effect" = effect.size)


    dir.create(output.path, showWarnings = FALSE)

    # Create hdf5 file

    h5createFile(output_h5)


    # creating structure of hdf5 file
    level1.groups <- c("phenotypes", "genotypes", "covars")

    for(l1 in 1:length(level1.groups)){
        h5createGroup(output_h5, level1.groups[l1])
        h5createGroup(output_h5, paste(level1.groups[l1], "col_info", sep = "/"))
        h5createGroup(output_h5, paste(level1.groups[l1], "row_info", sep = "/"))
    }

    h5createGroup(output_h5, "sim_info")
    h5createGroup(output_h5, "ROC_df")
    # write data from simulation to hdf5 file

    h5createDataset(output_h5, "genotypes/matrix", c(nrow(input.geno), ncol(input.geno)), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/matrix", c(nrow(pheno.list[["hf"]]), ncol(pheno.list[["hf"]])), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/noisefree", c(nrow(pheno.list[["hf"]]), ncol(pheno.list[["hf"]])), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/noise", c(nrow(pheno.list[["hf"]]), ncol(pheno.list[["hf"]])), chunk = NULL, level = 0)

    h5createDataset(output_h5, "sim_info/sig_map", c(ncol(sparse.effect), nrow(sparse.effect)), chunk = NULL, level = 0)

    h5write(colnames(pheno.list[["hf"]]), output_h5, "phenotypes/col_info/id")
    h5write(rownames(pheno.list[["hf"]]), output_h5, "phenotypes/row_info/id")
    h5write(colnames(input.geno), output_h5, "genotypes/col_info/id")
    h5write(rownames(input.geno), output_h5, "genotypes/row_info/id")

    h5write(pheno.list[["hf"]], output_h5, "phenotypes/matrix")
    h5write(pheno.list[["noisefree"]], output_h5, "phenotypes/noisefree")
    h5write(pheno.list[["noise"]], output_h5, "phenotypes/noise")

    h5write(input.geno, output_h5, "genotypes/matrix")

    h5write(t(sparse.effect), output_h5, "sim_info/sig_map")
    h5write(ground.truth, output_h5, "sim_info/ground_truth")
    h5write(eqtl.indexes, output_h5, "sim_info/generative_truth")
    h5write(simdetails, output_h5, "sim_info/sim_details")


}

