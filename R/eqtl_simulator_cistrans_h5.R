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
eqtl_simulator_cistrans_h5 <- function(input_h5 = NULL,
                                       n_pheno = 2000,
                                       n_geno = NULL,
                                       n_eqtl = 1500,
                                       cis_trans_ratio = 0.7,
                                       simulation.id = "sim1",
                                       output.path = "./",
                                       coeff_mean = 2.5, # now it is mean of normal distribution
                                       trans_nerf = 0.7,
                                       hot_spots = 8,
                                       hot_spot_impact = 0.1,
                                       hidden.factor = TRUE,
                                       factors = "sparse",
                                       factor_coeff = 2.1,
                                       hf_frac = 0.2,
                                       effect_size = 2,
                                       factor_type = "sparse"){


    output_h5 <- paste(simulation.id,".h5",sep = "")
    output_h5 <- paste(output.path, output_h5, sep = "")

    pheno.list <- list()


    geno_mx <- h5read(input_h5, "genotypes/matrix")
    sample_ids <- h5read(input_h5, "genotypes/row_info/id")
    n_sample <- nrow(geno_mx)

    if(is.null(n_geno)){
        n_geno <- ncol(geno_mx)
    } else if (!is.null(n_geno)){
        geno_subset <- sort(sample(1:ncol(geno_mx), n_geno,replace = FALSE), decreasing = FALSE)
        geno_mx <- geno_mx[,geno_subset]
    }
    geno_info <- as.data.frame(h5read(input_h5, "genotypes/col_info"), stringsAsFactors = FALSE)

    pheno_info <- as.data.frame(h5read(input_h5, "phenotypes/col_info"), stringsAsFactors = FALSE)

    subset_geno_info <- geno_info[geno_subset,]
    subset_pheno_info <- pheno_info[sort(sample(1:nrow(pheno_info), n_pheno,replace = FALSE), decreasing = FALSE),]
    geno_ids <- subset_geno_info$id
    # generating a matrix of cis and trans positions based on the geno and pheno info
    subset_cis_trans = create_cis_trans_mx(geno_info = subset_geno_info, pheno_info = subset_pheno_info, cis_threshold = 100000)

    # sample from phenotypes a subset
#    subset_pheno <- sort(sample(1:ncol(cis_trans_mx), n_pheno, replace = FALSE), decreasing = FALSE)

#    subset_pheno_info <- pheno_info[subset_pheno,]
    # subset the cis trans matrix
#    subset_cis_trans <- cis_trans_mx[,subset_pheno]

    geno_hot_spots <- sample(1:nrow(subset_cis_trans),hot_spots)
    dir.create(output.path, showWarnings = FALSE)

    ################################################################################
    #################              Simulate Phenotype              #################
    ################################################################################
    #pheno.list <- list()

    cis_count <- round(n_eqtl * cis_trans_ratio)

    trans_count <- n_eqtl - cis_count

    cis_indexes <- which(subset_cis_trans, arr.ind = TRUE)
    trans_indexes <- which(!subset_cis_trans, arr.ind = TRUE)

    sample_cis <- sort(sample(1:nrow(cis_indexes), cis_count, replace = FALSE), decreasing = FALSE)
    sample_trans <- sort(sample(1:nrow(trans_indexes), trans_count, replace = FALSE), decreasing = FALSE)

    sampled_cis_idx <- cis_indexes[sample_cis,]
    sampled_trans_idx <- trans_indexes[sample_trans,]
    # creat genotype phenotype pairs for eqtl

    # sample effect of cis eQTL
    cis_coeff <- sample(c(-1,1),1) * rnorm(cis_count, mean = coeff_mean, sd = 1)
    trans_coeff <- sample(c(-1,1),1) * rnorm(trans_count, mean = coeff_mean * trans_nerf, sd = 1)

    cis_eqtl_info <- data.frame("geno" = sampled_cis_idx[,1], "pheno" = sampled_cis_idx[,2],
                                "effect" = cis_coeff, "label" = rep("cis", nrow(sampled_cis_idx)))

    trans_eqtl_info <- data.frame("geno" = sampled_trans_idx[,1], "pheno" = sampled_trans_idx[,2],
                                  "effect" = trans_coeff, "label" = rep("trans", nrow(sampled_trans_idx)))

    # create ground truth dataframe which has all the genotype / phenotype
    # relationships and effectsize saved

    eqtl_indexes <- rbind(cis_eqtl_info, trans_eqtl_info)

    hot_spot_impact_max <- round(n_geno * hot_spot_impact)

    hot_spot_geno <- sort(sample(1:n_geno, hot_spots), decreasing = FALSE)
    for(h in hot_spot_geno){
      hot_subset <- trans_indexes[which(trans_indexes[,1] == h),]
      impact_pheno <- min(hot_spot_impact_max, nrow(hot_subset))
      hot_index <- hot_subset[sample(1:nrow(hot_subset), impact_pheno),]
      hot_coeff <- sample(c(-1,1),1) * rnorm(nrow(hot_index), mean = coeff_mean * trans_nerf, sd = 1)
      temp_df <- data.frame("geno" = hot_index[,1], "pheno" = hot_index[,2],
                            "effect" = hot_coeff, "label" = rep("trans_hot", nrow(hot_index)))
      eqtl_indexes <- rbind(eqtl_indexes, temp_df)
    }

    # create eQTL effect matrix
    eqtl_effect <- matrix(0, nrow = n_geno, ncol = n_pheno)

    # loop over trans.mx matrix to create trans relationships
    # and add it to the eqlt.indexes data frame

    for(i in 1:nrow(eqtl_indexes)){
        eqtl_effect[eqtl_indexes[i,"geno"], eqtl_indexes[i,"pheno"]] <- eqtl_indexes[i,"effect"]
    }

    pheno_list <- list()
    pheno_list[["noisefree"]] <- geno_mx %*% eqtl_effect


    colnames(pheno_list[["noisefree"]]) <- as.character(subset_pheno_info$id)
    rownames(pheno_list[["noisefree"]]) <- sample_ids

    ################################################################################
    #################              Add Gaussian Noise              #################
    ################################################################################

    pheno_list[["noise"]] <- pheno_list[["noisefree"]] +
                                        rnorm(length(pheno_list[["noisefree"]]),
                                              mean = 0,
                                              sd = 0.5)


    ################################################################################
    #################              Add Hidden factors              #################
    ################################################################################

    factor_details <- matrix(cbind(factors,factor_coeff),
                                   ncol = 2,
                                   byrow = FALSE)

    hf_list <- apply(factor_details, 1, function(x) hf_sim(n.genes = n_pheno,
                                                         n.samples = n_sample,
                                                         hf.type = x[1],
                                                         coeff.dist = x[2],
                                                         fraction.affected = hf_frac,
                                                         factor.effect.size = effect_size))

    # add all hidden factor effects
    hf_effect <- Reduce(`+`, lapply(hf_list, function(x) x$effect))

    pheno_list[["hf"]] <- pheno_list[["noise"]] + t(hf_effect)


    simdetails <- data.frame("N_pheno" = n_pheno, "N_eqtl" = n_eqtl, "eqtl_coeff" = coeff_mean,
                     "cis_trans_ratio" = cis_trans_ratio,
                     "trans_nerf" = trans_nerf, "N_hf" = length(factors), "HF_frac" = hf_frac,
                     "HF_effect" = effect_size, "HF_type" = factor_type)


    dir.create(output.path, showWarnings = FALSE)

    # Create hdf5 file
	  create_eqtl_input_h5(output_h5)
#    h5createFile(output_h5)


    # creating structure of hdf5 file
#    level1.groups <- c("phenotypes", "genotypes", "covars")

#    for(l1 in 1:length(level1.groups)){
#        h5createGroup(output_h5, level1.groups[l1])
#        h5createGroup(output_h5, paste(level1.groups[l1], "col_info", sep = "/"))
#        h5createGroup(output_h5, paste(level1.groups[l1], "row_info", sep = "/"))
#    }

    h5createGroup(output_h5, "sim_info")
    h5createGroup(output_h5, "ROC_df")
    # write data from simulation to hdf5 file

    h5createDataset(output_h5, "genotypes/matrix", c(n_sample, n_geno), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/matrix", c(n_sample, n_pheno), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/noisefree", c(n_sample, n_pheno), chunk = NULL, level = 0)
    h5createDataset(output_h5, "phenotypes/noise", c(n_sample, n_pheno), chunk = NULL, level = 0)

#    h5createDataset(output_h5, "sim_info/sig_map", c(ncol(sparse.effect), nrow(sparse.effect)), chunk = NULL, level = 0)

    h5write(colnames(pheno_list[["hf"]]), output_h5, "phenotypes/col_info/id")
    h5write(rownames(pheno_list[["hf"]]), output_h5, "phenotypes/row_info/id")
    h5write(as.character(subset_pheno_info$pheno_chr), output_h5, "phenotypes/col_info/pheno_chr")
    h5write(as.double(subset_pheno_info$pheno_start), output_h5, "phenotypes/col_info/pheno_start")
    h5write(as.double(subset_pheno_info$pheno_end), output_h5, "phenotypes/col_info/pheno_end")
    h5write(as.character(subset_geno_info$geno_chr), output_h5, "genotypes/col_info/geno_chr")
    h5write(as.double(subset_geno_info$geno_pos), output_h5, "genotypes/col_info/geno_pos")

    h5write(geno_ids, output_h5, "genotypes/col_info/id")
    h5write(sample_ids, output_h5, "genotypes/row_info/id")

    h5write(pheno_list[["hf"]], output_h5, "phenotypes/matrix")
    h5write(pheno_list[["noisefree"]], output_h5, "phenotypes/noisefree")
    h5write(pheno_list[["noise"]], output_h5, "phenotypes/noise")

    h5write(geno_mx, output_h5, "genotypes/matrix")

    #h5write(t(sparse.effect), output_h5, "sim_info/sig_map")
    h5createGroup(output_h5, "sim_info/ground_truth")
    h5write(as.numeric(eqtl_indexes$geno), output_h5, "sim_info/ground_truth/geno")
    h5write(as.numeric(eqtl_indexes$pheno), output_h5, "sim_info/ground_truth/pheno")
    h5write(as.numeric(eqtl_indexes$effect), output_h5, "sim_info/ground_truth/effect")
    h5write(as.character(eqtl_indexes$label), output_h5, "sim_info/ground_truth/cis_trans")
    #h5write(eqtl.indexes, output_h5, "sim_info/generative_truth")
    h5write(simdetails, output_h5, "sim_info/sim_details")

}



#' @export

create_cis_trans_mx <- function(geno_info = NULL, pheno_info = NULL, cis_threshold = 1000000){
  cis_threshold <- 100000
  cis_trans_mx <- matrix(NA, nrow = nrow(geno_info), ncol = nrow(pheno_info))

  for(i in 1:nrow(geno_info)){

    same_chr <- as.character(geno_info[i,"geno_chr"]) == as.character(pheno_info[,"pheno_chr"])

    distance <- cbind(abs(pheno_info[,"pheno_start"] - geno_info[i,"geno_pos"]), abs(pheno_info[,"pheno_end"] - geno_info[i,"geno_pos"]))

    min_distance <- apply(distance, 1, min)

    cis_distance <- min_distance <= cis_threshold

    cis_trans <- same_chr & cis_distance

    cis_trans_mx[i,] <- cis_trans

  }

  return(cis_trans_mx)
}
