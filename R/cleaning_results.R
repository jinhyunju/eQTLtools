#' @import gtools
#' @export
cleaning_results <- function(input.df, snp.info, gene.info, cis_distance){

    #input.df$geno_id <- snp.info[input.df$genotype.idx, "geno_id"]
    input.df <- merge(input.df, snp.info,
                      by.x = "genotype",by.y = "geno_id")
    input.df <- merge(input.df, gene.info,
                      by.x = "phenotype", by.y = "pheno_id")
    input.df$hit_id <- paste(input.df$genotype,
                             input.df$phenotype, sep = "_")
    if( "p.bh" %in% colnames(input.df)){
      result.df <- subset(input.df,
                          select= c(hit_id, pheno_symbol,
                                    geno_chr, geno_pos, pheno_chr,
                                    pheno_start,pheno_end, pval, p.bh))

    } else {
      result.df <- subset(input.df,
                          select= c(hit_id, pheno_symbol,
                                    geno_chr, geno_pos, pheno_chr,
                                    pheno_start,pheno_end, pval))

    }
    result.df$pheno_symbol <- as.character(result.df$pheno_symbol)
    result.df$pheno_chr <- as.character(result.df$pheno_chr)
    result.df$geno_chr <- as.character(result.df$geno_chr)

    result.df <- result.df[order(result.df$pval),]

    cis_trans_result <- apply(result.df, 1,
                               function(x) cis_labeling(x, cis_distance))



#    cis_trans_result2 <- sapply(1:nrow(result.df),
#                                 function(x) cis_labeling(result.df[x,], cis_distance))

    result.df$cis_trans <- cis_trans_result
    result.df$geno_chr <- factor(result.df$geno_chr,
                                 levels = mixedsort(unique(result.df$geno_chr)))
    result.df$pheno_chr <- factor(result.df$pheno_chr,
                                  levels = mixedsort(unique(result.df$pheno_chr)))

    return(result.df)
}
