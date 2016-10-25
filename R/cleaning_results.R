#' @import gtools
#' @export
cleaning_results <- function(input.df, snp.info, gene.info, cis_distance){

    #input.df$geno_id <- snp.info[input.df$genotype.idx, "geno_id"]
    input.df <- merge(input.df, snp.info,
                      by.x = "genotype",by.y = "geno_id")
    result.df <- merge(input.df, gene.info,
                      by.x = "phenotype", by.y = "pheno_id")
    result.df$hit_id <- paste(result.df$genotype,
                             result.df$phenotype, sep = "_")

    result.df$pheno_symbol <- as.character(result.df$pheno_symbol)
    result.df$pheno_chr <- as.character(result.df$pheno_chr)
    result.df$geno_chr <- as.character(result.df$geno_chr)

    result.df <- result.df[order(result.df$pval),]

    same_chr <- result.df$geno_chr == result.df$pheno_chr
    distance <- with(result.df, abs(pheno_start - geno_pos))
    cis_distance <- distance <= cis_distance
    cis_hit <- same_chr & cis_distance

    result.df$cis_trans <- ifelse(cis_hit, "cis","trans")

    return(result.df)
}
