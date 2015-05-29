#' @import gtools
cleaning_results <- function(input.df, snp.info, gene.info){

    #input.df$geno_id <- snp.info[input.df$genotype.idx, "geno_id"]
    input.df <- merge(input.df, snp.info, by.x = "genotype",by.y = "geno_id")
    input.df <- merge(input.df, gene.info, by = "phenotype")
    input.df$hit_id <- paste(input.df$genotype, input.df$phenotype, sep = "_")
    result.df <- subset(input.df, select= c(hit_id, pheno_symbol,geno_chr, geno_pos, pheno_chr, pheno_start,pheno_end, pval))
    result.df <- result.df[order(result.df$pval),]
    result.df$cis_trans <- apply(result.df, 1, cis_labeling)
    result.df$geno_chr <- factor(result.df$geno_chr, levels = mixedsort(unique(result.df$geno_chr)))
    result.df$pheno_chr <- factor(result.df$pheno_chr, levels = mixedsort(unique(result.df$pheno_chr)))

    return(result.df)
}
