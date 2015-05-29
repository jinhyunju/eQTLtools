#' @export
cis_labeling <- function(hit.info, cis_distance){
  pheno_start <- as.numeric(hit.info["pheno_start"])
  pheno_end <- as.numeric(hit.info["pheno_end"])
  geno_pos <- as.numeric(hit.info["geno_pos"])

  if( hit.info["geno_chr"] == hit.info["pheno_chr"]){
    distance <- min(abs(pheno_start - geno_pos), abs(pheno_end - geno_pos))
    if(distance <= cis_distance){
      label <- "cis"
    } else {
      label <- "trans"
    }
  } else {
    label <- "trans"
  }

  return(label)

}
