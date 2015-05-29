#' @export
sortBM_entrez <- function(ensg.input){
  #output <- getBM(attributes=c('entrezgene','ensembl_gene_id','hgnc_symbol',
  #                             'external_gene_id','gene_biotype','status',
  #                             'chromosome_name','start_position','end_position','description'),
  #                filters = 'entrezgene',values = ensg.input, mart = mart.ensembl)


  output <- getBM(attributes=c('entrezgene','hgnc_symbol',
                               'gene_biotype',#'status',
                               'chromosome_name','start_position','end_position'),
                  filters = 'entrezgene',values = ensg.input, mart = mart.ensembl)
  sorted.out <- output[order(match(output$entrezgene, ensg.input)),]
  return(sorted.out)

}
#' @export
sortBM <- function(ensg.input){
  # correcting for ENSGR gene ids
  ensg.input <- sub('R',0,ensg.input)
  #  test <- unlist(sapply(1:length(ensg.input),function (a) getBM(attributes=c('ensembl_gene_id'), filters = 'ensembl_gene_id',values = ensg.input[a], mart = ensembl)))
  output <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',
                               'gene_biotype',#'status',
                               'chromosome_name','start_position','end_position'),
                  filters = 'ensembl_gene_id',values = ensg.input, mart = mart.ensembl)
  sorted.out <- output[order(match(output$ensembl_gene_id, ensg.input)),]
  return(sorted.out)

}
