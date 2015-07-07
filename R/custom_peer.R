#' Run PEER
#'
#' Custom PEER function that generates factor coefficients and residuals.
#'
#' @param expr A gene expression matrix with dimensions N x g.
#' @param covars Known fixed effects that you wish to account for.
#' @param n.factors Number of factors to be estimated, if left unspecified 25% of the sample size
#'        is going to be used. (round(N * 0.25))
#'
#' @return List with residuals, weights and factors saved.
#' @keywords keywords
#'
#' @import peer
#' @export
custom_peer <- function(expr, n.factors = NULL, covars = NULL){
  model <- peer::PEER()
  if(is.null(n.factors)){
    n.factors <- round(nrow(expr) * 0.25)
  }

  peer::PEER_setNk(model, n.factors)
  peer::PEER_setPhenoMean(model, expr)

  if(!is.null(covars)){
    peer::PEER_setCovariates(model, covars)
  }

  PEER_update(model)

  peer.results <- list()
  # get factors
  peer.results[["factors"]] <- PEER_getX(model)
  colnames(peer.results[["factors"]]) <- paste("factor",c(1:ncol(peer.results[["factors"]])), sep = "_")
  rownames(peer.results[["factors"]]) <- rownames(expr)

  # get weights
  peer.results[["weights"]] <- PEER_getW(model)
  colnames(peer.results[["weights"]]) <- colnames(peer.results[["factors"]])
  rownames(peer.results[["weights"]]) <- colnames(expr)

  # get ARD parameters
  peer.results[["ARD"]] <- PEER_getAlpha(model)

  peer.results[["residuals"]] <- PEER_getResiduals(model)


  colnames(peer.results[["residuals"]]) <- colnames(expr)
  rownames(peer.results[["residuals"]]) <- rownames(expr)

  return(peer.results)

}
