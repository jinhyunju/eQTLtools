# Function that calculates the area under a ROC or PR curve
# trapezoid area calculations
# assumes that curves start at 0 and end at 1
#' @export
#'

auc_calc <- function(data.pr = NULL, type = c("ROC","PR")) {

  if ( is.null(data.pr) ) {
    return ( NA )
  }

  if(type == "ROC"){
    tpr <- data.pr$TPR
    fpr <- data.pr$FPR

    # Measure the area under the curve for the remaining points along the PR curve
    # trapezoid
    dr <- fpr[-1] - fpr[-length(fpr)]
    auc <- sum(dr * (tpr[-1]+tpr[-length(tpr)])/2)
  } else {
    rec <- data.pr$TPR
    pre <- data.pr$precision

    # Measure the area under the curve for the remaining points along the PR curve
    # trapezoid
    dr <- rec[-1] - rec[-length(rec)]
    auc <- sum(dr * (pre[-1]+pre[-length(pre)])/2)


  }
  return ( auc )
}
