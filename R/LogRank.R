#' binarize a vector by the mediate value
#'
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#' @return The binarized vector
#' 
binarize <- function(x=NULL, na.rm = TRUE)
{
  res = rep(0, length(x))
  m = median(x, na.rm = TRUE)
  index = which(x > m)
  res[index] = 1
  if (! na.rm) {
    index = which(is.na(x))
    res[index] = NA
  }
  return(res)
}

#' Evaluate the predicted risk scores from benchmark methods
#'
#' @import survival
#' @param data A numeric data frame or metrix indicating the predicted risk scores from benchmark methods
#' @param survival A data frame holds survival time and event status
#' @return The p-values of the benchmark methods
#' @examples
#' \dontrun{
#' res = LogRank(data=resMatrix[[i]], survival)
#' }
#' 
#' @seealso \code{\link{survdiff}}
#' @export
LogRank <- function(data, survival) {
  
  methods = toupper(colnames(data)[1:ncol(data)])
  
  samples = intersect(rownames(data), rownames(survival))
  data = data[samples, ]
  survival = survival[samples, ]
  
  cat("Calculating on",nrow(data), "samples \n")
  
  setT = survival$time
  setE = survival$event
  
  pvalues <- apply(X = as.matrix(data), MARGIN = 2, 
                  FUN = function(x, y, z) { 
                    group = binarize(x)
                    clusterNum=length(unique(group))
                    if (clusterNum>1){
                      sdf=survdiff(Surv(y, z) ~ group)
                      p_value=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
                    }
                    else
                    {
                      p_value=1
                    }
                    return(p_value)},
                  y = setT, z = setE )
  return(pvalues)
}
  
  