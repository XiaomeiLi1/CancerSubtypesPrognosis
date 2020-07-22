#' C-index calculation
#' 
#' Calculating Concordance Indices for the evaluation results of Cancer Prognosis methods
#' 
#' @import genefu
#' @import survcomp
#' @param data A dataframe object with rows for samples and columns represent corresponding methods. A return value from CancerPrognosis_xxx() function
#' @param survival A dataframe object which contains variables (columns) representing for survival time and event. The rows are the samples.
#' @param outputFolder (Optional) A desired folder to put the results or "output" by default
#' @param PValue if ture, return the object of function concordance.index().
#' 
#' @return
#' This function is used for its side-effect.
#' A plot for overall Concordance Index. For C-Index for each method, please refer in the outputFolder
#' 
#' 
#' @examples 
#' library("breastCancerMAINZ")
#' data("mainz")
#' methods <- c("AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx", 
#'             "TAMR13", "PIK3CAGS", "GENE70", "rorS", "RNAmodel", "Ensemble")
#' sampleInfo= pData(mainz)           
#' survival=data.frame(time=sampleInfo$t.dmfs,event=sampleInfo$e.dmfs, row.names=sampleInfo$samplename)
#' \dontrun{
#' res = CancerPrognosis_RNAData(data=mainz, platform="custom", methods=methods)
#' CIs = Cindex(data=res, survival,outputFolder="./mainz")
#' }
#' @export
#' 
Cindex <- function(data, survival, PValue = FALSE, outputFolder=NULL) {
  
  methods = toupper(colnames(data)[1:ncol(data)])
  
  samples = intersect(rownames(data), rownames(survival))
  data = data[samples, ]
  survival = survival[samples, ]
  
  cat("Calculating on",nrow(data), "samples \n")
  
  setT = survival$time
  setE = survival$event
    
  cindex <- apply(X = as.matrix(data), MARGIN = 2, 
                  FUN = function(x, y, z) { 
                    tt = survcomp::concordance.index(x = x, surv.time = y, surv.event = z, method = "noether", na.rm = TRUE)
                    if (PValue) return (tt)
                    else return(c("cindex" = tt$c.index, "cindex.se" = tt$se, "lower" = tt$lower, "upper" = tt$upper))},
                    y = setT, z = setE )
  if(! is.null(outputFolder)) {
    dir.create(file.path(outputFolder), showWarnings = FALSE)
    write.csv(cindex,file = paste(outputFolder, "Cindex.csv", sep = "/"),row.names = TRUE)
  }
  else {
    outputFolder = "./output"
    dir.create(file.path(outputFolder), showWarnings = FALSE)
    write.csv(cindex,file = paste(outputFolder, "Cindex.csv", sep = "/"),row.names = TRUE)
  }
  return(cindex)
}



