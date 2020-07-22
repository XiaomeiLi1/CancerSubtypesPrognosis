#' Evaluate Cancer Prognosis of miRNA expression data to compute the risk scores in 4 methods
#' 
#' Cancer Prognosis to compute the risk scores for miRNA RNA expression data based on existing 4 methods: 
#'            "hsa-miR-210", "hsa-miR-155", "hsa-miR-335", "miRNA10"
#' @import Biobase
#' @import stringr
#' @param data data to be computed for Cancer Prognosis risk scores; either a data matrix or ExpressionSet object. If it is a data matrix, rows= miRNAs and columns=terms/samples .  
#' @param methods A set of methods to be performed in the 4 methods: "hsa-miR-210", "hsa-miR-155", "hsa-miR-335", "miRNA10"
#' @examples 
#' data(TCGA500)
#' methods <- c("hsa-miR-210", "hsa-miR-155", "hsa-miR-335", "miRNA10")
#' res = CancerPrognosis_miRNAData(data=TCGA500, methods=methods)
#' @return 
#' A dataframe object with rows for samples and columns which represent dataset used and its corresponding methods
#' 
#' @references 
#' \itemize{
#'  \item \strong{hsa-miR-210}: Lee JA, Lee HY, Lee ES, Kim I, Bae JW. Prognostic implications of microRNA-21 overexpression in invasive ductal carcinomas of the breast. Journal of breast cancer. 2011;14(4):269-275. Yan LX, Huang XF, Shao Q, Huang MY, Deng L, Wu QL, et al. MicroRNA miR-21 overexpression in human breast cancer is associated with advanced clinical stage, lymph node metastasis and patient poor prognosis. Rna. 2008;. Markou A, Yousef GM, Stathopoulos E, Georgoulias V, Lianidou E. Prognostic significance of metastasis-related microRNAs in early breast cancer patients with a long follow-up. Clinical chemistry. 2013; p. clinchem-2013.
#'  \item \strong{hsa-miR-155}: Gasparini P, Cascione L, Fassan M, Lovat F, Guler G, Balci S, et al. microRNA expression profiling identifies a four microRNA signature as a novel diagnostic and prognostic biomarker in triple negative breast cancers. Oncotarget. 2014;5(5):1174.
#'  \item \strong{hsa-miR-335}: Andorfer CA, Necela BM, Thompson EA, Perez EA. MicroRNA signatures: clinical biomarkers for the diagnosis and treatment of breast cancer. Trends in molecular medicine. 2011;17(6):313-319.
#'  \item \strong{miRNA10}: Buffa FM, Camps C, Winchester L, Snell CE, Gee HE, Sheldon H, et al. microRNA associated progression pathways and potential therapeutic targets identified by integrated mRNA and microRNA expression profiling in breast cancer. Cancer research. 2011; p. canres-0489.
#' }
#' 
#' @export
#' 
#' 
CancerPrognosis_miRNAData <- function(data, methods) 
{
  if (!class(data) %in% c("matrix", "ExpressionSet")) {
    stop("data must be a matrix object or ExpressionSet (eset object)")
  }
  else
  {
    cat("***** Computing on risk scores ***** \n")
    
    if(class(data)=="ExpressionSet")
    {
      ddata = t(exprs(data))
      dannot = fData(data)
    }
    else
    {
      ddata = t(data)
      dannot = data.frame("Gene.symbol"=colnames(ddata)) 
    }
    
    rest <- NULL
    
    ##1. hsa-miR-210
    if ("hsa-miR-210" %in% methods) {
      cat("Processing on hsa-miR-210.... \n")
      single.sig <- "hsa-miR-210"
      #Coefficient from univariate cox regression
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) {
        single.sig <- c("hsa-miR-210-3p","hsa-miR-210-5p")
        index = which(toupper(dannot$Gene.symbol) %in% toupper(single.sig))
      }
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1){
        tdata = apply(ddata[,index],1,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for hsa-miR-210 signatures \n")
      }
      rest = cbind(rest, "miR-210" = tmp1)
    }
    
    ## 2. hsa-miR-155
    if ("hsa-miR-155" %in% methods) {
      cat("Processing on hsa-miR-155.... \n")
      single.sig <- "hsa-miR-155"
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) {
        single.sig <- c("hsa-miR-155-3p","hsa-miR-155-5p")
        index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      }
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1){
        tdata = apply(ddata[,index],1,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1))
      {
        cat("Failed ... Check your data for hsa-miR-155 signatures \n")
      }
      rest = cbind(rest, "miR-155" = tmp1)
    }
    
    ## 3. hsa-miR-335
    if ("hsa-miR-335" %in% methods) {
      cat("Processing on hsa-miR-335.... \n")
      single.sig <- "hsa-miR-335"
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) {
        single.sig <- c("hsa-miR-335-3p","hsa-miR-335-5p")
        index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      }
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1){
        tdata = apply(ddata[,index],1,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for hsa-miR-335 signatures \n")
      }
      rest = cbind(rest, "miR-335"=tmp1)
    }
    
    ## 4. miRNA10
    if ("miRNA10" %in% methods) {
      cat("Processing on miRNA10model .... \n")
      tmp1 = tryCatch({as.matrix(miRNA10model(data=t(ddata), annot=dannot))}, error = function(e) {NA})
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for miRNA10 signatures \n")
        tmp1 = rep(NA, nrow(dannot))
      } 
      else colnames(tmp1) = "miRNA10"
      rest = cbind(rest, "miRNA10"=tmp1)
    }
  }
  return(res = as.data.frame(rest))
}
