#' Evaluate Cancer Prognosis of Long non-coding RNA expression data to compute the risk scores 
#' 
#' Cancer Prognosis to compute the risk scores for long non-coding RNA data based on the 6 methods: "HOTAIR", "MALAT1", "DSCAM-AS1",
#'            "lncRNA5","lncRNA6","lncRNA12".
#' @param data data to be computed for cancer Prognosis risk scores; either a data matrix or ExpressionSet object. 
#' If it is a data matrix, rows= lncRNAs annotated with lncRNA symbols and columns=tems/samples .  
#' @param methods A set of methods to be performed in the 6 methods:"HOTAIR", "MALAT1", "DSCAM-AS1", "lncRNA5","lncRNA6","lncRNA12".
#' @examples 
#' data(TCGA500)
#' methods <- c("HOTAIR", "MALAT1", "DSCAM-AS1", "lncRNA12","lncRNA6","lncRNA5")
#' res = CancerPrognosis_LncRNAData(data=TCGA500, platform="custom", methods=methods)
#' @return 
#' A dataframe object with rows for samples and columns which represent dataset used and its corresponding methods
#' @references 
#' \itemize{
#'  \item \strong{HOTAIR}: Paw lowska E, Szczepanska J, Blasiak J. The Long Noncoding RNA HOTAIR in Breast Cancer: Does Autophagy Play a Role? International journal of molecular sciences. 2017;18(11):2317.
#'  \item \strong{MALAT1}: Wang Z, Katsaros D, Biglia N, Shen Y, Fu Y, Loo LW, et al. High expression of long non-coding RNA MALAT1 in breast cancer is associated with poor relapse-free survival. Breast cancer research and treatment. 2018; p. 1-11.
#'  \item \strong{DSCAM-AS1}: Niknafs YS, Han S, Ma T, Speers C, Zhang C, Wilder-Romans K, et al. The lncRNA landscape of breast cancer reveals a role for DSCAM-AS1 in breast cancer progression. Nature communications. 2016;7:12791.
#'  \item \strong{lncRNA12}: Zhou M, Zhong L, Xu W, Sun Y, Zhang Z, Zhao H, et al. Discovery of potential prognostic long non-coding RNA biomarkers for predicting the risk of tumor recurrence of breast cancer patients. Scientific reports. 2016;6:31038.
#'  \item \strong{lncRNA6}: Zhong L, Lou G, Zhou X, Qin Y, Liu L, Jiang W. A six-long non-coding RNAs signature as a potential prognostic marker for survival prediction of ER-positive breast cancer patients. Oncotarget. 2017;8(40):67861.
#'  \item \strong{lncRNA5}: Li J, Wang W, Xia P, Wan L, Zhang L, Yu L, et al. Identification of a five-lncRNA signature for predicting the risk of tumor recurrence in breast cancer patients. International journal of cancer. 2018.
#'  }
#'  
#' @export

# @param RNASig A dataframe represnts the customized Signatures provided for the LncRNAmodel methods. The columns should include "Ensembl.ID",  "Gene.symbol", "weight". Default NULL, the "lncRNA5model","lncRNA6model","lncRNA12model" methods will use the default Signature sets. 

CancerPrognosis_LncRNAData <- function(data,platform="custom",methods) 
{
  ## data loading
  res = NULL
  rest = NULL
  
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
    
    do.mapping = ifelse(platform == "affy", TRUE, FALSE)
    
    ## 1. HOTAIR
    if ("HOTAIR" %in% methods) {
      cat("Processing on HOTAIR.... \n")
      single.sig <- "HOTAIR"
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1) {
        tdata = apply(data[index,],2,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for HOTAIR signatures \n")
      }
      rest = cbind(rest, "HOTAIR"=tmp1)
    }

    ##2. MALAT1
    if ("MALAT1" %in% methods) {
      cat("Processing on MALAT1.... \n")
      single.sig <- "MALAT1"
      #Coefficient from univariate cox regression
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1) {
        tdata = apply(data[index,],2,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for MALAT1 signatures \n")
      }
      rest = cbind(rest, "MALAT1" = tmp1)
    }
    
    ## 3. DSCAM-AS1
    if ("DSCAM-AS1" %in% methods) {
      cat("Processing on DSCAM-AS1.... \n")
      single.sig <- "DSCAM-AS1"
      wight <- 1
      index = which(toupper(dannot$Gene.symbol) == toupper(single.sig))
      if(length(index) == 0) tmp1 = rep(NA, ncol(data))
      else if (length(index)>1) {
        tdata = apply(data[index,],2,mean)
        tmp1 <- wight*tdata
      }
      else tmp1 <- wight*ddata[,index]
      if (is.na(tmp1))
      {
        cat("Failed ... Check your data for DSCAM-AS1 signatures \n")
      }
      rest = cbind(rest, "DSCAM-AS1" = tmp1)
    }
    
    ## 4. lncRNA12model
    if ("lncRNA12" %in% methods) {
      cat("Processing on lncRNA12 model.... \n")
      tmp1 = tryCatch({as.matrix(lncRNA12model(data=t(ddata), annot=dannot, do.mapping = do.mapping))}, error = function(e) {NA})
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for lncRNA12 signatures \n")
        tmp1 = rep(NA, ncol(data))
      } 
      else colnames(tmp1) = "lncRNA12"
      rest = cbind(rest, "lncRNA12"=tmp1)
    }
    
    ## 5. lncRNA6model
    if ("lncRNA6" %in% methods) {
      cat("Processing on lncRNA6 model.... \n")
      tmp1 = tryCatch({as.matrix(lncRNA6model(data=t(ddata), annot=dannot, do.mapping = do.mapping))}, error = function(e) {NA})
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for lncRNA6 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      else colnames(tmp1) = "lncRNA6"
      rest = cbind(rest, "lncRNA6"=tmp1)
    }
    
    ## 6. lncRNA5model
    if ("lncRNA5" %in% methods) {
      cat("Processing on lncRNA5 model.... \n")
      tmp1 = tryCatch({as.matrix(lncRNA5model(data=t(ddata), annot=dannot, do.mapping = do.mapping))}, error = function(e) {NA})
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for lncRNA5 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      else colnames(tmp1) = "lncRNA5"
      rest = cbind(rest, "lncRNA5"=tmp1)
    }
    
  }
  return(res = as.data.frame(rest))
}
