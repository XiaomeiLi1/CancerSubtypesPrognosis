#' Evaluate Cancer Prognosis of coding RNA expression data to compute the risk scores in 13 methods
#' 
#' Cancer Prognosis to compute the risk scores for coding RNA expression data based on existing 13 methods: 
#'            "AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx","TAMR13", "PIK3CAGS", "GENE70", "rorS", "RNAmodel","Ensemble"
#' @import genefu   
#' @import survcomp   
#' @param data data to be computed for Cancer Prognosis risk scores; either a data matrix or ExpressionSet object. If it is a data matrix, rows= probes/genes annotated with Entrez ID and columns=terms/samples .  
#' @param platform The technical platform for data, "affy" or "agilent" or "custom". Defualt "custom" represnts unknow
#' @param methods A set of methods to be performed in the 13 methods: "AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx","TAMR13", "PIK3CAGS", "GENE70", "rorS", "RNAmodel","Ensemble",
#'               The "Ensemble" method is the average of the five methods, i.e. "GENIUS", "Endopredict", "OncotypeDx", "GENE70", "rorS". 
#' @param RNASig A dataframe represnts the customized signatures provided for the "RNAmodel" method. The columns should include "Gene.symbol","EntrezGene.ID", "weight". Default NULL, the default 30 RNA signatures will be used. 
#'
#' @return
#' A dataframe object with rows for samples and columns which represent dataset used and its corresponding methods.
# \itemize{
#  \item \strong{AURKA}: Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, et al. Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes. Clinical cancer research. 2008;14(16):5158-5165.
#  \item \strong{ESR1}: Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, et al. Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes. Clinical cancer research. 2008;14(16):5158-5165.
#  \item \strong{ERBB2}: Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, et al. Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes. Clinical cancer research. 2008;14(16):5158-5165.
#  \item \strong{GGI}: Sotiriou C, Wirapati P, Loi S, Harris A, Fox S, Smeds J, et al. Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis. Journal of the National Cancer Institute. 2006;98(4):262-272.
#  \item \strong{Endopredict}: Filipits M, Rudas M, Jakesz R, Dubsky P, Fitzal F, Singer CF, et al. A new molecular predictor of distant recurrence in ER-positive, HER2-negative breast cancer adds independent information to conventional clinical risk factors. Clinical Cancer Research. 2011;17(18):6012-6020.
#  \item \strong{OncotypeDx}: Paik S, Shak S, Tang G, Kim C, Baker J, Cronin M, et al. A multigene assay to predict recurrence of tamoxifen-treated, node-negative breast cancer. New England Journal of Medicine. 2004;351(27):2817-2826.
#  \item \strong{TAMR13}: Loi S, Haibe-Kains B, Desmedt C, Wirapati P, Lallemand F, Tutt AM, et al. Predicting prognosis using molecular profiling in estrogen receptor-positive breast cancer treated with tamoxifen. BMC genomics. 2008;9(1):239.
#  \item \strong{PIK3CAGS}: Loi S, Haibe-Kains B, Majjaj S, Lallemand F, Durbecq V, Larsimont D, et al. PIK3CA mutations associated with gene signature of low mTORC1 signaling and better outcomes in estrogen receptor{positive breast cancer. Proceedings of the  National Academy of Sciences. 2010;107(22):10208-10213.
#  \item \strong{GENE70}: Van't Veer LJ, Dai H, Van De Vijver MJ, He YD, Hart AA, Mao M, et al. Gene expression profiling predicts clinical outcome of breast cancer. nature. 2002;415(6871):530.
#  \item \strong{rorS}: Parker JS, Mullins M, Cheang MC, Leung S, Voduc D, Vickery T, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. Journal of clinical oncology. 2009;27(8):1160.
#  \item \strong{RNAmodel}: Volinia S, Croce CM. Prognostic microRNA/mRNA signature from the integrated analysis of patients with invasive breast cancer. Proceedings of the National Academy of Sciences. 2013;110(18):7413-7417.
#  \item \strong{Ensemble}: Ensemble calculates the average of predicted risk scores from 5 methods, GENIUS, EndoPredict, OncotypeDx, GENE70, and rorS. We chose these five methods because of their good performance.
#  }
#'  
#' @examples 
#' \dontrun{
#' library("breastCancerMAINZ")
#' data("mainz")
#' methods <- c("AURKA", "ESR1", "ERBB2", "GGI", "GENIUS", "Endopredict", "OncotypeDx", 
#'             "TAMR13", "PIK3CAGS", "GENE70", "rorS", "RNAmodel", "Ensemble")
#' res = CancerPrognosis_RNAData(data=mainz, platform="custom", methods=methods)
#'}
#'
#' @export
#' 
#' 
CancerPrognosis_RNAData <- function(data, platform="custom", methods, RNASig=NULL) 
{
  #data(genefu::scmgene.robust)
  scmgene.robust = genefu::scmgene.robust
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
      dannot = data.frame("EntrezGene.ID"=colnames(ddata)) 
    }
    
    # do.mapping TRUE if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise. default = "FALSE"
    
    do.mapping = ifelse(platform == "affy", FALSE, TRUE)
    
    ## 1. AURKA
    if ("AURKA" %in% methods) {
      cat("Processing on AURKA.... \n")
      modt  = as.data.frame(genefu::scmgene.robust$mod$AURKA)
      ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
      #due to transcriptome name of agilent platform (NM_)
      if (platform == "agilent")
        modt[, "probe"] = "NM_003600"
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({
          sig.score(x = modt,data = ddata, annot = dannot,do.mapping = domap1)$score}, error = function(e) {NA})
        if ((!is.na(tmp1)) || t == 2)
          break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for AURKA signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "AURKA" = tmp1)
    }
    
    ## 2. ESR1
    if ("ESR1" %in% methods) {
      cat("Processing on ESR1.... \n")
      modt = as.data.frame(genefu::scmgene.robust$mod$ESR1)
      # if agilent platform consider the probe published in Desmedt et al., CCR, 2008
      if (platform == "agilent")
        modt[, "probe"] = "NM_000125"
      
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({
          sig.score( x = modt,data = ddata,annot = dannot,do.mapping = domap1)$score}, error = function(e) {NA})
        if ((!is.na(tmp1)) || t == 2)
          break
        domap1 =  !domap1
        t = t + 1
      }
      if (is.na(tmp1))
      {
        cat("Failed ... Check your data for ESR1 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "ESR1" = tmp1)
    }
    
    ## 3. ERBB2
    if ("ERBB2" %in% methods) {
      ## ERBB2
      cat("Processing on ERBB2.... \n")
      modt = as.data.frame(genefu::scmgene.robust$mod$ERBB2)
      ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
      if(platform == "agilent") modt[ , "probe"] = "NM_004448"
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ sig.score(x=modt, data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1))|| t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for ERBB2 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "ERBB2"=tmp1)
    }
    
    ## 4. GGI
    if ("GGI" %in% methods) {
      cat("Processing on GGI.... \n")
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ ggi(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for GGI signatures \n")
        tmp1 = rep(NA, ncol(data))
      } 
      rest = cbind(rest, "GGI"=tmp1)
    }
    
    ## 5. GENIUS
    if ("GENIUS" %in% methods) {
      cat("Processing on GENIUS.... \n")
      
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ genius(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1))  || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for GENIUS signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "GENIUS"=tmp1)
    }
    
    # which(colnames(ddata) == "RBBP8")
  
    ## 6. ENDOPREDICT
    if ("Endopredict" %in% methods) {
      cat("Processing on EndoPredict.... \n")
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ endoPredict(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for Endopredict signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "Endopredict"=tmp1)
    }
    
    ## 7. oncotypedx
    if ("OncotypeDx" %in% methods) {
      cat("Processing on OncotypeDx.... \n")
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ oncotypedx(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for OncotypeDx signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "OncotypeDx"=tmp1)
    }
    
    
    ## 8. tamr13
    if ("TAMR13" %in% methods) {
      cat("Processing on TAMR13.... \n")
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ tamr13(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for TAMR13 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "TAMR13"=tmp1)
    }
    
    ## 9. PIK3CAGS
    if ("PIK3CAGS" %in% methods) {
      ## PIK3CAGS
      cat("Processing on PIK3CAGS.... \n")
      
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ pik3cags(data=ddata, annot=dannot, do.mapping=domap1) }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for PIK3CAGS signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "PIK3CAGS"=tmp1)
    }
    
    ## 10. GENE70
    if ("GENE70" %in% methods) {
      cat("Processing on GENE70.... \n")
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ gene70(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for GENE70 signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "GENE70"=tmp1)
    }
    
    ## 11. RORS
    if ("rorS" %in% methods) {
      cat("Processing on rorS.... \n")
      #please note to provide the corresponding gene names in the data matrix
      domap1 = do.mapping
      t = 1
      repeat {
        tmp1 = tryCatch({ rorS(data=ddata, annot=dannot, do.mapping=domap1)$score }, error = function(e) { NA })
        if ((!is.na(tmp1)) || t==2) break
        domap1 = !domap1
        t = t + 1
      }
      if (is.na(tmp1)) {
        cat("Failed ... Check your data for rorS signatures \n")
        tmp1 = rep(NA, ncol(data))
      }
      rest = cbind(rest, "rorS"=tmp1)
    }
    
    ## 12. ENSEMBLE
    if ("Ensemble" %in% methods) {
      cat("Processing on Ensemble method.... \n")
      temp <- scale(rest)
      comb = c("GENIUS", "Endopredict", "OncotypeDx", "GENE70","rorS")
      x = match(comb, colnames(rest))
      x = x[!is.na(x)]
      temp2 <- temp[, x] # Use good methods with relatively high C-index
      ensemble <- rowMeans(temp2, na.rm = T) #average of 5 methods
      if (is.na(ensemble)) {
        cat("Failed ... Check your data for Ensemble method \n")
        ensemble = rep(NA, ncol(data))
      }
      #else colnames(ensemble) = "Ensemble"
      rest <- cbind(rest, "Ensemble"=ensemble)
    }
    
    ## 13. RNAMODEL
    if ("RNAmodel" %in% methods) {
      cat("Processing on RNAmodel.... \n")
      # if (missing(RNASig)) {
      #   cat("Missing RNA Signature list")
      # } 
      # else {
        ## Prognostic microRNA/mRNA signature from the integrated analysis of patients with invasive breast cancer
        domap = ifelse(platform == "agilent", TRUE, FALSE)
        tmp1 = tryCatch({as.matrix(RNAmodel(data=t(ddata), annot=dannot,RNASignature=RNASig))}, error = function(e) {NA})
        if (is.na(tmp1)) {
          cat("Failed ... Check your data for RNAmodel signatures \n")
          tmp1 = rep(NA, ncol(data))
        }
        else colnames(tmp1) = "RNAmodel"
        rest = cbind(rest, "RNAmodel"= tmp1)
    #   }
     }
  }
  return(res = as.data.frame(rest))
}
