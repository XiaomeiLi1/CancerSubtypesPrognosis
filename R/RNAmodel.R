#' Evaluate Cancer Prognosis based on 37 miRNA/mRNA signatures
#' Cancer Prognosis to compute the risk scores for miRNA/mRNA data based on 37 microRNA/mRNA signatures
#' @import stringr
#' @param data data to be computed for cancer Prognosis risk scores; a data matrix, rows= miRNAs annotated with miRNA names and columns=terms/samples.
#' @param annot annot holds miRNA/mRNA names in the data.
#' @param RNASignature RNA signatures, if none, will import the 37 microRNA/mRNA signatures.
#' @examples 
#' data(TCGA500)
#' data = exprs(TCGA500)
#' annot = fData(TCGA500)
#' res =  RNAmodel(data,annot)
#' @details
#' The 7 miRNAs signatures are  "hsa-miR-103", "hsa-miR-1307", "hsa-miR-148b", "hsa-miR-328", "hsa-miR-484", "hsa-miR-874","hsa-miR-93".
#' The 30 mRNA signatures are "ACSL1" ,"ADAT1", "ANKRD52","BIRC6", "CPT1A", "CXCR7","DAAM1", "DIP2B", "FAM199X", "FAM91A1", "FRZB",  "GLA",  "GMCL1", "HRASLS",  "HSP90AA1", "MCM10", "ME1",  "NDRG1", "NOTCH2NL","OTUD6B","PDSS2", "PGK1",  "PIK3CA", "PTAR1", "SMG1",  "TRIM23", "TTC3",  "UBR5",  "UBXN7", "ZFC3H1".
#' Note that the miRbase version of miRNAs signatures is v16, please check the miRbase version of your own dataset.
#' @return 
#' A Numeric Vector
#' @references 
#' \itemize{
#'  \item \strong{RNAmodel}: Volinia S, Croce CM. Prognostic microRNA/mRNA signature from the integrated analysis of patients with invasive breast cancer. Proc Natl Acad Sci USA. 2013;110:7413‐7417.
#' }
#' @export
#' 
RNAmodel <- function(data,annot,RNASignature=NULL)
{
  if(is.null(RNASignature))
  {
    #data("sig.RNA37")
    RNASig = CancerSubtypesPrognosis::RNASig
  }
  else
    RNASig=RNASignature
  tdata = NULL
  SigW = NULL
  for(j in 1:nrow(RNASig)) {
    #if (toupper(do.dataset) %in% c("GSE6532", "GEO", "NKI")) 
    #  index = which(annot$EntrezGene.ID %in% RNASig$EntrezGene.ID[j])
    #else index = which(toupper(annot$Gene.symbol) %in% toupper(RNASig$Gene.symbol[j]))
    
    index = which(annot$EntrezGene.ID %in% RNASig$EntrezGene.ID[j])
    if(length(index!=0)) {
      ## Merge the same probles of RNAs with average value
      if (length(index)>1) tdata = cbind(tdata, apply(data[index,],2,mean))
      else
        tdata = cbind(tdata, data[index,])
      SigW = rbind(SigW, RNASig[j,])
    }
  }
  if(is.null(tdata) | is.null(SigW))
  {
    rs=NULL
  }
  else
  {
    colnames(tdata)=rownames(SigW)
    rs=tdata %*% SigW$weight
    rs <- rs - 8.877
  }
  return(rs)
}

#' Evaluate Cancer Prognosis based on 10 miRNA signatures
#' Cancer Prognosis to compute the risk scores for miRNA data based on 10 miRNA signatures
#' @import stringr
#' @param data data to be computed for cancer Prognosis risk scores; a data matrix, rows= miRNAs annotated with miRNA names and columns=terms/samples.
#' @param annot annot holds RNA names in the data.
#' @param RNASignature RNA signatures, if none, will import the 10 microRNA signatures.
#' @examples 
#' data(TCGA500)
#' data = exprs(TCGA500)
#' annot = fData(TCGA500)
#' res =  miRNA10model(data,annot)
#' @details
#' The 10 miRNAs signatures are  "hsa-miR-144", "hsa-miR-150", "hsa-miR-210", "hsa-miR-27b", "hsa-miR-30c", "hsa-miR-342"    "hsa-miR-128a", "hsa-miR-135a", "hsa-miR-767-3p", "hsa-miR-769-3p".
#' Note that the miRbase version of miRNAs signatures is v9_2, please check the miRbase version of your own dataset.
#' @return 
#' A Numeric Vector
#' @references 
#' \itemize{
#'  \item \strong{miRNA10model}: Buffa, F. M. et al. microRNA associated progression pathways and potential therapeutic targets identified by integrated mRNA and microRNA expression profiling in breast cancer. Cancer Res. 71, 5635 (2011).
#' }
#' @export
#' 
miRNA10model <- function(data,annot,RNASignature=NULL)
{
  if(is.null(RNASignature))
  {
    #data("sig.miRNA10")
    miRNA10 = CancerSubtypesPrognosis::miRNA10
  }
  else
    miRNA10=RNASignature
  tdata = NULL
  SigW = NULL
  for(j in 1:nrow(miRNA10)) {
    index = which(toupper(annot[,1]) %in% stringr::str_subset(toupper(annot[,1]),toupper(rownames(miRNA10)[j])))
    if(length(index)!=0) {
      ## Merge the same probles of RNAs with average value
      if (length(index)>1) 
        tdata = cbind(tdata, apply(data[index,],2,mean))
      else
        tdata = cbind(tdata, data[index,])
      SigW = rbind(SigW, miRNA10[j,])
    }
  }
  if(is.null(tdata) | is.null(SigW))
  {
    rs=NULL
  }
  else
  {
    colnames(tdata)=rownames(SigW)
    rs=tdata %*% SigW$weight
  }
  return(rs)
}

#' Evaluate Cancer Prognosis based on 12 lncRNA signatures
#' Cancer Prognosis to compute the risk scores for lncRNA data based on 12 lncRNA signatures
#' @import stringr
#' @param data data to be computed for cancer Prognosis risk scores; a data matrix, rows= lncRNA annotated with lncRNA names and columns=terms/samples.
#' @param annot annot holds RNA names in the data.
#' @param RNASignature RNA signatures, if none, will import the 12 lncRNA signatures.
#' @examples 
#' data(TCGA500)
#' data = exprs(TCGA500) 
#' annot = fData(TCGA500)
#' res =  lncRNA12model(data,annot)
#' @details
#' The 12 lncRNA signatures are "RP1-34M23.5","RP11-202K23.1", "RP11-560G2.1", "RP4-591L5.2","RP13-104F24.2", "RP11-506D12.5","ERVH48-1","RP4-613B23.1",  "RP11-360F5.1" ,"CTD-2031P19.5", "RP11-247A12.8", "SNHG7".
#' @return 
#' A Numeric Vector
#' @references 
#' \itemize{
#'  \item \strong{lncRNA12model}: Zhou M, Zhong L, Xu W, Sun Y, Zhang Z, Zhao H, Yang L, Sun J. Discovery of potential prognostic long non-coding RNA biomarkers for predicting the risk of tumor recurrence of breast cancer patients. Sci Rep. 2016;6:31038.
#' }
#' @export
#' 
lncRNA12model <- function(data,annot,RNASignature=NULL) #"TCGA500"
{
  if(is.null(RNASignature))
  {
    #data("sig.lncRNA12")
    lncRNA12 = CancerSubtypesPrognosis::lncRNA12
  }
  else
    lncRNA12=RNASignature
  tdata = NULL
  SigW = NULL
  for(j in 1:nrow(lncRNA12)) {
    index = which(annot$Gene.symbol %in% lncRNA12$Gene.symbol[j])
    if(length(index)!=0) {
      ## Merge the same probles of RNAs with average value
      if (length(index)>1) 
        tdata = cbind(tdata, apply(data[index,],2,mean))
      else
        tdata = cbind(tdata, data[index,])
      SigW = rbind(SigW, lncRNA12[j,])
    }
  }
  
  if(is.null(tdata) | is.null(SigW))
  {
    rs=NULL
  }
  else
  {
    colnames(tdata)=rownames(SigW)
    rs=tdata %*% SigW$weight
  }
  return(rs)
}

#' Evaluate Cancer Prognosis based on 5 lncRNA signatures
#' Cancer Prognosis to compute the risk scores for lncRNA data based on 5 lncRNA signatures
#' @import stringr
#' @param data data to be computed for cancer Prognosis risk scores; a data matrix, rows= lncRNA annotated with lncRNA names and columns=terms/samples.
#' @param annot annot holds RNA names in the data.
#' @param RNASignature RNA signatures, if none, will import the 5 lncRNA signatures.
#' @examples 
#' data(TCGA500)
#' data = exprs(TCGA500) 
#' annot = fData(TCGA500)
#' res =  lncRNA5model(data,annot)
#' @details
#' The 5 lncRNA signatures are "RP11-524D16-A.3", "HOTAIR", "AL645608.1", "TSPOAP1-AS1", "RP11-13L2.4".
#' @return 
#' A Numeric Vector
#' @references 
#' \itemize{
#'  \item \strong{lncRNA5model}: Li J, Wang W, Xia P, et al. Identification of a five‐lncRNA signature for predicting the risk of tumor recurrence in patients with breast cancer. Int J Cancer. 2018; 143: 2150‐2160.
#' }
#' @export
#' 
lncRNA5model <- function(data,annot,RNASignature=NULL)
{
  if(is.null(RNASignature))
  {
    #data("sig.lncRNA5")
    lncRNA5 = CancerSubtypesPrognosis::lncRNA5
  }
  else
    lncRNA5=RNASignature
  tdata = NULL
  SigW = NULL
  for(j in 1:nrow(lncRNA5)) {
    index = which(annot$Gene.symbol %in% lncRNA5$Gene.symbol[j])
    if(length(index)!=0) {
      ## Merge the same probles of RNAs with average value
      if (length(index)>1) 
        tdata = cbind(tdata, apply(data[index,],2,mean))
      else
        tdata = cbind(tdata, data[index,])
      SigW = rbind(SigW, lncRNA5[j,])
    }
  }
  
  if(is.null(tdata) | is.null(SigW))
  {
    rs=NULL
  }
  else
  {
    colnames(tdata)=rownames(SigW)
    rs=tdata %*% SigW$weight
  }
  return(rs)
}

#' Evaluate Cancer Prognosis based on 6 lncRNA signatures
#' Cancer Prognosis to compute the risk scores for lncRNA data based on 6 lncRNA signatures
#' @import stringr
#' @param data data to be computed for cancer Prognosis risk scores; a data matrix, rows= lncRNA annotated with lncRNA names and columns=terms/samples.
#' @param annot annot holds RNA names in the data.
#' @param RNASignature RNA signatures, if none, will import the 6 lncRNA signatures.
#' @examples 
#' data(TCGA500)
#' data = exprs(TCGA500) 
#' annot = fData(TCGA500)
#' res =  lncRNA6model(data,annot)
#' @details
#' The 6 lncRNA signatures are "HAGLR","STK4-AS1", "DLEU7-AS1", "LINC00957","LINC01614", "ITPR1-AS1".
#' @return 
#' A Numeric Vector
#' @references 
#' \itemize{
#'  \item \strong{lncRNA6model}: Zhong L, Lou G, Zhou X, Qin Y, Liu L, Jiang W. A six-long non-coding RNAs signature as a potential prognostic marker for survival prediction of ER-positive breast cancer patients. Oncotarget. 2017;8(40):67861.
#' }
#' @export
#' 
lncRNA6model <- function(data,annot,RNASignature=NULL)
{
  if(is.null(RNASignature))
  {
    #data("sig.lncRNA6")
    lncRNA6 = CancerSubtypesPrognosis::lncRNA6
  }
  else
    lncRNA6=RNASignature
  tdata = NULL
  SigW = NULL
  for(j in 1:nrow(lncRNA6)) {
    index = which(annot$Gene.symbol %in% lncRNA6$Gene.symbol[j])
    if(length(index)!=0) {
      ## Merge the same probles of RNAs with average value
      if (length(index)>1) 
        tdata = cbind(tdata, apply(data[index,],2,mean))
      else
        tdata = cbind(tdata, data[index,])
      SigW = rbind(SigW, lncRNA6[j,])
    }
  }
  
  if(is.null(tdata) | is.null(SigW))
  {
    rs=NULL
  }
  else
  {
    colnames(tdata)=rownames(SigW)
    rs=tdata %*% SigW$weight
  }
  return(rs)
}

