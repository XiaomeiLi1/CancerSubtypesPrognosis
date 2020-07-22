#' Dataset: Gene expression
#'
#'A glioblastoma (GBM) gene expression dataset downloaded from TCGA. This is a small dataset with 1500 genes 
#'and  100 cancer samples extracted from gene expression data for examples.
#' 
#'
#'\itemize{
#'  \item Rows are genes
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name GeneExp
#'@examples
#' data(GeneExp)
NULL


#' Dataset: miRNA expression
#'
#'A glioblastoma (GBM) miRNA expression dataset downloaded from TCGA. This is a small miRNA expression dataset with 470 miRNAs and 
#' 100 cancer samples extracted from miRNA expression data for examples.
#'
#'\itemize{
#'  \item Rows are miRNAs
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name miRNAExp
#'@examples
#' data(miRNAExp)
NULL


#' Dataset: Survival time
#'
#'\itemize{
#'  \item A vector representing the right censored survival time (days) for GBM cancer patients matched with the "GeneExp" and "miRNAExp"
#' datasets.
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A numeric vector
#'@name time
#'@examples
#' data(time)
NULL

#' Dataset: Survival status
#'
#'\itemize{
#'  \item A vector representing the survival status for GBM cancer patients matched with the "GeneExp" and "miRNAExp"
#'. 0=alive or censored, 1=dead
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A numeric vector
#'@name status
#'@examples
#' data(status)
NULL

#' Dataset: A default ranking of features for the fuction ExecuteWSNF()
#' 
#' A dataframe represents the regulatory ranking for features(mRNA,TF,miRNA) caculated based on the miRNA-TF-miRNA regulatory network which was promoted in our published work: 
#' Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data(PLos One,2016).
#'
#'\itemize{
#'  \item  mRNA_TF_miRNA_ID : ENTREZID for genes(mRNA,TF)  and miRBase Accession ID for miRNAs.
#'  \item  mRNA_TF_miRNA.v21._SYMBOL: gene symbol and miRNA names(miRBase Version 21)
#'  \item  feature_ranking:  the numeric values represents regulatory ranking for each feature.
#'}
#'
#'@docType data
#'@keywords datasets
#'@format dataframe
#'@name Ranking
#'@references 
#' Xu, T., Le, T. D., Liu, L., Wang, R., Sun, B., & Li, J. (2016). Identifying cancer subtypes from mirna-tf-mrna regulatory networks and expression data. PloS one, 11(4), e0152792.
#'@examples
#' data(Ranking)
NULL

#' Dataset: TCGA500 ExpressionSet
#'
#'The TCGA breast cancer datasets is downloaded from the TCGA data portal and the datasets consist of level 3 mRNA and miRNA, lncRNA expression data from multiple platforms.
#' 
#'
#'\itemize{
#'  \item ExpressionSet
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A ExpressionSet
#'@name TCGA500
#'@examples
#' data(TCGA500)
NULL

#' Dataset: RNA signatures
#'
#' The 7 miRNAs signatures are  "hsa-miR-103", "hsa-miR-1307", "hsa-miR-148b", "hsa-miR-328", "hsa-miR-484", "hsa-miR-874","hsa-miR-93".
#' The 30 mRNA signatures are "ACSL1" ,"ADAT1", "ANKRD52","BIRC6", "CPT1A", "CXCR7","DAAM1", "DIP2B", "FAM199X", "FAM91A1", "FRZB",  "GLA",  "GMCL1", "HRASLS",  "HSP90AA1", "MCM10", "ME1",  "NDRG1", "NOTCH2NL","OTUD6B","PDSS2", "PGK1",  "PIK3CA", "PTAR1", "SMG1",  "TRIM23", "TTC3",  "UBR5",  "UBXN7", "ZFC3H1".
#' Note that the miRbase version of miRNAs signatures is v16, please check the miRbase version of your own dataset.
#'
#'\itemize{
#'  \item Rows are RNAs
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data frame
#'@name RNASig
#'@examples
#' data(sig.RNA37)
NULL

#' Dataset: miRNA signatures
#'
#' The 10 miRNAs signatures are  "hsa-miR-144", "hsa-miR-150", "hsa-miR-210", "hsa-miR-27b", "hsa-miR-30c", "hsa-miR-342"    "hsa-miR-128a", "hsa-miR-135a", "hsa-miR-767-3p", "hsa-miR-769-3p".
#' Note that the miRbase version of miRNAs signatures is v9_2, please check the miRbase version of your own dataset.
#'
#'\itemize{
#'  \item Rows are miRNAs
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data frame
#'@name miRNA10
#'@examples
#' data(sig.miRNA10)
NULL

#' Dataset: lncRNA signatures
#'
#'12 lncRNA signatures
#'
#'\itemize{
#'  \item Rows are lncRNAs
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data frame
#'@name lncRNA12
#'@examples
#' data(sig.lncRNA12)
NULL

#' Dataset: lncRNA signatures
#'
#'5 lncRNA signatures
#'
#'\itemize{
#'  \item Rows are lncRNAs
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data frame
#'@name lncRNA5
#'@examples
#' data(sig.lncRNA5)
NULL

#' Dataset: lncRNA signatures
#'
#'6 lncRNA signatures
#'
#'\itemize{
#'  \item Rows are lncRNAs
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data frame
#'@name lncRNA6
#'@examples
#' data(sig.lncRNA6)
NULL