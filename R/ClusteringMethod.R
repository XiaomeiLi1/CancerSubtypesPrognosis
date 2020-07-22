#' Execute Consensus Clustering
#'
#' This function is based on the R package "ConsensusClusterPlus". 
#' We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The parameters are compatible to the original R package "ConsensusClusterPlus" function "ConsensusClusterPlus()".\cr
#' Please note: we add a new parameter "clusterNum" which represents the result with cancer subtypes group we want to return. 
#' 
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @param clusterNum A integer representing the return cluster number, this value should be less
#' than maxClusterNum(maxK). This is the only additional parameter in our function compared to the original
#' R package "ConsensusClusterPlus". All the other parameters are compatible to the function "ConsensusClusterPlus().
#' @param d data to be clustered; either a data matrix where columns=items/samples and rows are features. For example, a gene expression matrix of genes in rows and microarrays in columns, or ExpressionSet object, or a distance object (only for cases of no feature resampling)  
#' 
#' Please Note: We add a new data type (list) for this parameter. Please see details and examples.
#' @param maxK integer value. maximum cluster number  for Consensus Clustering Algorithm to evaluate.
#' @param reps  integer value. number of subsamples(in other words, The iteration number of each cluster number) 
#' @param clusterAlg character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
#' @param distance character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param  title character value for output directory. This title can be an absolute or relative path
#' @param pItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param pFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param plot Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param innerLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param finalLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param writeTable Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param verbose Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param corUse Please refer to the "ConsensusClusterPlus" package for detailed information.

#' @return A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{originalResult} : The clustering result of the original function "ConsensusClusterPlus()"
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  
#'  }
#'  
#' @details 
#'  If the data is a list containing the matched mutli-genomics  data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   we use "z-score" to normalize features for each data matrix first. Then all the normalized data matrices from the data list are concatenated
#'   according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.
#'   
#' @seealso \code{ConsensusClusterPlus}
#'   
#' @examples
#' ### The input dataset is a single gene expression matrix.
#' data(GeneExp)
#' data(miRNAExp)
#' result1=ExecuteCC(clusterNum=3,d=GeneExp,maxK=10,clusterAlg="hc",distance="pearson",title="GBM")
#' result1$group
#' 
#' ### The input dataset is multi-genomics data as a list
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result2=ExecuteCC(clusterNum=3,d=GBM,maxK=5,clusterAlg="hc",distance="pearson",title="GBM")
#' result2$group
#' 
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., Golub, T. (2003) Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Machine Learning, 52, 91-118.
#' @export
#'
ExecuteCC<-function(clusterNum,
                    d,maxK=10,clusterAlg="hc",
                    distance="pearson",title="ConsensusClusterResult",
                    reps=500, pItem=0.8, pFeature=1,plot="png",
                    innerLinkage="average", finalLinkage="average",
                    writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,
                    verbose=FALSE,corUse="everything")
{
  start = Sys.time()
  if(is.list(d))
  {
    temp=NULL
    for(i in 1: length(d))
    {
      temp=rbind(temp,d[[i]])
    }
    temp=t(scale(t(temp)))
  }
  else
   temp=d
  originalResult=ConsensusClusterPlus(
      temp, maxK=maxK,clusterAlg=clusterAlg,
      distance=distance,title=title,
      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
      innerLinkage=innerLinkage, finalLinkage=finalLinkage,
      writeTable=writeTable,weightsItem=weightsItem,weightsFeature=weightsFeature,
      verbose=verbose,corUse=corUse)
  group=originalResult[[clusterNum]][["consensusClass"]]
  distanceMatrix=originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix,'class')="Similarity"
  #icl=calcICL(result,title =fileName,plot="png" )
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=originalResult,timing=time.taken)
  result
}



#' Execute iCluster (Integrative clustering of multiple genomic data)
#'
#' Shen (2009) proposed a latent variable regression with a lasso constraint for joint modeling of multiple omics 
#' data types to identify common latent variables that can be used to cluster patient samples into biologically and clinically relevant disease subtypes.\cr
#' This function is based on the R package "iCluster". 
#' The R package "iCluster" should be installed. 
#' We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The parameters is compatible to the original R package "iCluster" function "iCluster2()".\cr
#' Please note: The data matrices are transposed in our function comparing to the original R package "iCluster" on the behalf of the unified input format with other functions.
#' We try to build a standardized flow for cancer subtypes analysis and validation.
#'
#' @importFrom iCluster iCluster2 plotiCluster
#' @param datasets A list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' In order to unify the input parameter with other clustering methods, the data matrices are transposed comparing to the definition in the original "iCluster" package.
#' @param k Number of subtypes for the samples
#' @param lambda Penalty term for the coefficient matrix of the iCluster model
#' @param scalar Logical value. If true, a degenerate version assuming scalar covariance matrix is used.
#' @param max.iter  maximum iteration for the EM algorithm
#' @param scale Logical value. If true, the genomic features in the matrix is centered.

#' @return A list with the following elements.
#'\itemize{
#'   \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{originalResult} : The clustering result of the original function "iCluster2()".
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  
#'  }
#'  
#'  
#' @details 
#'  For iCluster algorithm, it cannot process high-dimensional data, otherwise it is very very time-consuming or reports a mistake.  
#'  Based on test, it could smoothly run for the matrix with around 1500 features. Normally it need feature selection step first to reduce feature number.
#' @references
#' Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. Bioinformatics 25, 2906-2912.\cr
#' Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. PLoS ONE 7, e35236
#' 
#' @seealso \code{\link{iCluster2}}
#' 
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' data1=FSbyVar(GeneExp, cut.type="topk",value=500)
#' data2=FSbyVar(miRNAExp, cut.type="topk",value=100)
#' GBM=list(GeneExp=data1,miRNAExp=data2)
#' result=ExecuteiCluster(datasets=GBM, k=3, lambda=list(0.44,0.33,0.28))
#' result$group
#' @export
#'
ExecuteiCluster<-function(datasets, k, lambda=NULL, scale=TRUE, scalar=FALSE, max.iter=10)
{
  start = Sys.time()
  data1=list()
  for(i in 1:length(datasets))
  {
    data1[[i]]=t(datasets[[i]])
  }
  
  fit=iCluster2(datasets=data1, k=k, lambda=lambda, scale=scale, scalar=scalar, max.iter=10) 
  
#  plotiCluster(fit=fit, label=rownames(data1[[1]]))
  plotiCluster(fit=fit, label=NULL)
  group=fit$clusters
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,originalResult=fit,timing=time.taken)
  result
}


#' Execute SNF(Similarity Network Fusion )
#'
#' SNF is a multi-omics data processing method that constructs a fusion patient similarity network
#' by integrating the patient similarity obtained from each of the genomic data types. 
#' SNF calculates the similarity between patients using each single data type separately. The similarities
#' between patients from different data types are then integrated by a cross-network diffusion process to construct the fusion patient similarity matrix. 
#' Finally, a clustering method is applied to the fusion patient similarity matrix to cluster patients into different groups, which imply different cancer subtypes.
#' This function is based on the R package "SNFtool". 
#' The R package "SNFtool" should be installed. 
#' We write a function to integrate the clustering process and unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation.\cr
#' Please note: The data matrices are transposed in our function comparing to the original R package "SNFtools".
#' We try to build a standardized flow for cancer subtypes analysis and validation.
#'
#' @importFrom SNFtool dist2 affinityMatrix SNF spectralClustering displayClusters
#' @param datasets A list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' @param clusterNum A integer representing the return cluster number
#' @param K Number of nearest neighbors
#' @param alpha Variance for local model
#' @param t Number of iterations for the diffusion process
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#' @return A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{originalResult} : The clustering result of the original SNF algorithm"
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  }
#' @references
#' B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014
#' 
#' @seealso \code{\link{affinityMatrix}} \code{\link{SNF}}
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#' result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#' result$group
#' @export
#'
ExecuteSNF<-function(datasets, clusterNum, K=20, alpha=0.5, t=20,plot=TRUE)
{
  start = Sys.time()
  # library("SNFtool")
  #browser()
  W_temp=list()
  for(i in 1:length(datasets))
  {
    distance=dist2(as.matrix(t(datasets[[i]])), as.matrix(t(datasets[[i]])))
    W_temp[[i]] = affinityMatrix(distance, K, alpha)
  }
  if(length(W_temp) == 1) {
    W = W_temp[[1]]
  }
  else W = SNFtool::SNF(W_temp, K=K, t=t)
  group =spectralClustering(W,clusterNum)
  
  diag(W)=0
  diag(W)=max(W)
  distanceMatrix=W
  attr(distanceMatrix,'class')="Similarity"
  
  if(plot)
    displayClusters(W, group)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=group,timing=time.taken)
  result
}


#' Execute the combined SNF (Similarity Network Fusion) and Consensus clustering
#'
#' This function is a combined process of SNF and Consensus Clustering for cancer subtypes identification.
#' First it applied SNF to get the fusion patients similarity matrix. Then use this 
#' fusion patients similarity matrix as the sample distance for Consensus Clustering.
#'
#' @importFrom SNFtool dist2 affinityMatrix SNF 
#' @param datasets A list containing data matrices. For each data matrix, 
#' the rows represent genomic features, and the columns represent samples. Same as ExecuteSNF
#' @param clusterNum A integer representing the return cluster number. Same as ExecuteSNF
#' @param K Number of nearest neighbors. Same as ExecuteSNF
#' @param alpha Variance for local model. Same as ExecuteSNF
#' @param t Number of iterations for the diffusion process. Same as ExecuteSNF
#' @param maxK integer value. maximum cluster number for Consensus Clustering Algorithm to evaluate. Same as ExecuteCC.
#' @param pItem Same as ExecuteCC
#' @param reps integer value. number of subsamples(in other words, The iteration number of each cluster number). Same as ExecuteCC
#' @param title character value for output directory. This title can be an absolute or relative path. Same as ExecuteCC
#' @param plot Same as ExecuteCC
#' @param finalLinkage Same as ExecuteCC
#' 
#' @return Same as the ExecuteCC(). A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{originalResult} : The clustering result of the original function "ConsensusClusterPlus()"
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @seealso \code{\link{ExecuteSNF}} \code{\link{ExecuteCC}}
#' @examples
#' 
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp,miRNAExp)
#' result=ExecuteSNF.CC(GBM, clusterNum=3, K=20, alpha=0.5, t=20,
#'                     maxK = 5, pItem = 0.8,reps=500, 
#'                     title = "GBM", plot = "png", 
#'                     finalLinkage ="average")
#' result$group
#' @export
#'
ExecuteSNF.CC<-function(datasets, clusterNum, K=20, alpha=0.5, t=20,
                        maxK = 10, pItem = 0.8,reps=500, 
                        title = "ConsensusClusterResult", plot = "png", 
                        finalLinkage ="average")
{
  start = Sys.time()
  W_temp=list()
  for(i in 1:length(datasets))
  {
    distance=dist2(as.matrix(t(datasets[[i]])), as.matrix(t(datasets[[i]])))
    W_temp[[i]] = affinityMatrix(distance, K, alpha)
  }
  if(length(W_temp) == 1) {
    W = W_temp[[1]]
  }
  else W = SNFtool::SNF(W_temp, K=K, t=t)
  W=as.dist(W)
  result=ExecuteCC(clusterNum=clusterNum,d=W,maxK=maxK,
                   clusterAlg="spectralAlg",title=title,reps=reps,
                   pItem=pItem,plot=plot,
                   finalLinkage=finalLinkage)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result$timing = time.taken
  result
}


#' Execute the WSNF(Weighted Similarity Network Fusion)
#'
#' WSNF is a caner subtype identificaton method with the assistance of the gene regulatory network information. The basic idea of the WSNF is
#' to set the different regulatory importance(ranking) for each feature. In the WSNF manuscript, WSNF makes use of the miRNA-TF-mRNA regulatory 
#' network to take the importance of the features into consideration.
#' 
#' @importFrom SNFtool affinityMatrix SNF 
#' @param datasets A list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' @param feature_ranking A list containing numeric vetors. The length of the feature_ranking list should equal to the length of datasets list.
#' For each numeric vetor represents the ranking of each feature in the corresponding data matrix. The order of the ranking should also mathch 
#' the order of the features in the corresponding data matrix.  
#' We proive a ranking list for most mRNA, TF(transcription factor) and miRNA features. The ranking for features caculated based on the miRNA-TF-miRNA 
#' regulatory network which was promoted in our published work: Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and 
#' Expression Data(PLos One,2016).
#' @param beta A tuning parameter for the feature_ranking contributes the weight of each feature. \cr
#' A linear model is applied to integrate feature_ranking and MAD(median absolute deviation) to generated the final weight for each feature using 
#' for the algorithm. The final weight is cauculated as the formula below:\cr
#' Weight(f_i)=beta * feature_ranking + (1-beta) MAD(f_i)
#' @param clusterNum A integer representing the return cluster number
#' @param K Number of nearest neighbors
#' @param alpha Variance for local model
#' @param t Number of iterations for the diffusion process
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#' 
#' @return A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{originalResult} : The clustering result of the original SNF algorithm.
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @seealso \code{\link{ExecuteSNF}}
#' @references 
#' Xu, T., Le, T. D., Liu, L., Wang, R., Sun, B., & Li, J. (2016). Identifying cancer subtypes from mirna-tf-mrna regulatory networks and expression data. PloS one, 11(4), e0152792.
#' @examples
#' data(GeneExp)
#' data(miRNAExp)
#' GBM=list(GeneExp,miRNAExp)
#' ###1. Use the defualt ranking in the package.
#' data(Ranking)
#' ####Retrieve the feature ranking for genes
#' gene_Name=rownames(GeneExp)
#' index1=match(gene_Name,Ranking$mRNA_TF_miRNA.v21_SYMBOL)
#' gene_ranking=data.frame(gene_Name,Ranking[index1,],stringsAsFactors=FALSE)
#' index2=which(is.na(gene_ranking$ranking_default))
#' gene_ranking$ranking_default[index2]=min(gene_ranking$ranking_default,na.rm =TRUE)
#' 
#' ####Retrieve the feature ranking for miRNAs
#' miRNA_ID=rownames(miRNAExp)
#' index3=match(miRNA_ID,Ranking$mRNA_TF_miRNA_ID)
#' miRNA_ranking=data.frame(miRNA_ID,Ranking[index3,],stringsAsFactors=FALSE)
#' index4=which(is.na(miRNA_ranking$ranking_default))
#' miRNA_ranking$ranking_default[index4]=min(miRNA_ranking$ranking_default,na.rm =TRUE)
#' ###Clustering
#' ranking1=list(gene_ranking$ranking_default ,miRNA_ranking$ranking_default)
#' result1=ExecuteWSNF(datasets=GBM, feature_ranking=ranking1, beta = 0.8, clusterNum=3, 
#'                    K = 20,alpha = 0.5, t = 20, plot = TRUE)
#' 
#' ###2. User input ranking
#' # Fabricate a ranking list for demonstrating the examples.
#' ranking2=list(runif(nrow(GeneExp), min=0, max=1),runif(nrow(miRNAExp), min=0, max=1))
#' result2=ExecuteWSNF(datasets=GBM, feature_ranking=ranking2, beta = 0.8, clusterNum=3, 
#'                    K = 20,alpha = 0.5, t = 20, plot = TRUE)
#' 
#' @export
#'
ExecuteWSNF<-function(datasets,feature_ranking,beta=0.8,clusterNum,K=20, alpha=0.5, t=20,plot=TRUE)
{
  start = Sys.time()
  if(is.list(feature_ranking))
  {
    if(length(feature_ranking)==length(datasets))
    {
      W_temp=list()
      for(i in 1:length(datasets))
      {
        mads=apply(datasets[[i]],1,mad)
        mads=mads/sum(mads)
        feature_ranking1=as.numeric(feature_ranking[[i]])
        feature_ranking1=feature_ranking1/sum(feature_ranking1)
        weight=feature_ranking1*beta+(1-beta)*mads
        distance=.distanceWeighted2(as.matrix(t(datasets[[i]])),weight)
        W_temp[[i]] = affinityMatrix(distance, K, alpha)
      }
      if(length(W_temp) == 1) {
        W = W_temp[[1]]
      }
      else W = SNFtool::SNF(W_temp, K=K, t=t)
      group =spectralClustering(W,clusterNum)
      
      diag(W)=0
      diag(W)=max(W)
      distanceMatrix=W
      attr(distanceMatrix,'class')="Similarity"
      
      if(plot)
        displayClusters(W, group)
      time.taken = as.numeric(Sys.time() - start, units='secs')
      result=list(group=group,distanceMatrix=distanceMatrix,originalResult=group,timing=time.taken)
    }
    else
    {
      stop("The length of feature_ranking(list) is not equal to lengh of datasets(list)")
    }
  }
result
}



#' Execute Consensus NMF (Nonnegative matrix factorization)
#'
#' Brunet applied nonnegative matrix factorization (NMF) to analyze the Gene MicroArray dataset in 2004. In the original paper, the author
#' proved that NMF is an efficient method for distinct molecular patterns identification and provides a powerful method 
#' for class discovery. This method was implemented in an R package "NMF". Here we applied the "NMF" package to 
#' conduct the cancer subtypes identification. We write a shell to unify the input and output format.
#' It is helpful for the standardized flow of cancer subtypes analysis and validation. 
#' The R package "NMF" should be installed. 
#' 
#' @importFrom NMF nmf predict
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' If the matrices have negative values, first the negative values will be set to zero to get a matrix 1;
#' all the positive values will be set to zero to get the matrix 2; then a new matrix with all positive values will be
#' get by concatenating matrix1 and -maxtrix2.
#'  
#' 
#' @param clusterNum Number of subtypes for the samples
#' @param nrun Number of runs to perform NMF. A default of 30 runs are performed, allowing the computation of a consensus matrix that is used in selecting the best result for cancer subtypes identification as Consensus Clustering method.
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{originalResult} : A NMFfitX class from the result of function "nmf()".
#'  
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' [1] Brunet, Jean-Philippe, Pablo Tamayo, Todd R Golub, and Jill P Mesirov. "Metagenes and Molecular Pattern Discovery Using Matrix Factorization." Proceedings of the National Academy of Sciences 101, no. 12 (2004):4164-69.
#' 
#' [2] Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for Nonnegative Matrix Factorization." BMC Bioinformatics 11 (2010): 367. doi:10.1186/1471-2105-11-367.
#' @seealso \code{\link{nmf}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' #To save the  execution time, the nrun is set to 5, but the recommended value is 30.
#' result=ExecuteCNMF(GeneExp,clusterNum=3,nrun=5)
#' result$group
#' 
#' @export
#'
ExecuteCNMF<-function(datasets, clusterNum,nrun=30 )
{
  start = Sys.time()
  if(is.list(datasets))
  {
    temp=NULL
    for(i in 1: length(datasets))
    {
      temp=rbind(temp,datasets[[i]])
    }
  }
  else
    temp=datasets
  ## change all value to positive
  data1=rbind(pmax(temp,0),-pmin(temp,0))
  index=which(rowSums(data1)==0)
  data1=data1[-index,]
  res=nmf(data1,rank=clusterNum,nrun=nrun)
  
  distanceMatrix=slot(res,"consensus")
  attr(distanceMatrix,'class')="Similarity"
  
  group=as.numeric(as.vector(predict(res)))
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=res,timing=time.taken)
  result
}

#' Execute PAM50 classifier
#'
#' PAM50 is a gene-based method to classify samples to the five subtypes: Basal, Luminal A, Luminal B, Her2-enriched and Normal-like. 
#' PAM50 constructs a centroid-based predictor by using the Prediction Analysis of Microarray (PAM) algorithm on 50 gene signatures.
#' The R package "genefu" should be installed. 
#' 
#' @importFrom genefu molecular.subtyping
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#'  
#' 
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' Joel S Parker, Michael Mullins, Maggie CU Cheang, Samuel Leung, David Voduc, Tammi Vickery, Sherri Davies, Christiane Fauron, Xiaping He, Zhiyuan Hu, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. Journal of clinical oncology, 27(8):1160, 2009.
#' 
#' @seealso \code{\link{molecular.subtyping}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' result=ExecutePAM50(GeneExp)
#' result$group
#' 
#' @export
#'
ExecutePAM50<-function(datasets)
{
  start = Sys.time()
  if(is.list(datasets))
  {
    temp=NULL
    for(i in 1: length(datasets))
    {
      temp=rbind(temp,datasets[[i]])
    }
  }
  else
    temp=datasets
  
  mRNA_Anno = data.frame(Gene.Symbol=rownames(datasets))

  res <- molecular.subtyping(sbt.model = "pam50",data=t(datasets),
                                 annot=mRNA_Anno,do.mapping=FALSE)
  group = as.numeric(res$subtype)

  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,timing=time.taken)
  result
}

#' Execute IntClust
#' IntClust is a integrative method to classify samples to ten breast cancer subtypes. 
#' IntClust applies iCluster on a matched mRNA-CNV breast cancer dataset with 997 samples and identifies ten breast cancer subtypes (so-called integrative subtypes). Similar to PAM50, IntClust builds three centroid-based predictors based on 612 cis-eQTLs gene drivers by using PAM.
#' The R package "genefu" should be installed. 
#' 
#' @importFrom genefu molecular.subtyping
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#'  
#' 
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' [1] Christina Curtis, Sohrab P Shah, Suet-Feung Chin, Gulisa Turashvili, Oscar M Rueda, Mark J Dunning, Doug Speed, Andy G Lynch, Shamith Samarajiwa, Yinyin Yuan, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature, 486(7403):346, 2012.
#' 
#' [2] H Raza Ali, Oscar M Rueda, Suet-Feung Chin, Christina Curtis, Mark J Dunning, Samuel AJR Aparicio, and Carlos Caldas. Genome-driven integrated classification of breast cancer validated in over 7,500 samples. Genome biology, 15(8):431, 2014.
#' @seealso \code{\link{molecular.subtyping}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' result=ExecuteIntClust(GeneExp)
#' result$group
#' 
#' @export
#'
###IC10 classifier for copy number data and/or gene expression data
ExecuteIntClust<-function(datasets)
{
  start = Sys.time()
  if(is.list(datasets))
  {
    temp=NULL
    for(i in 1: length(datasets))
    {
      temp=rbind(temp,datasets[[i]])
    }
  }
  else
    temp=datasets
  
  mRNA_Anno = data.frame(Gene.Symbol=rownames(datasets))
  
  res = molecular.subtyping(sbt.model = "intClust",data=t(datasets),
                      annot=mRNA_Anno,do.mapping=FALSE)
  
  group = as.numeric(res$subtype)
  
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,timing=time.taken)
  result
}

#' Execute CIMLR (Cancer Integration via Multikernel Learning)
#'
#' CIMLR calculates the similarity between patients in multi-omic data by combining a set of Gaussian kernels for each single-omic data.
#' 
#' The R package "CIMLR" should be installed. 
#' 
#' @importFrom CIMLR CIMLR
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#'  
#' 
#' @param clusterNum Number of subtypes for the samples
#' @param k tuning parameter in CIMLR. A default of 10 is performed.
#' @param ncore ratio of the number of cores to be used when computing the multi-kernel
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' Daniele Ramazzotti, Avantika Lal, Bo Wang, Serafim Batzoglou, and Arend Sidow. Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. Nature communications, 9(1):4453, 2018.
#' @seealso \code{\link{CIMLR}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' result=ExecuteCIMLR(GeneExp,clusterNum=5)
#' result$group
#' 
#' @export
#'
ExecuteCIMLR<-function(datasets, clusterNum, k = 10, ncore = 0, plot=TRUE)
{
  start = Sys.time()
  if(is.list(datasets))
  {
    res = CIMLR(X = datasets, c = clusterNum, k, cores.ratio = ncore)
  }
  else
    res = CIMLR(X = list(datasets), c = clusterNum, k, cores.ratio = ncore)
  
  group =as.numeric(res$y$cluster)
  
  distanceMatrix=res$S
  attr(distanceMatrix,'class')="Similarity"
  
  if(plot)
    displayClusters(distanceMatrix, group)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,distanceMatrix=distanceMatrix,timing=time.taken)
  result
}

#' Execute PINS (Perturbation clustering for data INtegration and disease Subtyping)
#'
#' PINS initially clusters patients based on each omic data separately using perturbation clustering and outputs the optimal number of clusters ki, original connectivity matrix Ci and the perturbed connectivity matrix Ai, where i is the index of i-th omic data.
#' The R package "PINSPlus" should be installed. 
#' 
#' @importFrom PINSPlus PerturbationClustering
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#' 
#' @param clusterNum Number of subtypes for the samples
#' @param ncore Number of cores that the algorithm should use. Default value is 1.
#' @param clusteringMethod The name of built-in clustering algorithm that PerturbationClustering will use. Currently supported algorithm are kmeans, pam and hclust. Default value is "kmeans".
#' @param perturbMethod The clustering algorithm function that will be used instead of built-in algorithms.
#' @param iterMin The minimum number of iterations. Default value is 20.
#' @param iterMax The maximum number of iterations. Default value is 200.
#' @param madMin The minimum of Mean Absolute Deviation of AUC of Connectivity matrix for each k. Default value is 1e-03.
#' @param msdMin The minimum of Mean Square Deviation of AUC of Connectivity matrix for each k. Default value is 1e-06.
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' Tin Nguyen, Rebecca Tagett, Diana Diaz, and Sorin Draghici. A novel approach for data integration and disease subtyping. Genome research, 27(12):2025-2039, 2017.
#' @seealso \code{\link{PerturbationClustering}}
#' 
#' 
#' @examples
#' data(GeneExp)
#' result=ExecutePINS(GeneExp,clusterNum=5)
#' result$group
#' 
#' @export
#'
ExecutePINS <- function(datasets, clusterNum, ncore = 1, clusteringMethod = "kmeans", 
                        perturbMethod = "noise",iterMin = 20, iterMax = 200, 
                        madMin = 0.001, msdMin = 1e-06) 
{
  start = Sys.time()
  
  if (!is.list(datasets)) {
    pins.ret = PINSPlus::PerturbationClustering(data=t(datasets),
                                                kMax = clusterNum, ncore = ncore, 
                                                clusteringMethod = clusteringMethod,  
                                                perturbMethod = perturbMethod,
                                                iterMin = iterMin, iterMax = iterMax, 
                                                madMin = madMin, msdMin = msdMin)
    group = pins.ret$cluster
    
  } 
  else {
    omics.transposed = lapply(datasets, t)
    pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed,
                                            kMax = clusterNum, ncore = ncore, 
                                            clusteringMethod = clusteringMethod,  
                                            perturbMethod = perturbMethod,
                                            iterMin = iterMin, iterMax = iterMax, 
                                            madMin = madMin, msdMin = msdMin)
    group = as.numeric(as.factor(pins.ret$cluster2))
  }
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list("group"=group, timing=time.taken))
}

#' Execute NEMO (NEighborhood based Multi-Omics clustering)
#'
#' NEMO is a similarity-based multi-omic clustering method based on the radial basis of function kernel and spectral clustering method. 
#' Different from other multi-omic clustering methods, NEMO can apply to partial data that some omics may not be measured for some patients. NEMO is simple, and faster than other multi-omics clustering methods but achieved the comparable performance to others. 
#' 
#' @param datasets A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.
#'  
#' 
#' @param clusterNum Number of subtypes for the samples
#' @param num.neighbors The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#' @return 
#' A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'   
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group 
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'   
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'   
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'   
#'  \item \strong{timing} : The running time.
#'  }
#'  
#' @details
#'  If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.

#' @references
#' Nimrod Rappoport and Ron Shamir. Nemo: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356, 2019.
#' @examples
#' data(GeneExp)
#' result=ExecuteNEMO(GeneExp,clusterNum=3)
#' result$group
#' 
#' @export
#'
ExecuteNEMO<-function(datasets, clusterNum, num.neighbors=50,plot=TRUE)
{
  start = Sys.time()
  if(is.list(datasets))
  {
    res = nemo.clustering(datasets, num.clusters=clusterNum, num.neighbors=num.neighbors)
  }
  else
    res = nemo.clustering(list(datasets), num.clusters=clusterNum, num.neighbors=num.neighbors)

  group =as.numeric(res$clustering)
  
  distanceMatrix=res$graph
  attr(distanceMatrix,'class')="Similarity"
  
  if(plot)
    displayClusters(distanceMatrix, group)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(group=group,distanceMatrix=distanceMatrix,timing=time.taken)
  result
}
