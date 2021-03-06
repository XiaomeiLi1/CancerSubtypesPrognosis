% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusteringMethod.R
\name{ExecutePAM50}
\alias{ExecutePAM50}
\title{Execute PAM50 classifier}
\usage{
ExecutePAM50(datasets)
}
\arguments{
\item{datasets}{A data matrix or a list containing data matrices. For each data matrix, the rows represent genomic features, and the columns represent samples.}
}
\value{
A list with the following elements.
\itemize{
 \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
  
  This is the most important result for all clustering methods, so we place it as the first component. The format of group 
  is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
  of group is also compatible with the K-means result and the hclust (after using the cutree() function).
  
 \item \strong{timing} : The running time.
 }
}
\description{
PAM50 is a gene-based method to classify samples to the five subtypes: Basal, Luminal A, Luminal B, Her2-enriched and Normal-like. 
PAM50 constructs a centroid-based predictor by using the Prediction Analysis of Microarray (PAM) algorithm on 50 gene signatures.
The R package "genefu" should be installed.
}
\details{
If the data is a list containing the matched mutli-genomics data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
  The data matrices in the list are concatenated according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
  Our purpose is to make convenient comparing the different method with same dataset format. See examples.
}
\examples{
data(GeneExp)
result=ExecutePAM50(GeneExp)
result$group

}
\references{
Joel S Parker, Michael Mullins, Maggie CU Cheang, Samuel Leung, David Voduc, Tammi Vickery, Sherri Davies, Christiane Fauron, Xiaping He, Zhiyuan Hu, et al. Supervised risk predictor of breast cancer based on intrinsic subtypes. Journal of clinical oncology, 27(8):1160, 2009.
}
\seealso{
\code{\link{molecular.subtyping}}
}
