#' Evaluate Cancer CancerSubtyping methods of mRNA, miRNA or multiomics data based on the running time, the average Silhouette score, and the p-value of Logrank test
#' 
#' CancerSubtypes to compute the running time, the average Silhouette score, and the p-value of Logrank test for mRNA, miRNA or multiomics data based on existing 11 methods: 
#'            "PAM50","IntClust","CC","CNMF","iCluster","SNF","SNF-CC","WSNF","CIMLR","PINS","NEMO"
#' @import genefu   
#' @import  ConsensusClusterPlus
#' @import iCluster
#' @import SNFtool
#' @import NMF
#' @import doParallel
#' @import foreach
#' @import CIMLR
#' @import PINSPlus
#' @importFrom  iC10 matchFeatures
#' @param dn datasets to be computed for cancer subtypes; a string vector. 
#' @param omics The types of data, "mRNA" or "miRNA" or "multiomics". 
#' @param methods A set of methods to be performed in the 11 methods: "PAM50","IntClust","CC","CNMF","iCluster","SNF","SNF-CC","WSNF","CIMLR","PINS","NEMO".
#' @param fileFolder A file folder name provided for saving results. 
#' @param logFile A file name provided for saving log information. 
#' 
#' @return 
#' a list contains timeTable matrix, silTable matrix, and pvalueTable matrix
#' 
#' @examples 
#' \dontrun{
#' dn = c("TCGA","UK","HEL","GSE19783")
#' methods = c("PAM50","IntClust","CC","CNMF","iCluster","SNF","SNF-CC","WSNF","CIMLR","PINS","NEMO")
#' omics = "mRNA"
#' res = CancerSubtypes(dn, omics, methods, fileFolder=omics, logFile = omics)
#' }
#'
#' @export
#' 
CancerSubtypes <- function(dn = NULL, omics = NULL, methods = NULL, 
                         fileFolder=NULL, logFile = NULL) 
{
  #### define log file
  fpath = getFilePath(fileFolder,logFile)
  sink(paste0(fpath,".txt"),append = TRUE,split = TRUE)
  
  ####load expression data
  loadData("data",dn)
  
  ####define number of cores
  ncore = detectCores()
  
  ###summary of results
  ###run time
  timeTable = matrix(data = NA, nrow = length(dn), ncol = length(methods))
  rownames(timeTable) = dn
  colnames(timeTable) = methods
  
  ###average Silhouette
  silTable = matrix(data = NA, nrow = length(dn), ncol = length(methods))
  rownames(silTable) = dn
  colnames(silTable) = methods
  
  ###p value of Logrank test
  pvalueTable = matrix(data = NA, nrow = length(dn), ncol = length(methods))
  rownames(pvalueTable) = dn
  colnames(pvalueTable) = methods
  
  #### get the expression data
  data = list()
  data1 = list()
  data2 = list()
  for(i in dn) {
    if(omics == "multiomics") {
      dn1 = paste(i,"mRNA",sep='_')
      dn2 = paste(i,"miRNA",sep='_')
      data = list("mRNA" = get(dn1),"miRNA" = get(dn2))
      
      ## feature selction
      data2$miRNA=FSbyMAD(data$miRNA, cut.type = "cutoff", 0.001)
      data2$mRNA=FSbyMAD(data$mRNA, cut.type = "topk", 2000)
      ##normalized multiomics data
      data1$mRNA=data.normalization(data2$mRNA, type = "feature_zscore")
      data1$miRNA=data.normalization(data2$miRNA, type = "feature_zscore")
    }
    else if(omics == "mRNA") {
      dn1 = paste(i,omics,sep='_')
      data = list("mRNA" =  get(dn1))
      data1$mRNA=FSbyMAD(data$mRNA, cut.type = "topk", 2000)
      data2$mRNA = data1$mRNA
    }
    else if(omics == "miRNA") {
      dn1 = paste(i,omics,sep='_')
      data = list("miRNA" =  get(dn1))
      data1$miRNA=FSbyMAD(data$miRNA, cut.type = "cutoff", 0.01)
      data2$miRNA = data1$miRNA
    }
    else {
      cat("Need to provide gene expression data.\n")
    }
    
    ###clinical information, first two columns are time and status
    dclin = get(paste(i,"clinical",sep='_'))
    
    ###PAM50
    if ("PAM50" %in% methods) {
      cat("run PAM50...\n")
      if (omics == "mRNA") {
        #### load PAM50 data set
        #### gene symbols
        data("pam50")
        gene.symbol = rownames(pam50$centroids)
        index = which(rownames(data$mRNA) %in% gene.symbol)
        temp=data$mRNA[index,]
        
        #####calculate similarity matrix
        ddist <- as.matrix(dist2(as.matrix(t(temp)), as.matrix((t(temp)))))
        distanceMatrix = affinityMatrix(ddist, 20, 0.5)
        attr(distanceMatrix,'class')="Similarity"
        
        PAM50_res = ExecutePAM50(data$mRNA)
        fpath = getFilePath(fileFolder, paste(i,"PAM50",sep = "_"))
        save(PAM50_res, file = fpath)
        
        timeTable[i,"PAM50"] = PAM50_res$timing
        silTable[i,"PAM50"] = getMeanSilhouette(PAM50_res$group, distanceMatrix)
        pvalueTable[i,"PAM50"] = getPvalue(dclin[,1],dclin[,2],PAM50_res$group)
        
        cat("finish PAM50.\n")
      }
      else {
        cat("Need to provide at least one source (mRNA or miRNA) of data.\n")
      }
    }
    
    ###IntClust
    if ("IntClust" %in% methods) {
      cat("run IntClust...\n")
      if (omics == "mRNA") {
        #### load iC50 library
        feat <- matchFeatures(Exp = data$mRNA, Exp.by.feat = "Gene.Symbol")
        #### gene symbols
        gene.symbol = rownames(feat$Exp)
        index = which(rownames(data$mRNA) %in% gene.symbol)
        temp=data$mRNA[index,]
        
        #####calculate similarity matrix
        ddist <- as.matrix(dist2(as.matrix(t(temp)), as.matrix((t(temp)))))
        distanceMatrix = affinityMatrix(ddist, 20, 0.5)
        attr(distanceMatrix,'class')="Similarity"
        
        IC10_res = ExecuteIntClust(data$mRNA)
        fpath = getFilePath(fileFolder, paste(i,"IC10",sep = "_"))
        save(IC10_res, file = fpath)
        
        timeTable[i,"IntClust"] = IC10_res$timing
        silTable[i,"IntClust"] = getMeanSilhouette(IC10_res$group, distanceMatrix)
        pvalueTable[i,"IntClust"] = getPvalue(dclin[,1],dclin[,2],IC10_res$group)
        
        cat("finish IntClust.\n")
      }
      else {
        cat("Need to provide gene expression data.\n")
      }
    }
    
    ###CC
    if ("CC" %in% methods) {
      cat("run CC ...\n")
      CC_res=ExecuteCC(clusterNum=5,d=data1,maxK=5,clusterAlg="hc",
                       distance="pearson",title="CC_res")
      fpath = getFilePath(fileFolder, paste(i,"CC",sep = "_"))
      save(CC_res, file = fpath)
      
      timeTable[i,"CC"] = CC_res$timing
      silTable[i,"CC"] = getMeanSilhouette(CC_res$group, CC_res$distanceMatrix)
      pvalueTable[i,"CC"] = getPvalue(dclin[,1],dclin[,2],CC_res$group)
      
      cat("finish CC.\n")
    }
    
    ###CNMF
    if ("CNMF" %in% methods) {
      cat("run CNMF ...\n")
      CNMF_res=ExecuteCNMF(data2,clusterNum=5,nrun=30)
      fpath = getFilePath(fileFolder, paste(i,"CNMF",sep = "_"))
      save(CNMF_res, file = fpath)
      
      timeTable[i,"CNMF"] = CNMF_res$timing
      silTable[i,"CNMF"] = getMeanSilhouette(CNMF_res$group, CNMF_res$distanceMatrix)
      pvalueTable[i,"CNMF"] = getPvalue(dclin[,1],dclin[,2],CNMF_res$group)
      
      cat("finish CNMF.\n")
    }
    
    ###iCluster
    if ("iCluster" %in% methods) {
      cat("run iCluster ...\n")
      
      temp=NULL
      if(length(data1) > 1)
      {
        for(j in 1: length(data1))
        {
          temp=rbind(temp,data1[[j]])
        }
      } else temp=data1[[1]]
      
      
      ddist <- as.matrix(dist2(as.matrix(t(temp)), as.matrix((t(temp)))))
      distanceMatrix = affinityMatrix(ddist, 20, 0.5)
      attr(distanceMatrix,'class')="Similarity"
      
      iCluster_res=ExecuteiCluster(datasets=data1, k=5, lambda=list(0.44,0.33,0.28))
      fpath = getFilePath(fileFolder, paste(i,"iCluster",sep = "_"))
      save(iCluster_res, file = fpath)
      
      timeTable[i,"iCluster"] = iCluster_res$timing
      silTable[i,"iCluster"] = getMeanSilhouette(iCluster_res$group, distanceMatrix)
      pvalueTable[i,"iCluster"] = getPvalue(dclin[,1],dclin[,2],iCluster_res$group)
      
      cat("finish iCluster.\n")
    }
    
    ###SNF
    if ("SNF" %in% methods) {
      cat("run SNF ...\n")
      SNF_res=ExecuteSNF(datasets=data,clusterNum=5)
      fpath = getFilePath(fileFolder, paste(i,"SNF",sep = "_"))
      save(SNF_res, file = fpath)
      
      timeTable[i,"SNF"] = SNF_res$timing
      silTable[i,"SNF"] = getMeanSilhouette(SNF_res$group, SNF_res$distanceMatrix)
      pvalueTable[i,"SNF"] = getPvalue(dclin[,1],dclin[,2],SNF_res$group)
      
      cat("finish SNF.\n")
    }
    
    ###SNF-CC
    if ("SNF-CC" %in% methods) {
      cat("run SNF-CC ...\n")
      SNF_CC_res=ExecuteSNF.CC(data1, clusterNum=5, K=20, alpha=0.5, t=20,
                               maxK = 5, pItem = 0.8,reps=500, 
                               title = "SNF_CC_res", plot = "png", 
                               finalLinkage ="average")
      fpath = getFilePath(fileFolder, paste(i,"SNF_CC",sep = "_"))
      save(SNF_CC_res, file = fpath)
      
      timeTable[i,"SNF-CC"] = SNF_CC_res$timing
      silTable[i,"SNF-CC"] = getMeanSilhouette(SNF_CC_res$group, SNF_CC_res$distanceMatrix)
      pvalueTable[i,"SNF-CC"] = getPvalue(dclin[,1],dclin[,2],SNF_CC_res$group)
      
      cat("finish SNF-CC.\n")
    }
    
    if ("WSNF" %in% methods) {
      cat("run WSNF ...\n")
      data(Ranking)
      if (omics == "mRNA"){
        ####Retrieve the feature ranking for genes
        gene_Name=rownames(data$mRNA)
        index1=match(gene_Name,Ranking$mRNA_TF_miRNA.v21_SYMBOL)
        gene_ranking=data.frame(gene_Name,Ranking[index1,],stringsAsFactors=FALSE)
        index2=which(is.na(gene_ranking$ranking_default))
        gene_ranking$ranking_default[index2]=min(gene_ranking$ranking_default,na.rm =TRUE)
        
        ###Clustering
        ranking1=list(gene_ranking$ranking_default)
      }
      
      else if (omics == "miRNA"){
        
        ####Retrieve the feature ranking for miRNAs
        miRNANames=rownames(data$miRNA)
        version=checkMiRNAVersion(miRNANames)
        Accessions=miRNA_NameToAccession(miRNANames, version=version)
        index3=match(Accessions$Accession,Ranking$mRNA_TF_miRNA_ID)
        
        miRNA_ranking=data.frame(Accessions$Accession,Ranking[index3,],stringsAsFactors=FALSE)
        index4=which(is.na(miRNA_ranking$ranking_default))
        miRNA_ranking$ranking_default[index4]=min(miRNA_ranking$ranking_default,na.rm =TRUE)
        ###Clustering
        ranking1=list(miRNA_ranking$ranking_default)
      }
      
      else {
        ####Retrieve the feature ranking for genes
        gene_Name=rownames(data$mRNA)
        index1=match(gene_Name,Ranking$mRNA_TF_miRNA.v21_SYMBOL)
        gene_ranking=data.frame(gene_Name,Ranking[index1,],stringsAsFactors=FALSE)
        index2=which(is.na(gene_ranking$ranking_default))
        gene_ranking$ranking_default[index2]=min(gene_ranking$ranking_default,na.rm =TRUE)
        
        ####Retrieve the feature ranking for miRNAs
        miRNANames=rownames(data$miRNA)
        version=checkMiRNAVersion(miRNANames)
        Accessions=miRNA_NameToAccession(miRNANames, version=version)
        index3=match(Accessions$Accession,Ranking$mRNA_TF_miRNA_ID)
        
        miRNA_ranking=data.frame(Accessions$Accession,Ranking[index3,],stringsAsFactors=FALSE)
        index4=which(is.na(miRNA_ranking$ranking_default))
        miRNA_ranking$ranking_default[index4]=min(miRNA_ranking$ranking_default,na.rm =TRUE)
        ###Clustering
        ranking1=list(gene_ranking$ranking_default, miRNA_ranking$ranking_default)
      }
      
      WSNF_res=ExecuteWSNF(datasets=data, feature_ranking=ranking1,
                           beta = 0.8, clusterNum=5, 
                           K = 20,alpha = 0.5, t = 20, plot = TRUE)
      fpath = getFilePath(fileFolder, paste(i,"WSNF",sep = "_"))
      save(WSNF_res, file = fpath)
      
      timeTable[i,"WSNF"] = WSNF_res$timing
      silTable[i,"WSNF"] = getMeanSilhouette(WSNF_res$group, WSNF_res$distanceMatrix)
      pvalueTable[i,"WSNF"] = getPvalue(dclin[,1],dclin[,2],WSNF_res$group)
      
      cat("finish WSNF.\n")
    }
    
    if ("CIMLR" %in% methods) {
      cat("run CIMLR ...\n")
      CIMLR_res=ExecuteCIMLR(datasets=data,clusterNum=5,ncore = ncore-2)
      fpath = getFilePath(fileFolder, paste(i,"CIMLR",sep = "_"))
      save(CIMLR_res, file = fpath)
      
      timeTable[i,"CIMLR"] = CIMLR_res$timing
      silTable[i,"CIMLR"] = getMeanSilhouette(CIMLR_res$group, CIMLR_res$distanceMatrix)
      pvalueTable[i,"CIMLR"] = getPvalue(dclin[,1],dclin[,2],CIMLR_res$group)
      
      cat("finish CIMLR.\n")
    }
    
    if ("PINS" %in% methods) {
      cat("run PINS ...\n")
      PINS_res=ExecutePINS(datasets=data,clusterNum=5,ncore = ncore-1)
      fpath = getFilePath(fileFolder, paste(i,"PINS",sep = "_"))
      save(PINS_res, file = fpath)
      
      ###calculate similarity matrix
      temp=NULL
      if(length(data) > 1)
      {
        for(j in 1: length(data))
        {
          temp=rbind(temp,data[[j]])
        }
      }
      else
        temp=data[[1]]
      
      ddist <- as.matrix(dist2(as.matrix(t(temp)), as.matrix((t(temp)))))
      distanceMatrix = affinityMatrix(ddist, 20, 0.5)
      attr(distanceMatrix,'class')="Similarity"
      
      timeTable[i,"PINS"] = PINS_res$timing
      silTable[i,"PINS"] = getMeanSilhouette(PINS_res$group, distanceMatrix)
      pvalueTable[i,"PINS"] = getPvalue(dclin[,1],dclin[,2],PINS_res$group)
      
      cat("finish PINS.\n")
    }
    
    ####NEMO
    if ("NEMO" %in% methods) {
      cat("run NEMO ...\n")
      NEMO_res=ExecuteNEMO(datasets=data,clusterNum=5)
      fpath = getFilePath(fileFolder, paste(i,"NEMO",sep = "_"))
      save(NEMO_res, file = fpath)
      
      timeTable[i,"NEMO"] = NEMO_res$timing
      silTable[i,"NEMO"] = getMeanSilhouette(NEMO_res$group, NEMO_res$distanceMatrix)
      pvalueTable[i,"NEMO"] = getPvalue(dclin[,1],dclin[,2],NEMO_res$group)
      
      cat("finish NEMO.\n")
    }
    
    fpath = getFilePath(fileFolder, paste0("timeTable",".csv"))
    write.csv(timeTable, file = fpath)
    
    fpath = getFilePath(fileFolder, paste0("silTable",".csv"))
    write.csv(silTable, file = fpath)
    
    fpath = getFilePath(fileFolder, paste0("pvalueTable",".csv"))
    write.csv(pvalueTable, file = fpath)
  }
  
  sink()
  return(list("timeTable" = timeTable, "silTable" = silTable, "pvalueTable" = pvalueTable))
}
