#' get file path based on a given folder name and a given file name
#'
#' @param foldername A given folder name; if NULL, will set to the current working directory.
#' @param filename A given file name
#' @return The full path of the file
#' 
getFilePath<-function(foldername=NULL,filename=NULL)
{
  mainDir<-getwd()
  filename1=filename
  if(!is.null(foldername))
  {
    subDir=foldername
    dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
    fpath=file.path(mainDir, subDir,filename1)
  }
  else
    fpath=file.path(mainDir,filename1)  
  fpath
}

#' load R data from a given directory
#'
#' @param filefolder A number
#' @param dataSets names of datasets; a string vector.
#' @examples
#' \dontrun{
#' dn = c("TCGA","UK","HEL","GSE19783")
#' loadData("data",dn)
#' }
#' 
loadData <- function(filefolder=NULL, dataSets = NULL){
  if(!is.na(dataSets)) {
    for(i in dataSets) {
      fpath = getFilePath(filefolder, i)
      load(fpath,.GlobalEnv)
    }
  }
  else {
    cat("Need to provide name(s) of data set(s).\n")
  }
}
