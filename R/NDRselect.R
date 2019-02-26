#' @useDynLib NDRindex
#' @importFrom Rcpp sourceCpp
##NULL

TMM<-function(Data)
{
  sec<-calcNormFactors(as.matrix(Data),method="TMM")
  res<-log2(t(t(Data)/sec)+1)
  tsec<-sec
  for (i in 1:length(sec))
  {
    if (sec[i]>0)
      tsec[i]=1
    else
      tsec[i]=-1
  }
  sec=abs(sec)
  res<-log2(t(t(Data)/sec)+1)
  res<-t(t(Data)/tsec)
  return (res)
}

scarn<-function(Data)
{
  sec <- computeSumFactors(as.matrix(Data))
  tsec<-sec
  for (i in 1:length(sec))
  {
    if (sec[i]>0)
      tsec[i]=1
    else
      tsec[i]=-1
  }
  sec=abs(sec)
  res<-log2(t(t(Data)/sec)+1)
  res<-t(t(Data)/tsec)
  return(res)
}

pca<-function(Data)
{
  res <- prcomp(t(Data))
  res <- predict(res)
  res2 <- res[,1:2]
  return(res2)
}

seurat<-function(Data)
{
  pbmc <- CreateSeuratObject(raw.data = Data, min.cells = 0, min.genes = 0)
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize")
  return(pbmc@raw.data)
}

#' Test all combiantion of normalization and dimensional reduction tools and find the best one
#'
#' @param Norm A character array contains commands to call normalization tools. The arguments that represents the input data should be replaced by 'Data'.
#' @param Red A character array contains commands to call dimensional reduction tools. The arguments that represents the input data should be replaced by 'Data'.
#' @param NormReturn A character array contains the command that can extract data from normalization results and change it to right form that following dimension reduction tools could analysis.Use 'Data' to represent normalization results.
#' @param RedReturn A character array contains the command that can extract data from dimensional reduction tools results and change it to right form that following NDRindex could analysis.Use 'Data' to represent dimensional reduction results.
#' @param NormName A character array contains the name of normalization tools.
#' @param RedName A character array contains the name of dimensional reduction tools.
#' @param ClusterName A string represent the clustering methods you want to use. 'kmeans' and 'hclust' are available
#' @param ClusterNumber A number represent the number of clusters you want.
#' @return
#' A list with three elements.
#' $reslist contains all combinations' name and the NDRindex of them.
#' $resData contains the best producing data of all combinations.
#' $cluster contians the clustering result.
#' @examples
#' NDRselect(Data,
#'          Norm=c('betweenLaneNormalization(as.matrix(Data), which=\"full\")','scale(Data)','Linnorm(as.matrix(Data))','TMM(Data)','scarn(Data)','seurat(Data)'),
#'          Red=c('Rtsne(t(Data),dim=2,perplexity = 15,check_duplicates=FALSE)','pca(Data)','Data<-sammon(d=dist(t(Data),method=\"euclidean\"),k=2)'),
#'          NormReturn=c('as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)'),
#'          RedReturn=c('Data$Y','as.matrix(Data)','Data$points'),
#'          NormName=c('EDASeq','scale','Linnorm','TMM','scarn','seurat'),
#'          RedName=c('tSne','pca','sammon'),
#'          ClusterName='kmeans',
#'          ClusterNumber=4)
#' @export
NDRselect<- function(Data,
                    Norm=c('betweenLaneNormalization(as.matrix(Data), which=\"full\")','scale(Data)','Linnorm(as.matrix(Data))','TMM(Data)','scarn(Data)','seurat(Data)'),
                    Red=c('Rtsne(t(Data),dim=2,perplexity = 15,check_duplicates=FALSE)','pca(Data)','Data<-sammon(d=dist(t(Data),method=\"euclidean\"),k=2)'),
                    NormReturn=c('as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)','as.matrix(Data)'),
                    RedReturn=c('Data$Y','as.matrix(Data)','Data$points'),
                    NormName=c('EDASeq','scale','Linnorm','TMM','scarn','seurat'),
                    RedName=c('tSne','pca','sammon'),
                    ClusterName='kmeans',
                    ClusterNumber=4)
{
  index=c(NA)
  NormN=c(NA)
  RedN=c(NA)
  NDRselectreslist<-data.frame(index,NormN,RedN)
  nowrow=1
  maxindex=0
  tData<-Data
  rewData<-Data
  NDRselectres<-list()
  for (i in 1:length(Norm))
    for (j in 1:length(Red))
    {
      print(Norm[i])
      print(NormReturn[i])
      print(Red[j])
      print(RedReturn[j])
      Data<-eval(parse(text=Norm[i]))
      Data<-eval(parse(text=NormReturn[i]))
      Data<-eval(parse(text=Red[j]))
      Data<-eval(parse(text=RedReturn[j]))
      NDRselectreslist[nowrow,1]=NDRindex(Data,nrow(Data),ncol(Data))
      NDRselectreslist[nowrow,2]=NormName[i]
      NDRselectreslist[nowrow,3]=RedName[j]
      if (NDRselectreslist[nowrow,1]>maxindex)
      {
        resData<-Data
        maxindex-NDRselectreslist[nowrow,1]
      }
      nowrow=nowrow+1
      Data<-tData
    }
  NDRselectreslist<-dplyr::arrange(NDRselectreslist,dplyr::desc(index))
  NDRselectres$reslist<-NDRselectreslist
  NDRselectres$resData<-resData
  if (ClusterName=='kmeans')
    NDRselectres$cluster<-kmeans(resData,ClusterNumber)$cluster
  if (ClusterName=='hclust')
  {
    idis<-dist(resData,method = "euclidean")
    clust2_1norm<-hclust(idis,method="average")
    NDRselectres$cluster<-cutree(clust2_1norm,k=ClusterNumber)
  }
  return(NDRselectres)
}
