#It produce a distribution isoforms number in different miRNAs from raw datasets.
isomiR.summary<-function(X){
  # extract miRNA names  from selected dataset  
  isomiR<-matrix(unlist(strsplit(rownames(X),"\\|")),ncol=3,byrow=TRUE)
  miRNA<-names(table(isomiR[,1]))
  isoNo<-as.matrix(table(isomiR[,1]))
  rownames(isoNo)<-gsub("\\_","",miRNA)
  colnames(isoNo)<-"n.isoform"
  barplot(table(isoNo[,1]),xlab="Number of isoform in each miRs",las=2,main="isomiRs Histogram",ylab="Number of miRs")
 }