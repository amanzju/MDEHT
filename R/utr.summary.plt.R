#It will produce a distribution of 3 prime and 5 prime isomiRs compare to architype miRNA

utr.summary<-function(X){
  # extract miRNA names
  isomiRs<-matrix(unlist(strsplit(rownames(X),"\\|")),ncol=3,byrow=TRUE)
  miRNAs<-names(table(isomiRs[,1]))
 # plot for the distribution of 3 UTR and 5 UTR at both side of architype miRNA
  p3UTR<-isomiRs[grep("3p",isomiRs[,1]),]
  p5UTR<-isomiRs[grep("5p",isomiRs[,1]),]

  mp1<-intersect(names(table(p3UTR[,2])),names(table(p5UTR[,2])))
  leftmatrix<-rbind(table(p3UTR[,2])[mp1],table(p5UTR[,2])[mp1])
  leftmatrix<- leftmatrix[, order(as.integer(colnames(leftmatrix)))]

  mp2<-intersect(names(table(p3UTR[,3])),names(table(p5UTR[,3])))
  rightmatrix<-rbind(table(p3UTR[,3])[mp2],table(p5UTR[,3])[mp2])
  rightmatrix<- rightmatrix[, order(as.integer(colnames(rightmatrix)))]

  # extract miRNA names from full dataset
  IsoNo<-as.matrix(table(isomiRs[,1]))
  rownames(IsoNo)<-gsub("\\_","",miRNAs)
  colnames(IsoNo)<-"n.isoform"

  par(mfrow=c(2,1))
  barplot(leftmatrix,beside=T,col=c("greenyellow","darkorchid1"),ylab="No. of isomiRs",
	xlab="Genomic position relative to archetype miRNA",main="left side of isomiRs",las=2)
  legend("topright",c("3' endpoints","5' endpoints"),fill = c("greenyellow","darkorchid1"),bty = "n")

  barplot(rightmatrix,beside=T,col=c("greenyellow","darkorchid1"),ylab="No. of isomiRs",
	xlab="Genomic position relative to archetype miRNA",main="right side of isomiRs",las=2)
  legend("topright",c("3' endpoints","5' endpoints"),fill = c("greenyellow","darkorchid1"),bty = "n")

 }