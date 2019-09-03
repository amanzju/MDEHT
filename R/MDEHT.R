##-------------------------------------------------------------------
## 	Name: MDETH.R
## 	"A multivariate approach for detecting differential expression of microRNA isoforms in RNA sequencing studies"
## 	Date: ...........
## 	Contact: .........

## Input:
#X:  a matrix of one group (treat as tumor sample group) of RNA-Seq isomiRs counts (isomiRs are rows, samples are columns)
#Y:  a matrix of another group (treat as normal sample group) of RNA-Seq isomiRs counts (isomiRs are rows, samples are columns)
#    when Y =NULL then it will perform one sample Hotelling's test
#we have provided log2(cpm(...)+1) normalized data for X and Y
#inohf.isomiRs: 	a logical parameter. if you want include only high frequency isoforms in a miRNA when it contains large number of isomiRs. Default FALSE
#inohf.isomiRs.no: only used when inohf.isomiRs is TRUE and default value 25
#merge.olf.isomRs: a logical parameter. if you want merge only low frequency isoforms as a new isoform from a miRNA when it contains large number of isomiRs.
#   Default FALSE
#merge.olf.isomRs.no: only used when merge.olf.isomRs is TRUE and default value 25
#adjp.method:	p values adjusted method. Defaults to "BH".


## Output:
#tumor.mean: a vector of mean values from tumor group
#normal.mean: a vector of mean values from normal group
#T.value: a vector of Hotellig's T-statistic values
#p.value: a vector containing the raw p-values testing differential expression for each miRNA
#adj.pval: a vector containing the p-values after adjusting for multiple testing using the method of Benjamini-Hochberg
#FC:		a vector containing the estimated log fold changes for each miRNA



################# Function start    ######################

MDEHT<-function(X, Y, inohf.isomiRs = FALSE, inohf.isomiRs.no = 25,
  merge.olf.isomRs = FALSE, merge.olf.isomRs.no = 25,adjp.method="BH"){
  if(is.null(Y)){print("Peforming one sample Hotelling's T-square test")
          } else{print("Peforming two sample Hotelling's T-square test")}

 ################  Two sample test  #########################
 if(is.matrix(Y) | is.data.frame(Y)){

  #intersect isoforms
  XintYiso<-intersect(rownames(X),rownames(Y))
  X <- X[XintYiso,]
  Y <- Y[XintYiso,]

  # extract miRNAs
  isoNamesplit<-strsplit(rownames(X),"\\|")
  miRs<-unique(sapply(isoNamesplit,function(x)x[1]))

  lenmiR<-NULL
  for(i in 1:length(miRs)){
    echmiR<-miRs[i]
    lenmiR[i]<-length(grep(echmiR,rownames(X)))}

  mxmiR<-max(lenmiR)
  TSS<-(ncol(X)+ncol(Y))
  assum<-mxmiR<(TSS-2)
  if(!isTRUE(inohf.isomiRs)&!isTRUE(merge.olf.isomRs)&!isTRUE(assum))stop("Your maxmum number of isoform in a miRNA is ",mxmiR,
      "\n  Total number of sample is ",TSS,
      "\n  So, you must be provide maximum number of isoform in miRNA less than ",(TSS-2),
      "\n  Or you can set 'TRUE' for parameter 'inohf.isomiRs' or 'merge.olf.isomRs' to overcome this problem.")

   #############  Fold change value  ################
   Xbind<-NULL
   Ybind<-NULL

   for(j in 1:length(miRs)){
    mir<-miRs[j]
    No.iso<-grep(mir,rownames(X))
    if(length(No.iso)>1){
      Xbind<-rbind(Xbind,colSums(as.matrix(X[grep(mir,rownames(X)),]))) # sum up all isoforms for a miRNA in tumor group
      Ybind<-rbind(Ybind,colSums(as.matrix(Y[grep(mir,rownames(Y)),]))) # sum up all isoforms for a miRNA in normal group
     }else {
        Xbind<-rbind(Xbind,colSums(t(as.matrix(X[grep(mir,rownames(X)),])))) # sum up all isoforms for a miRNA in tumor group
        Ybind<-rbind(Ybind,colSums(t(as.matrix(Y[grep(mir,rownames(Y)),])))) # sum up all isoforms for a miRNA in normal group
       }
    }
   rownames(Xbind)<-rownames(Ybind)<-miRs

   TM<-apply(Xbind,1,function(x)mean(x))         # calculate the mean values of all miRNA in tumor group
   NM<-apply(Ybind,1,function(x)mean(x))         # calculate the mean values of all miRNA in normal group
   fc<-apply(Xbind,1,mean)/apply(Ybind,1,mean)   # fold change values

   ts<-ncol(Xbind)        # Tumor sample size
   ns<-ncol(Ybind)        # Normal sample size
   DataMatrix<-cbind(X,Y) # bind both group datasets
   group<-as.factor(c(rep("tumor",ts),rep("normal",ns))) # create group labels

 ########## Hoteling Test  ###########
    pvalue<-NULL
    Tvalue<-NULL

    for(k in 1:length(miRs)){
       mi<-miRs[k]
       dm<-as.matrix(DataMatrix[grep(mi,rownames(DataMatrix)),])
       dm<-if(dim(dm)[2]==1)(dm=t(dm))else(dm=dm)
     # exclude an isoform if it contain 90 percent 0 across the all sample
       if(nrow(dm)>1){
          keep<-(rowSums(dm==0)/ncol(dm))<=0.9
          #keep<-apply(dm ,1,function(x){(sum(x==0)/ncol(dm))<=0.9})
          if(sum(keep)==0){dm <-dm
           }else
            {dm <-dm[keep,]
             if(is.vector(dm)){dm<-t(as.matrix(dm))}
             }
        }else{dm<-dm}

      #remove low friquency isoforms
       if(inohf.isomiRs){
        if(dim(dm)[1]>inohf.isomiRs.no){
         isoM<-as.vector(rowMeans(dm))
         dm<-as.data.frame(cbind(dm,isoM))
	   dm<-dm[order(dm$isoM,decreasing=T),]
         dm<-dm[1:inohf.isomiRs.no,]
         dm$isoM<-NULL
        }else{dm<-dm}
       }


      #merge low frequency isoforms as a single isoform
       if(merge.olf.isomRs){
        if(dim(dm)[1]>merge.olf.isomRs.no){
	   isoM<-as.vector(rowMeans(dm))
         dm<-as.data.frame(cbind(dm,isoM))
	   dm<-dm[order(dm$isoM,decreasing=T),]
         hed<-dm[1:(merge.olf.isomRs.no-1),]
         led<-colSums(dm[-c(1:(merge.olf.isomRs.no-1)),])
         dm<-rbind(hed,led)
         dm$isoM<-NULL
        }else{dm<-dm}
       }


  if(dim(dm)[1]>(dim(dm)[2]-2))stop("Number of observation p is greater than n-2.",
       "\n  Reduce the number in parameter 'inohf.isomiRs' or 'merge.olf.isomRs' if any of them is TRUE.")

     # Hotelling's test start
       dataframe<-as.data.frame(cbind(group,t(dm)))
       split.data = split(dataframe[,-1],dataframe[,1])
       x1 <- as.matrix(split.data[[1]])
       x2 <- as.matrix(split.data[[2]])
       p <- ncol(x1)
       n1 <- nrow(x1)
       n2 <- nrow(x2)
       n <- n1 + n2
       const<-((n1*n2)/n)
       xbar1 <- colMeans(x1)
       xbar2 <- colMeans(x2)
       dbar <- as.vector(xbar2-xbar1)
       vcm <- ((n1-1)*cov(x1)+(n2-1)*cov(x2))/(n-2)
       dtm<-det(vcm)
       tvalue<-tryCatch(
         if (dtm==0){
          t2 <- const*(dbar%*%ginv(vcm)%*%dbar)
          }else
           {t2 <-  const*(dbar %*% solve(vcm)%*%dbar)}
           ,error=function(e){return(NA)}
         )

	     Tvalue[k]<-tvalue
        if(is.na(tvalue)){pvalue[k]<-NA
         }else
         { test<-as.vector((n-p-1)*tvalue/((n-2)*p))
           pvalue[k]<-pf(test,p,n-p-1,lower.tail=FALSE)
         }
     }

 pvalAdj<-p.adjust(pvalue,method=adjp.method)
  output<-cbind(tumor.mean=TM,normal.mean=NM,T.value=Tvalue,p.value=pvalue,adj.pval=pvalAdj,FC=fc)
  rownames(output)<-gsub("\\_","",miRs)
  }  ###### end loop for two sample test



###########One sample test##################

if(is.null(Y)){
  # extract miRNAs
  isoNamesplit<-strsplit(rownames(X),"\\|")
  miRs<-unique(sapply(isoNamesplit,function(x)x[1]))

  lenmiR<-NULL
  for(i1 in 1:length(miRs)){
    echmiR<-miRs[i1]
    lenmiR[i1]<-length(grep(echmiR,rownames(X)))}

  mxmiR<-max(lenmiR)
  TSS<-ncol(X)
  assum<-mxmiR<TSS

  if(!isTRUE(merge.olf.isomRs)&!isTRUE(inohf.isomiRs)&!isTRUE(assum))stop("Your maxmum number of isoform in a miRNA is ",mxmiR,
      "\n  Total number of sample is ",TSS,
      "\n  So, you must be provide maximum number of isoform in miRNA less than ",(TSS-2),
      "\n  Or you can set 'TRUE' for parameter 'inohf.isomiRs' or 'merge.olf.isomRs' to overcome this problem")

    ########## Hoteling Test  ###########

    pvalue<-NULL
    Tvalue<-NULL

    for(j1 in 1:length(miRs))
      {
       mi<-miRs[j1]
       dm<-as.matrix(X[grep(mi,rownames(X)),])
       dm<-if(dim(dm)[2]==1)(dm=t(dm))else(dm=dm)
     # exclude an isoform if it contain 90 percent 0 across the all sample
      if(nrow(dm)>1){
           keep<-(rowSums(dm==0)/ncol(dm))<=0.9
           #keep<-apply(dm ,1,function(x){(sum(x==0)/ncol(dm))<=0.9})
           if(sum(keep)==0){dm <-dm
           }else
            {dm <-dm[keep,]
             if(is.vector(dm)){dm<-t(as.matrix(dm))}
             }
      }else{dm<-dm}



      #remove low friquency isoforms
       if(inohf.isomiRs){
        if(dim(dm)[1]>inohf.isomiRs.no){
         isoM<-as.vector(rowMeans(dm))
         dm<-as.data.frame(cbind(dm,isoM))
	   dm<-dm[order(dm$isoM,decreasing=T),]
         dm<-dm[1:inohf.isomiRs.no,]
         dm$isoM<-NULL
        }else{dm<-dm}
       }


      #merge low frequency isoforms as a single isoform
       if(merge.olf.isomRs){
        if(dim(dm)[1]>merge.olf.isomRs.no){
	   isoM<-as.vector(rowMeans(dm))
         dm<-as.data.frame(cbind(dm,isoM))
	   dm<-dm[order(dm$isoM,decreasing=T),]
         hed<-dm[1:(merge.olf.isomRs.no-1),]
         led<-colSums(dm[-c(1:(merge.olf.isomRs.no-1)),])
         dm<-rbind(hed,led)
         dm$isoM<-NULL
        }else{dm<-dm}
       }


  if(dim(dm)[1]>(dim(dm)[2]-1))stop("Number of observation p is greater than n-1.",
       "\n  Reduce the number in parameter 'inohf.isomiRs' or 'merge.olf.isomRs' if any of them is TRUE.")


    #################################
       x <- as.matrix(t(dm))
       m <- colMeans(x)
       s<-cov(x)
       n <- nrow(x)
       p <- ncol(x)
       const<- (n*(n-p))/((n-1)*p)
       dtm<-det(s)
       test<-tryCatch(
         if (dtm==0){
          t2 <- as.vector(const*(m%*%ginv(s)%*%m))
          }else
           {t2 <-  as.vector(const*(m%*%solve(s,m)))}
           ,error=function(e){return(NA)}
         )

	   Tvalue[j1]<-test
        if(is.na(test)){pvalue[j1]<-NA
         }else
         {pvalue[j1]<-pf(test,p,n-p,lower.tail=FALSE)}
     }

  pvalAdj<-p.adjust(pvalue,method="BH")
  output<-cbind(T.value=Tvalue,p.value=pvalue,adj.pval=pvalAdj)
  rownames(output)<-gsub("\\_","",miRs)
  }  ###### end loop for one sample test

  output<-as.data.frame(output)
  return(output)

 }

 ########### End main function loop
