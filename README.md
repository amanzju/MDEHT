# Description
The MDEHT package executes multivariate Hotelling’s T-square test in isomiRs RNA-Seq data to detect differential expression miRNA. Where Hotelling’s T-square test is a generalization of Student's t- statistic and is widely used for testing the difference between two multivariate means. Currently, the MDEHT package performs both one or two samples multivariate test. In one sample test, the null hypothesis defined as the vector of the mean value of multiple variables is equal to zero. Under this package, we have provided a tow sample random dataset to implements MDEHT function.


## Access
MDEHT and isomiRseeker are free software tools that can be used by academic groups for non-commercial purpose. 

## System requirements
> Linux/Unix or Mac OS platform  
> R (https://www.r-project.org/)  
isomiRseeker required the following software to be installed in the user’s environment.  
> samtools (http://www.htslib.org/download/)  
> bedtools (https://bedtools.readthedocs.io/en/stable/content/installation.html)  
> bwa (https://github.com/lh3/bwa)  
> Perl 5 (http://www.perl.org/) 

## Usage isomiRseeker 
To generate an expression of isomiR from reading alignments in bam file format, run the following code:

perl /path/isomiRseeker.pl -b -p output_file_name input_file.bam

**Note:** Output from the above code will produce two files for each sample one is count format and another one is rpm format on your working directory. And finally, process the count/cpm data as in a matrix format for all sample for a specific disease type.

### Example

perl /path/isomiRseeker.pl -b -p TCGA-2F-A9KQ-01A-11R-A38M-13_mirna TCGA-2F-A9KQ-01A-11R-A38M-13_mirna.bam

## Usage MDEHT
The normalized of isomiRs read count expression data should be a matrix with isomiRs in rows and samples in columns. The row names should be same as our isomiRseeker tool generated isomiRs symbols. We suggest, most widely used CMP normalization method to apply in isomiRs data before run MDEHT method.  


### Example
#required the following packages
> install.packages("devtools")  
> library(devtools)  
> install_github("amanzju/MDEHT")  
> library(edgeR)  
> library (MDEHT) 

#load an example dataset
> data("ExampleData") 

#use CPM normalization
> NormalizedData = log2(cpm(ExampleData)+1) 

#separate sample group  
> TumorData = NormalizedData[,grep("01A", colnames(NormalizedData))]  
> NormalData = NormalizedData[,grep("11A", colnames(NormalizedData))] 

#view the summary of isomiRs and UTR in plot  
> isomiRs<-rownames(TumorData)    
> isomiR.summary(x = isomiRs)   
> utr.summary(x = isomiRs)   

> results = MDEHT (X = TumorData, Y = NormalData)  
> head(results)   

## Output
The output of MDEHT is a summary of putative miRNAs, including the mean values for two groups, the value of the T-statistic, p-value and FDR adjusted p-value and log2 fold change value. Similarly, for one sample test, the output will provide the value of the T-statistic, p-value, and FDR adjusted p-value.


## Restriction of MDEHT
MDEHT will stop when the number of isoforms in miRNA greater than n1 + n2 – 2, where n1 and n2 are sample sizes of two different group. This restriction can be released by using
the technical parameter inohf.isomiRs and merge.olf.isomRs in MDEHT function. Where inohf.isomiRs: represent include only high-frequency isoforms in a miRNA when it contains a large number of isomiRs and merge.olf.isomRs: represent merge only low-frequency isoforms as a new isoform from a miRNA when it contains a large number of isomiRs. 
