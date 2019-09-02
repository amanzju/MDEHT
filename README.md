# Description
Package implementing MDEHT method to perform a multivariate differential expression by Hotelling’s test using isomiRs in small RNA-Seq data. Hotelling’s test is a generalization of Student's t- statistic and is widely used for testing the difference between two multivariate means.

## Access
MDEHT and isomiRseeker are free software tools that can be used by academic groups for non-commercial purpose. 

## System requirements
> Linux/Unix or Mac OS platform

> R (https://www.r-project.org/) 

> isomiRseeker required the following software to be installed in the user’s environment.  

> samtools (http://www.htslib.org/download/)

> bedtools (https://bedtools.readthedocs.io/en/stable/content/installation.html)

> bwa (https://github.com/lh3/bwa)

> Perl 5 (http://www.perl.org/)


## Usage isomiRseeker 
To generate expression of isomiR from read alignments in bam file format, run the following code:

> perl /path/isomiRseeker.pl -b -p output_file_name input_file.bam

### Note:
output from the above code will produce two files for each sample one is count format and another one is rpm format on your working directory. And finally, process the count/cpm data as in a matrix format for all sample for a specific disease type.

### Example 
perl /path/isomiRseeker.pl -b -p TCGA-2F-A9KQ-01A-11R-A38M-13_mirna TCGA-2F-A9KQ-01A-11R-A38M-13_mirna.bam

## Usage MDEHT
The input parameters of the main function of MDEHT as two matrices (two-sample test) or a matrix (one-sample test) of the CPM normalized of isomiRs read count data. Each data matrix row represents the isomiRs and column represent the sample. Firstly, load the MDEHT function as source code then run MDEHT function.

### Example 
#required the following packages

> install.packages("devtools")

> library(devtools)

> install_github("amanzju/MDEHT")

> library(edgeR)

> library(MDEHT)


#load an example dataset

> data("ExampleData")

#use CPM normalization

> NormalizedData = log2(cpm(ExampleData)+1)

#separate sample group

> TumorData = NormalizedData[,grep("01A", colnames(NormalizedData))]

> NormalData = NormalizedData[,grep("11A", colnames(NormalizedData))]

#view the summary of isomiRs and UTR in plot

> isomiR.summary(X = TumorData)

> utr.summary(X = TumorData)

> results = MDEHT (X = TumorData, Y = NormalData)

> head(results)

By default, all parameter under MDETH function is FALSE. the output is a data frame that contains the various measurement of all miRNAs.

## Output
The output of MDEHT is a summary of putative miRNAs, including the mean values for two groups, the value of the T-statistic, p-value and FDR adjusted p-value and log2 fold change value. Similarly, for one sample test, the output will provide the value of the T-statistic, p-value, and FDR adjusted p-value.

## Restriction of MDEHT
MDEHT will stop when the number of isoforms in miRNA greater than n1 + n2 – 2, where n1 and n2 are sample sizes of two different group. This restriction can be released by using
the parameter inohf.isomiRs and merge.olf.isomRs in MDEHT function. Where inohf.isomiRs: represent include only high-frequency isoforms in a miRNA when it contains a large number of isomiRs and merge.olf.isomRs: represent merge only low-frequency isoforms as a new isoform from a miRNA when it contains a large number of isomiRs.
