# ieITH

A set of algorithms for measuring ITH in tumors using one of the multi-omics profiles. These algorithms are designed based on  information entropy. 

## Description

The package is used to calculate IntraTumor Heterogeneity (ITH) based on "SNP6" files, "maf" files, DNA methylation, mRNA expression, miRNA expression, lncRNA expression, or protein expression. It can also be used to calculate integrated scores based on the different types of ITH scores described above. 

## Details

- The  function `ieITH()` containing 3 parameters: data_tumor, data_normal and data_type.
    + "data_tumor" is a dataframe or matrix containing one of the multi-omics profiles in tumors: 1) For DNA methylation, mRNA expression, miRNA expression, lncRNA expression, or protein expressionin, the row name is probes or RNAs or proteins and the column name is tumor sample ID; 2) For "maf" files, it at least includes two columns ("sample" and "vaf"); 3) For "SNP6" files,it at least includes 2 columns ("sample" and "value").
    + "data_normal" is a dataframe or matrix containing matched normal sample profiles for DNA methylation, mRNA expression, miRNA expression, lncRNA expression, or protein expression, the row names of "data_normal" should be matched with the row names of "data_tumor". If the "data_normal" is NULL, calculate the ITH score by with tumor samples.
    + "data_type" is a character belonging to one of the c("cnv", "mut", "met", "lnc", "mir", "mrn", "pro").
- The input of the function `intITH()` should be a dataframe containing 3 columns ("sample", "ieITH_score" and "score_type").



## Installation

You can install the released version of **ieITH** with:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("WangX-Lab/ieITH/ieITH")
```

## Examples

```R
## Calculate the ITH score for each tumor sample based on profiles of mRNA expression.
path1 <- system.file("extdata", "example_mrnITH_tumor.txt", package = "ieITH", mustWork = TRUE)
path2 <- system.file("extdata", "example_mrnITH_normal.txt", package = "ieITH", mustWork = TRUE)
data_tumor <- read.table(path1, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
data_normal <- read.table(path2, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
mrnITH <- ieITH(data_tumor, data_normal, "mrn")

## Calculate the CNAs-based ITH score for each tumor sample with the input of "SNP6" files.
path <- system.file("extdata", "example_cnvITH.txt", package = "ieITH", mustWork = TRUE)
data_tumor <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
cnvITH <- ieITH(data_tumor, data_normal = NULL, "cnv")

## Calculate the somatic mutations-based ITH score for each tumor sample with the input of "maf" files.
path <- system.file("extdata", "example_mutITH.txt", package = "ieITH", mustWork = TRUE)
data_tumor <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
mutITH <- ieITH(data_tumor, data_normal = NULL, "mut")

## Integrate the different types of IE-based ITH scores. 
path <- system.file("extdata", "example_intITH.txt", package = "ieITH", mustWork = TRUE)
input_data_int <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
ieITH_score <- intITH(input_data_int)
```

### Contact

E-mail any questions to Hongjing Ai <hongjingai@stu.cpu.edu.cn>, Dandan Song <dandan.song@stu.cpu.edu.cn>, Xiaosheng Wang <xiaosheng.wang@cpu.edu.cn>

