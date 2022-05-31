#' @title  ieITH
#'
#' @description
#' \code{ieITH} ieITH evaluates the ITH level of each bulk tumor sample based on one of the multi-omics profiles in tumors. 
#' 
#' @details The multi-omics profiles could be "SNP6" files, "maf" files, DNA methylation, mRNA expression, miRNA expression, lncRNA expression, and protein expression.
#'
#' @param data_tumor A dataframe or matrix containing one of the multi-omics profiles in tumors: 
#' 1) For DNA methylation, mRNA expression, miRNA expression, lncRNA expression, and protein expressionin, the row name is probes or RNAs or proteins and the column name is tumor sample ID;
#' 2) For "maf" files, it at least includes two columns ("sample" and "vaf");
#' 3) For "SNP6" files,it at least includes two columns ("sample" and "value").
#' @param data_normal A dataframe or matrix containing matched normal sample profiles for DNA methylation, mRNA expression, miRNA expression, lncRNA expression, and protein expressionin, 
#' the row names of data_normal should be matched with the row names of data_tumor.If the data_normal is NULL, calculate the IE-based ITH score by with tumor samples.
#' 
#' @param data_type A character belonging to one of the c("cnv", "mut", "met", "lnc", "mir", "mrn", "pro").

#' @return A dataframe with 3 columns:
#' \item{sample}{Tumor samples to be calculated.}
#' \item{ieITH_score}{The ieITH score of each sample.}
#' \item{score_type}{The type of IE-based ITH score.}
#' @export
#'
#' @examples
#' #example1
#' path1 <- system.file("extdata", "example_mrnITH_tumor.txt", package = "ieITH", mustWork = TRUE)
#' path2 <- system.file("extdata", "example_mrnITH_normal.txt", package = "ieITH", mustWork = TRUE)
#' data_tumor <- read.table(path1, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
#' data_normal <- read.table(path2, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "", row.names = 1)
#' mrnITH <- ieITH(data_tumor, data_normal, "mrn")
#' #example2
#' path <- system.file("extdata", "example_cnvITH.txt", package = "ieITH", mustWork = TRUE)
#' data_tumor <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
#' cnvITH <- ieITH(data_tumor, data_normal = NULL, "cnv")
#' #example3
#' path <- system.file("extdata", "example_mutITH.txt", package = "ieITH", mustWork = TRUE)
#' data_tumor <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
#' mutITH <- ieITH(data_tumor, data_normal = NULL, "mut")

ieITH <- function (data_tumor, data_normal = NULL, data_type) {

  # Check arguments -----------------------------------------------------------------------------
  if (missing(data_type) || !(data_type %in% c("cnv", "mut", "met", "lnc", "mir", "mrn", "pro")))
    stop("'data_type' is missing or incorrect")
  if (is.null(data_normal) && data_type %in% c("met", "lnc", "mir", "mrn", "pro")) {
    message("'data_normal' is missing, use tumor samples as control")
    data_normal <- data_tumor
  }
  if (missing(data_tumor) || !class(data_tumor) %in% c("matrix", "data.frame"))
    stop("'data_tumor' is missing or incorrect")
  if (!class(data_normal) %in% c("matrix", "data.frame") && !(data_type %in% c("cnv", "mut")))
    stop("'data_normal' is missing or incorrect")
  
  # Calculate the ITH score for each tumor sample based on profiles of DNA methylation, lncRNA, miRNA, mRNA, or protein expression.
  if (data_type %in% c("met","lnc","mir","mrn","pro")) {
    normal_mean = apply(data_normal, 1, mean, na.rm = T)
    tumor_abs = abs(data_tumor - normal_mean)

    maxRow <- apply(tumor_abs, 1, max, na.rm=T)
    minRow <- apply(tumor_abs, 1, min, na.rm=T)
    
    tumor_scale <- (tumor_abs - minRow) / (maxRow - minRow)

    group <- 20
    gap <- 1 / group

    score <- c()
    for(j in 1 : ncol(tumor_scale)) {
      col_temp <- tumor_scale[, j] / gap
      col_group <- floor(col_temp)
      col_group[col_group == group] <- group - 1
      p_i <- table(col_group) / sum(table(col_group))
      entropy <- sum(-p_i * log2(p_i))
      score <- c(score, entropy)
    }

    ieITH_score <- data.frame(sample = colnames(data_tumor), ieITH_score = score, score_type = paste0(data_type, "ITH"))
    return(ieITH_score)
  }
  
  # Calculate the CNAs-based ITH score for each tumor sample with the input of "SNP6" files.
  if (data_type == "cnv") {
    group <- 20  
    gap <- 6 / group

    sample_cnv <- unique(data_tumor[, "sample"])

    score <- c()
    for (j in 1 : length(sample_cnv)) {
      sample_cnv0 <- data_tumor[data_tumor[ ,"sample"] == sample_cnv[j], ]
      sample_cnv0 <- sample_cnv0[sample_cnv0[ ,"value"] <= 3 & sample_cnv0[ ,"value"] >= -3, ]
      col_temp <- (sample_cnv0[ , "value"] + 3) / gap
      col_group <- floor(col_temp)
      col_group[col_group == group] <- group - 1
      p_i <- table(col_group) / sum(table(col_group))  
      entropy <- sum(-p_i * log2(p_i))
      score <- c(score, entropy)
    }

    ieITH_score <- data.frame(sample = sample_cnv, ieITH_score = score, score_type = paste0(data_type, "ITH"))
    return(ieITH_score)
  }

  # Calculate the somatic mutations-based ITH score for each tumor sample with the input of "maf" files.
  if (data_type == "mut") {
    group <- 20  
    gap <- 1 / group

    sample_mut <- unique(data_tumor[, "sample"])

    score <- c()
    for (j in 1:length(sample_mut)) {
      sample_mut0 <- data_tumor[data_tumor[ , "sample"] == sample_mut[j], ]
      col_temp <- sample_mut0[, "vaf"] / gap
      col_group <- floor(col_temp)
      col_group[col_group == group] <- group - 1
      p_i <- table(col_group) / sum(table(col_group))  
      entropy <- sum(-p_i * log2(p_i))
      score <- c(score, entropy)
    } 

    ieITH_score <- data.frame(sample = sample_mut, ieITH_score = score, score_type = paste0(data_type, "ITH"))
    return(ieITH_score)
  }
}
