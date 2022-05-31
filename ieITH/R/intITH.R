#' @title  intITH
#'
#' @description
#' \code{intITH} intITH is the average of different IE-based ITH scores after scaling into [0-1] within each type of ITH scores. 
#'
#' @param input_data_int  A dataframe containing 3 columns: "sample", "ieITH_score" and "score_type".
#'
#' @return A dataframe with 2 columns:
#' \item{sample}{Tumor samples to be calculated.}
#' \item{intITH_score}{The intITH score of each sample.}
#' @export
#' @importFrom reshape2 dcast
#'
#' @examples
#' path <- system.file("extdata", "example_intITH.txt", package = "ieITH", mustWork = TRUE)
#' input_data_int <- read.table(path, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, sep = "\t", quote = "")
#' ieITH_score <- intITH(input_data_int)


intITH <- function (input_data_int) {

	# Check arguments ------------------------------------------------------------
	if (missing(input_data_int) || class(input_data_int) != "data.frame")
  stop("'input_data_int' is missing or incorrect")

  pos <- match("ieITH_score", colnames(input_data_int))
  colnames(input_data_int)[pos] <- "value"
  wide_data <- reshape2::dcast(input_data_int,score_type~sample)
  rownames(wide_data) <- wide_data[,1]
  wide_data <- wide_data[,-1]

  # scale ITH scores to [0, 1]
	maxRow <- apply(wide_data, 1, max, na.rm=T)
  minRow <- apply(wide_data, 1, min, na.rm=T)	
  scale_data <- (wide_data - minRow) / (maxRow - minRow)
  	 
  scale_data <- as.data.frame(t(scale_data))
  mean_score <- rowMeans(scale_data, na.rm = TRUE) # calculate the mean score
  intITH_score <- data.frame(sample = rownames(scale_data), intITH_score = as.numeric(mean_score))
  return(intITH_score)
}