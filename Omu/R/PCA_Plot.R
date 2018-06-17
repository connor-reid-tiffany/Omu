#' PCA_plot

#' Performs an ordination and outputs a PCA plot
#' @param data Metabolomics count data
#' @param colData Metabolomics meta data
#' @param variable The independent variable you wish to compare and contrast
#' @param color String of what you want to color by. usually should be the same as variable
#' @param response_variable String of the response_variable, usually should be "Metabolite"
#' @import ggfortify
#' @importFrom ggplot2 autoplot
#' @importFrom stats prcomp
#' @examples
#' PCA_plot(data = c57_nos2KO_mouse_countDF, colData = c57_nos2KO_mouse_metadata,
#' variable = "Treatment", color = "Treatment", response_variable = "Metabolite")
#' @export

PCA_plot <- function(data, colData, variable, color, response_variable){

  rownames(data) <- data[,response_variable]
  data[,response_variable] <- NULL

  data_Int <- data[sapply(data, function(x) is.numeric(x))]


  #Transform for 'normalization' and T test
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
  data_Transpose[, variable] = colData[, variable][match(colData$Sample, data_Transpose$Sample)]
  data_Numeric <- data_Transpose[, !names(data_Transpose) %in% c(variable, 'Sample')]
  data_Numeric <- data.frame(lapply(data_Numeric, function(x) as.numeric(as.character(x))),check.names=F, row.names = rownames(data_Numeric))

  Plot <- autoplot(prcomp(data_Numeric), data = data_Transpose, colour = color, frame = TRUE, frame.type = 'norm')


  return(Plot)
}
