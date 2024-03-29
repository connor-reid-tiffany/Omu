#' Create a PCA plot

#' @description Performs an ordination and outputs a PCA plot using a metabolomics
#' count data frame and metabolomics metadata
#' @param count_data Metabolomics count data
#' @param metadata Metabolomics metadata
#' @param variable The independent variable you wish to compare and contrast
#' @param color String of what you want to color by. Usually should be the same as variable.
#' @param response_variable String of the response_variable, usually should be "Metabolite"
#' @param label TRUE or FALSE, whether to add point labels or not
#' @param size An integer for point size.
#' @param ellipse TRUE or FALSE, whether to add confidence interval ellipses or not.
#' @import ggfortify
#' @importFrom ggplot2 autoplot
#' @importFrom stats prcomp
#' @examples
#' PCA_plot(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' variable = "Treatment", color = "Treatment", response_variable = "Metabolite")
#' @export

PCA_plot <- function(count_data, metadata, variable, color, response_variable="Metabolite", label = FALSE, size = 2, ellipse = FALSE){

  if(any(names(metadata) %in% variable)==FALSE){

    "variable not found in metadata. did you make a typo?"

  }

  if(any(names(metadata) %in% color)==FALSE){

    "color variable not found in metadata. did you make a typo?"

  }

  if(any(names(count_data) %in% response_variable)==FALSE){

    stop("metabolomics data are missing the response variable column. Did you make a typo?")

  }

  if(identical(sort(as.character(colnames(count_data)[unlist(lapply(count_data, is.numeric))])), sort(as.character(metadata$Sample))==FALSE)){

    stop("Sample names in count_data and metadata do not match.")

  }

  if(any(colnames(metadata)=="Sample")==FALSE){

    stop("metadata is missing Sample column")

  }



  rownames(count_data) <- count_data[,response_variable]
  count_data[,response_variable] <- NULL

  data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]


  #Transform for 'normalization' and T test
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
  data_Transpose <- data_Transpose[order(data_Transpose$Sample),]
  metadata <- metadata[order(metadata$Sample),]
  data_Transpose[, variable] = metadata[, variable][match(metadata$Sample,
                                                          data_Transpose$Sample)]
  data_Numeric <- data_Transpose[, !names(data_Transpose) %in% c(variable, 'Sample')]
  data_Numeric <- data.frame(lapply(data_Numeric, function(x) as.numeric(as.character(x))),check.names=F, row.names = rownames(data_Numeric))
  if(ellipse==FALSE && label==FALSE){

  Plot <- autoplot(prcomp(data_Numeric), data = data_Transpose, colour = color, size = size)

}else if(ellipse==TRUE && label==TRUE){

    Plot <- autoplot(prcomp(data_Numeric), data = data_Transpose, colour = color, frame = TRUE,label = label, label.size = 3.5, label.repel = T, frame.type = 'norm', size = size)

  }else if(ellipse==TRUE && label==FALSE){

      Plot <- autoplot(prcomp(data_Numeric), data = data_Transpose, colour = color, frame = TRUE, frame.type = 'norm', size = size)

    }else if(ellipse==FALSE && label==TRUE){

        Plot <- autoplot(prcomp(data_Numeric), data = data_Transpose, colour = color,label = label, label.size = 3.5, label.repel = T, size = size)

      }

  return(Plot)
}
