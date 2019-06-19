#' Create a box plot
#' @description Takes a metabolomics count data frame and creates boxplots. It is recommended to
#' either subset, truncate, or agglomerate by hierarchical metadata.
#' @param count_data A metabolomics count data frame, either from read_metabo or omu_summary
#' @param metadata The descriptive meta data for the samples
#' @param Factor The column name for the experimental variable
#' @param log_transform TRUE or FALSE. Recommended for visualization purposes. If true data is
#' transformed by the natural log
#' @param fill_list Colors for the plot which is colored by Factor, in the form of c("")
#' @param aggregate_by Hierarchical metadata value to sum metabolite values by, i.e. "Class"
#' @param response_variable The response variable for the data, i.e. "Metabolite"
#' @importFrom utils stack
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 scale_fill_manual
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#'
#' plot_boxplot(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' log_transform = TRUE, Factor = "Treatment", response_variable = "Metabolite",
#' aggregate_by = "Subclass_2", fill_list = c("darkgoldenrod1", "dodgerblue2"))

#' @export

plot_boxplot <- function(count_data, metadata, aggregate_by, log_transform,
                         Factor, response_variable, fill_list){
  Abundance = NULL
  Sample = NULL
  hm_df <- count_data
  #separate response variable from data frame and make a character vector for
  #subsetting the samples
  response_variable_vect <- as.vector(count_data[,response_variable])
  keep_columns <- as.character(metadata[, "Sample"])
  #replicate the response variable vector by the number of samples
  response_variable_df <- as.data.frame(rep(response_variable_vect,
                                            length(keep_columns)))
  #rename the column in the new replicated dataframe to response variable
  colnames(response_variable_df)[1] <- response_variable
  #subset out the sample columns only
  count_data <- count_data[,keep_columns]
  #log_transform
  if (log_transform==TRUE){
    count_data = log(count_data)
  } else if (log_transform==FALSE){
    count_data = count_data
  }
  #use stack to transform into a sample column and abundance column
  count_data <- stack(count_data)
  #rename columns to Sample and Abundance
  colnames(count_data)[1] <- "Abundance"
  colnames(count_data)[2] <- "Sample"
  #bind response variable dataframe to new dataframe
  count_data <- cbind(response_variable_df, count_data)
  #match desired factor data to new dataframe
  count_data[,Factor] <- metadata[,Factor][match(count_data[,"Sample"],
                                                   metadata[,"Sample"])]
  colnames(count_data)[4] <- "Factor"
  if (missing(aggregate_by)){
    count_data = count_data
  }else if (!missing(aggregate_by)){
    count_data$Class <- hm_df$Class[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_1 <- hm_df$Subclass_1[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_2 <- hm_df$Subclass_2[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_3 <- hm_df$Subclass_3[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_4 <- hm_df$Subclass_4[match(count_data[,response_variable], hm_df[,response_variable])]

  }

  #plot_box and whiskers plot

  if (!missing(aggregate_by)){
    plot <- ggplot(count_data, aes(x=Factor, y=Abundance, fill = Factor)) +
      geom_boxplot()+
      facet_wrap(~ count_data[,aggregate_by], strip.position = "bottom", scales = "free_y") +
      scale_fill_manual(values = fill_list) + theme_bw()
  }else if (missing(aggregate_by)){
    plot <- ggplot(count_data, aes(x=Factor, y=Abundance, fill = Factor)) +
      geom_boxplot()+
      facet_wrap(~ count_data[,response_variable], strip.position  = "bottom", scales = "free_y")+
      scale_fill_manual(values = fill_list) + theme_bw()
  }

  return(plot)
}
