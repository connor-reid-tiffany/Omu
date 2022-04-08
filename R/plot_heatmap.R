#' Create a heatmap
#' @description Takes a metabolomics count data frame and creates a heatmap. It is recommended to
#' either subset, truncate, or agglomerate by metabolite metadata to improve legibility.
#' @param count_data A metabolomics count data frame.
#' @param metadata The descriptive meta data for the samples.
#' @param Factor The column name for the independent variable in your metadata.
#' @param log_transform TRUE or FALSE. Recommended for visualization purposes. If true data is
#' transformed by the natural log.
#' @param high_color Color for high abundance values
#' @param low_color Color for low abundance values
#' @param aggregate_by Hierarchical metadata value to sum metabolite values by, i.e. "Class"
#' @param response_variable The response variable for the data, i.e. "Metabolite"
#' @importFrom utils stack
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#'
#' plot_heatmap(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' log_transform = TRUE, Factor = "Treatment", response_variable = "Metabolite",
#' aggregate_by = "Subclass_2", high_color = "darkgoldenrod1", low_color = "dodgerblue2")

#' @export

plot_heatmap <- function(count_data, metadata, Factor, response_variable,
                         log_transform=FALSE, high_color, low_color, aggregate_by){

  if(any(colnames(metadata)=="Sample")==FALSE){

          stop("metadata is missing Sample column")

                                                    }
  if(any(colnames(metadata)==Factor)==FALSE){

        stop("metadata is missing Factor column. Did you make a typo?")

      }

  count_data_sample <- count_data[,names(count_data) %in% metadata$Sample]

  if(identical(colnames(count_data_sample), as.character(metadata$Sample))==FALSE){

        stop("Sample names in count_data and metadata do not match.")

      }
  Abundance = NULL
  Sample = NULL
  #create dataframe to match hierarchical metadata to the new dataframe
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
  if (missing(Factor)){
    count_data <- count_data
  } else if (!missing(Factor)){

  #match desired factor data to new dataframe
  count_data[,Factor] <- metadata[,Factor][match(count_data[,"Sample"],
                                                 metadata[,"Sample"])]

  }
  if (missing(aggregate_by)){
    count_data = count_data
  }else if (!missing(aggregate_by)){

    if(any(names(hm_df) %in% aggregate_by)==FALSE){

      stop("Metabolomics data are missing metadata columns. Did you forget to use assign_hierarchy?")

    }

    count_data$Class <- hm_df$Class[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_1 <- hm_df$Subclass_1[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_2 <- hm_df$Subclass_2[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_3 <- hm_df$Subclass_3[match(count_data[,response_variable], hm_df[,response_variable])]
    count_data$Subclass_4 <- hm_df$Subclass_4[match(count_data[,response_variable], hm_df[,response_variable])]

  }

  #plot the heatmap
  if (!missing(aggregate_by) & !missing(Factor)){
  plot <- ggplot(data = count_data, mapping = aes(x = Sample, y =count_data[,aggregate_by], fill = Abundance)) +
    geom_tile() + xlab(label = "Sample") + scale_fill_gradient(name = "Abundance",
                                                 low = low_color,
                                                 high = high_color) +
    facet_grid(~ count_data[,Factor], switch = "x", scales = "free_x", space = "free_x") +
    theme_bw() + theme(panel.grid = element_blank())

  }else if (!missing(aggregate_by) & missing(Factor)){
    plot <- ggplot(data = count_data, mapping = aes(x = Sample, y = count_data[,aggregate_by], fill = Abundance)) +
      geom_tile() + xlab(label = "Sample") + scale_fill_gradient(name = "Abundance",
                                                                 low = low_color,
                                                                 high = high_color) +
      theme_bw() + theme(panel.grid = element_blank())

  }else if (missing(aggregate_by) & !missing(Factor)){
    plot <- ggplot(data = count_data, mapping = aes(x = Sample, y = count_data[,response_variable], fill = Abundance)) +
      geom_tile() + xlab(label = "Sample") + scale_fill_gradient(name = "Abundance",
                                                                 low = low_color,
                                                                 high = high_color) +
      facet_grid(~ count_data[,Factor], switch = "x", scales = "free_x", space = "free_x") +
      theme_bw() + theme(panel.grid = element_blank())

  }else if (missing(aggregate_by) $ missing(Factor)){
    plot <- ggplot(data = count_data, mapping = aes(x = Sample, y = count_data[,response_variable], fill = Abundance)) +
      geom_tile() + xlab(label = "Sample") + scale_fill_gradient(name = "Abundance",
                                                                 low = low_color,
                                                                 high = high_color) +
      theme_bw() + theme(panel.grid = element_blank())

  }

  return(plot)

}
