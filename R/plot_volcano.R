#' Create a volcano plot
#'
#' @description Creates a volcano plot as ggplot2 object using the output of omu_summary
#' @param count_data The output file from the omu_summary function.
#' @param column The column with metadata you want to highlight points in the plot with,
#' i.e. "Class"
#' @param strpattern A character vector of levels of the column you want the plot to focus on,
#' i.e. strpattern = c("Carbohydrates", "Organicacids")
#' @param fill A character vector of colors you want your points to be. Must be of length
#' 1 + length(strpattern) to account for points not in strpattern.
#' Levels of a factor are organzed alphabetically. All levels not in the
#' strpattern argument will be set to NA.
#' @param sig_threshold An integer. Creates a horizontal dashed line for a significance threshold.
#' i.e. sig_threshold = 0.05. Defaut value is 0.05
#' @param alpha A character vector for setting transparency of factor levels.Must be of length
#' 1 + length(strpattern) to account for points not in strpattern.
#' @param shape A character vector for setting the shapes for your column levels.
#' Must be of length 1 + length(strpattern) to account for points not in strpattern.
#' See ggplot2 for an index of shape values.
#' @param color A character vector of colors for the column levels.
#' Must be of length 1 + length(strpattern) to account for points not in strpattern.
#' If you choose to use shapes
#' with outlines, this list will set the outline colors.
#' @param size Size of the points in the plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom stringr str_replace_na
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#' \dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:20,]}
#' t_test_df <-  omu_summary(count_data = c57_nos2KO_mouse_countDF,
#' metadata = c57_nos2KO_mouse_metadata, numerator = "Strep", denominator = "Mock",
#' response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH", test_type = "welch")
#'
#' plot_volcano(count_data = t_test_df, column = "Class", strpattern = c("Carbohydrates"),
#' fill = c("firebrick2", "white"), sig_threshold = 0.05, alpha = c(1,1),
#' shape = c(1,24), color = c("black", "black"), size = 2)
#'
#' plot_volcano(count_data = t_test_df, sig_threshold = 0.05, size = 2)
#' @export




plot_volcano <- function(count_data, column, size, strpattern, fill, sig_threshold, alpha, shape,
  color){

 if(is.null(count_data$padj)==TRUE){

   stop("count_data must be the output of omu_summary or omu_anova and have a padj column
   and a log2FoldChange column")

 }



 log2FoldChange <- padj <- Metabolite <- NULL

  if (missing(sig_threshold)){
     sig_threshold = 0.05
  } else {sig_threshold = sig_threshold}
  if (missing(column)){
    count_data[, "Color"] <- NA
    count_data$Color = count_data$padj <= sig_threshold
      volcano_plot <- ggplot(count_data, aes(x = log2FoldChange,
        y = -log10(padj), text = paste("Metabolite:", Metabolite))) +
        geom_point(size = size, aes(color = count_data$Color, alpha = count_data$Color,
          shape = count_data$Color, fill = count_data$Color)) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "black")) +
        scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 1)) +
        scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21)) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white")) +
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed")

  }else{

    if(any(names(count_data) %in% column)==FALSE){

      stop("count_data is missing metabolite metadata. Did you forget to use assign_hierarchy?")

    }

    if(all(strpattern %in% count_data[,column])==FALSE){

      stop("One or more elements in strpattern are missing from column. Did you make a typo?")

    }

    if(isTRUE(length(fill)!=1 + length(strpattern))==FALSE){

      stop("the length of the fill argument must be 1 greater than the strpattern argument")

    }

    if(isTRUE(length(alpha)!=1 + length(strpattern))==FALSE){

      stop("the length of the alpha argument must be 1 greater than the strpattern argument")

    }

    if(isTRUE(length(shape)!=1 + length(strpattern))==FALSE){

      stop("the length of the shape argument must be 1 greater than the strpattern argument")

    }

    if(isTRUE(length(color)!=1 + length(strpattern))==FALSE){

      stop("the length of the color argument must be 1 greater than the strpattern argument")

    }
    data2 = count_data
    data2[,column] <- sapply(data2[,column], function(x) replace(x, x %in% strpattern, NA))
    data2[,column] <- factor(data2[,column])
    to_remove = levels(data2[,column])
    count_data[,column] = sapply(count_data[,column], function(x) replace(x, x %in% to_remove, NA))
    factor = count_data[,column]
    factor = str_replace_na(factor, replacement = "NA")
    count_data[,column] = factor(count_data[,column])
    count_data[sapply(count_data, is.character)] <- lapply(count_data[sapply(count_data, is.character)],
    as.factor)
    volcano_plot <- ggplot(count_data, aes(x = log2FoldChange, y = -log10(padj),
    text = paste("Metabolite:", Metabolite))) +
        geom_point(size = size, aes(fill = factor(factor), alpha = factor(factor),
        shape = factor(factor))) +
        scale_fill_manual(values = fill) +
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed") +
        scale_alpha_manual(values = alpha) +
        scale_shape_manual(values = shape) +
        scale_color_manual(values = color)

      }
return(volcano_plot)
  }
