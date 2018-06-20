#' plot_volcano
#'
#' This function takes an input data frame from t_test and creates a volcano plot as a ggplot2 object.
#' @param data The output file from the Clean_DESeq_results() function.
#' @param column The factor as a string, i.e. column = "Class"
#' @param strpattern A character vector of levels of the factor you want the plot to focus on, i.e. strpattern = c("Carbohydrates", "Organicacids")
#' @param fill A character vector of colors you want your points to be. Levels of a factor are organzed alphabetically. All levels not in the strpattern argument will be set to NA.
#' @param sig_threshold An integer. Creates a horizontal dashed line for a significance threshold. i.e. sig_threshold = 0.05. Defaut value is 0.05
#' @param alpha A character vector for setting transparency of factor levels. i.e. alpha_list = (1, 0.5, 1)
#' @param shape A character vector for setting the shapes for your factor levels. See ggplot2 for an index of shapes
#' @param color A character vector of colors for the factor levels. If you choose to use shapes with outlines, this list will set the outline colors.
#' @param size Size of points in plot
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
#'
#' t_test_df <-  omu_summary(data = c57_nos2KO_mouse_countDF, colData = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH")
#'
#' plot_volcano(data = t_test_df, column = "Class", strpattern = c("Carbohydrates", "Lipids"),
#' fill = c("firebrick2", "dodgerblue2", "white"), sig_threshold = 0.05, alpha = c(1,1,1),
#' shape = c(1,21,24), color = c("black", "black", "black"), size = 2)
#'
#' plot_volcano(data = t_test_df, sig_threshold = 0.05, size = 2)
#' @export




plot_volcano <- function(data, column, size, strpattern, fill, sig_threshold, alpha, shape, color){
  if (missing(sig_threshold)){
     sig_threshold = 0.05
  } else {sig_threshold = sig_threshold}
  if (missing(column)){
    data[, "Color"] <- NA
    data$Color = data$padj <= sig_threshold
      ggplot(data, aes(x = log2FoldChange, y = -log10(padj), text = paste("Metabolite:", Metabolite))) +
        geom_point(size = size, aes(color = data$Color, alpha = data$Color, shape = data$Color, fill = data$Color)) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "black")) +
        scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 1)) +
        scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 21)) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white"))+
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed")
  }else{
    data2 = data
    data2[,column] <- sapply(data2[,column], function(x) replace(x, x %in% strpattern, NA))
    data2[,column] <- factor(data2[,column])
    to_remove = levels(data2[,column])
    data[,column] = sapply(data[,column], function(x) replace(x, x %in% to_remove, NA))
    factor = data[,column]
    factor = str_replace_na(factor, replacement = "NA")
    data[,column] = factor(data[,column])
    data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)],
                                               as.factor)
    ggplot(data, aes(x = log2FoldChange, y = -log10(padj), text = paste("Metabolite:", Metabolite))) +
        geom_point(size = size, aes(fill = factor(factor), alpha = factor(factor), shape = factor(factor))) +
        scale_fill_manual(values = fill) +
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed") +
        scale_alpha_manual(values = alpha) +
        scale_shape_manual(values = shape) +
        scale_color_manual(values = color)}
  }
