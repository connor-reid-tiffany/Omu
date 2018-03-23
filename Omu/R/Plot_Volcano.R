#'plot_volcano
#'
#'This function takes a data frame from T_Test and creates a Volcano Plot as a ggplot2 graphical object.
#'@param data The output file from the T_Test().
#'@param factor A column of metadata from the data file, i.e. data$Class. If you provide the function with a factor parameter, every other parameter except sig_threshold is required. Choosing no factor parameter will produce a volcano plot with significant points in red, and not significant points in black.
#'@param column A character vector, i.e. column = c("Class")
#'@param strpattern A character vector of levels of the factor you want the plot to focus on, i.e. strpattern = c("Carbohydrates", "Organicacids")
#'@param fill A character vector of colors you want your points to be. Levels of a factor are organzed alphabetically. All levels not in the strpattern argument will be set to NA.
#'@param sig_threshold Numeric. Creates a horizontal dashed line for a significance threshold. i.e. sig_threshold = 0.05. Defaut value is 0.05
#'@param alpha A character vector for setting transparency of factor levels. i.e. alpha = (1, 0.5, 1)
#'@param shape A character vector for setting the shapes for your factor levels. See ggplot2 for an index of shapes
#'@param color A character vector of colors for the factor levels. If you choose to use shapes with outlines, this list will set the outline colors.
#'@example plot_volcano(data = data, column = "Subclass_2", strpattern = c("SugaralcoholsFig", "SugaracidsFig"),fill = c("black", "hotpink", "cyan"), sig_threshold = 1.301029996, alpha = c(0.25, 1, 1), shape = c(1, 24, 22),color = c("black", "black", "black"))
#'@export



plot_volcano <- function(data, column, strpattern, fill, sig_threshold,  alpha, shape, color){
  if (missing(sig_threshold)) sig_threshold = 0.05
  else sig_threshold = sig_threshold
  if (missing(strpattern)){
    data[, "Color"] <- NA
    data[,"Color"] = data[,"padj"] <= sig_threshold
      ggplot(data, aes(x = log2FoldChange, y = -log10(padj), text = paste("Metabolite:", Metabolite))) +
        geom_point(size = 3.5, aes(color = data[,"Color"], alpha = data[,"Color"], shape = data[,"Color"], fill = data[,"Color"])) +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "black")) +
        scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1)) +
        scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1)) +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "white"))+
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed")
  }else{
    data[, column] <- sapply(data[, column], as.character)
    Vector = data[,column]
    Match <- str_match(pattern = strpattern, string = Vector)
    Match = str_replace_na(Match, replacement = "NA")
    data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)],
                                               as.factor)
    ggplot(data, aes(x = log2FoldChange, y = -log10(padj), text = paste("Metabolite:", Metabolite))) +
        geom_point(size = 3.5, aes(fill = factor(Match), alpha = factor(Match), shape = factor(Match))) +
        scale_fill_manual(values = fill) +
        geom_hline(aes(yintercept = -log10(sig_threshold)), linetype = "dashed") +
        scale_alpha_manual(values = alpha) +
        scale_shape_manual(values = shape) +
        scale_color_manual(values = color)}
  }
