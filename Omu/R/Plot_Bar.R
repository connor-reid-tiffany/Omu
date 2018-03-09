#' plot_bar
#' Creates a ggplot2 object using the output file from the Count_Fold_Changes function
#'@param data The output file from Count_Fold_Changes
#'@param fill A list of length 2 containing colors for filling the bars, Factors are alphanumeric, so the first color is for the "Decrease" bar while the second is for "Increase"
#'@param color A list of length 2 containing colors for the bar outlines
#'@param size A list of 2 numbers for the size of the bar outlines.
#'@keywords metabo
#'@export
#'@examples plot_bar(data = data, fill = levels, color = levels, size = levels)

plot_bar <- function(data, fill, size, color){
colnames(data)[1] <- "Class"
    ggplot(data = data, aes(x=reorder(Class, -Significant_Changes), y = Significant_Changes)) +
           geom_bar(stat = "identity", aes(fill = colour, color = colour, size = colour)) +
           theme_bw() +
           theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size = 10), axis.text.y = element_text(size = 10)) +
           theme(panel.border = element_blank(), axis.line = element_line(color = "Black", size=1, lineend = "square")) +
           theme(plot.title = element_text(hjust = 0.5)) +
           scale_size_manual(values = size) +
           scale_fill_manual(values = fill) +
           scale_color_manual(values = color)
}
