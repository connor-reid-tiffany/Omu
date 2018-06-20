#' plot_bar
#' Creates a ggplot2 object using the output file from the count_fold_changes function
#' @param data The output file from Count_Fold_Changes
#' @param fill A character vector of length 2 containing colors for filling the bars, Factors are alphanumeric, so the first color is for the "Decrease" bar while the second is for "Increase"
#' @param color A character vector of length 2 containing colors for the bar outlines
#' @param size A character vector of 2 numbers for the size of the bar outlines.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_size_manual
#' @importFrom ggplot2 scale_color_manual
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#' t_test_df <- omu_summary(data = c57_nos2KO_mouse_countDF, colData = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH")
#' fold_change_counts <- count_fold_changes(data = t_test_df, "Class",
#' column = "Class", sig_threshold = 0.05, keep_unknowns = FALSE)
#' plot_bar(data = fold_change_counts, fill = c("firebrick2", "dodgerblue2"),
#' color = c("black", "black"), size = c(1,1))
#' @export

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
